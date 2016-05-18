/* libCSAM
 	Copyright (C)2013-2016 Rodrigo Canovas

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
#include "./CRPSPreMF.h"

CRPSPreMF::CRPSPreMF(cds_word sample,  bool create){
	buffRPS = NULL;
	buffQual = NULL;
	buffField = NULL;
	lastFieldBlock = (cds_word)-1;
	initData(sample, create);
}

CRPSPreMF::~CRPSPreMF(){
	if(huffA != NULL){  //if one exist, all off them exist
		delete huffA; delete huffC; delete huffG; delete huffT;
	} 
	if(buffQual != NULL)
		delete [] buffQual;
	if(buffField != NULL)
		delete [] buffField;
	if(buffRPS != NULL)
		delete [] buffRPS;
}

void CRPSPreMF::CreateFile(ifstream &SamFile, ofstream &fp, int SamOrRps, vector<cds_word>* IndexClusterToLines){
	cds_word *index;
	size_t index_pos = 0;
	getValuesPreSeq(posWin, WINSIZE); // get for last time the values from winSeq
	for(cds_word j = 0; j < WINSIZE; j++) //free winSeq
		delete [] WinSeq[j];
	delete [] WinSeq;
	if(tmp_file_string.length() > 0){
		fileSeq << tmp_file_string;  //write to fileQual (it must exist)
		tmp_file_string = "";
	}
	ClusterIndex[ClusterIndex.size() - 1] = ClusterReads;
	if(overlapcluster){ //add extra block if there is at least one read overlap with the next block
		ClusterIndex.push_back((cds_word)-1);
		if(IndexClusterToLines != NULL)
			IndexClusterToLines->push_back(numberOfLines);
		ClusterOverlap.push_back(overlapcluster);
	}
	fileSeq.close();
	huffA = new Huffman(occ_notA, 249, true);
	huffC = new Huffman(occ_notC, 249, true);
	huffG = new Huffman(occ_notG, 249, true);
	huffT = new Huffman(occ_notT, 249, true);
	SaveValue(fp, ROW);
	cout << "ClusterSize: " << ClusterOverlap.size() << endl;
	cds_word novacio = 0;
	for(cds_word i = 0; i < ClusterIndex.size()  ; i++){
		if(ClusterIndex[i] != (cds_word)-1)
			novacio++;
	}
	cout << "No vacio: " << novacio << endl;		
	index = &ClusterIndex[0]; //Index data
	if(ClusterIndex.size() != ClusterOverlap.size())
		cout << "Error: Number of blocks created is not correct" << endl;
	SaveValue(fp, (cds_word)ClusterIndex.size());
	index_pos = fp.tellp();
	SaveValue(fp, index, ClusterIndex.size());
	SaveValue(fp, variableSize);
	//if(variableSize){
		index = &ClusterOverlap[0];
		SaveValue(fp, index, ClusterOverlap.size());
	//}	
	cds_word Variables[6];
	Variables[0] = numberOfLinesNR;	Variables[1] = numberOfLines;	Variables[2] = sizeLine.size();
	Variables[3] = preSeqLength;	Variables[4] = (cds_word)(RNames.size()); Variables[5] = IndexRate;
	SaveValue(fp, Variables, 6);
	SaveValue(fp, &sizeLine[0], sizeLine.size());
	for(cds_word i = 0; i < RNames.size(); i++){
		SaveValue(fp, (cds_word)(RNames[i].length()));
		SaveValue(fp, (char *)(RNames[i].c_str()), RNames[i].length());
		SaveValue(fp, RNamesCluster[i]);
	}
	huffA->Save(fp); huffC->Save(fp); huffG->Save(fp); huffT->Save(fp);
	init_preseq = fp.tellp();
	init_preseq += 3 * sizeof(size_t);
	SaveValue(fp, init_preseq);
	SaveValue(fp, init_seq);
	SaveValue(fp, end_seq);
	buffer = init_buffer();
	buffer_use = 0;
	CompressPressumeSeq(fp, &buffer, &buffer_use);
	init_seq = fp.tellp();
	//RE-READ SAM FILE AND CREATE NEW FILE
	ComputeCSeq(SamFile, fp, SamOrRps);
	end_seq = fp.tellp();
	fp.seekp(index_pos);
	index = &ClusterIndex[0];
	SaveValue(fp, index, ClusterIndex.size());
	fp.seekp(init_preseq - 2 * sizeof(size_t));
	SaveValue(fp, init_seq);
	SaveValue(fp, end_seq);
	remove(tmp_name.c_str());
	fp.seekp(end_seq);
}

void CRPSPreMF::ComputeCSeq(ifstream &SamFile, ofstream &fp, int SamOrRps){
	string line, Ref, Seq, PSeq;
	ifstream fileTemp;
	vector<string> tokens;
	cds_word location = 0;
	size_t pos_input = 0, end_input = 0;
	cds_word posPre = 0,  length_block = 0, Pos = 0, totalPos = 0;
	cds_word cluster_ref = 0, last_cluster = (cds_word)-1, next = 0;
	bitsLength = (cds_word)(ceil(log2(sizeLine.size()))); //bits used to store the line lengths
	buffer_use = last_pos = farther_pos = 0;
	last_ref = "*";
	tmp_file_string = "";
	fileTemp.open(tmp_name.c_str()); //OPEN PRESUMED SEQUENCE
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	getline(SamFile, line);
	while (SamFile.good()){
		if(line[0] != '@'){
			Tokenize(line, tokens, "\t");
			if(SamOrRps){
				Pos = atoi(tokens[1].c_str());
				Ref = tokens[0];
				Seq = tokens[2];
			}
			else{
				Pos = atoi(tokens[3].c_str());      
				Ref = tokens[2];      
				Seq = tokens[9];
			}
			if((tmp_file_string.length() - posPre) < (farther_pos - last_pos + Seq.length())){
				pos_input += posPre;
				fileTemp.seekg(pos_input);
				length_block = (end_input - pos_input);
				if(length_block > MAXRAMBYTES)
					length_block = MAXRAMBYTES;
				tmp_file_string = "";
				tmp_file_string = LoadValue<char>(fileTemp, length_block);
				posPre = 0;
			}
			if(Ref.compare("*")){ //if ref is not '*'
				if(!last_ref.compare(Ref)){
					if(farther_pos  > Pos){
						posPre += (Pos - last_pos);
						totalPos += (Pos - last_pos);
					}
					else{
						posPre += farther_pos - last_pos; 
						totalPos += farther_pos - last_pos; 
					}
				}
				else{ //NEW REFERENCE
					posPre += farther_pos - last_pos;
					totalPos += farther_pos - last_pos;
					cluster_ref = getFirstCluster(Ref, &next);
				}
				length_block = cluster_ref + Pos / IndexRate; //compute the Index position
				if(length_block != last_cluster){ //add block info and pointer to the index
					ClusterReads = ClusterIndex[length_block];
					ClusterIndex[length_block] = location;
					last_cluster = length_block;
					SetVarField(buffer, buffer_use,  buffer_use + kWordSize - 1, ClusterReads); //Save number of read in the block
					buffer_use += kWordSize;
					SetVarField(buffer, buffer_use,  buffer_use + kWordSize - 1, totalPos); //Save pointer to pressume seq
					buffer_use += kWordSize;
					location += 2 * kWordSize;
				}
			}
			else //NOT REFERENCE READ
				cout << "Error: Not referenced read are ignored so far" << endl;
			PSeq = tmp_file_string.substr(posPre, Seq.size());
			location += CompressLine(Ref, Pos, Seq, PSeq, &buffer, &buffer_use);
			check_buffer(buffer, &buffer_use, fp);
		}
		tokens.clear();
		getline(SamFile,line);
	}	
	if(buffer_use != 0)
		SaveValue(fp, buffer,  (buffer_use + kWordSize - 1) / kWordSize);
	delete [] buffer;
	buffer = NULL;
}


void CRPSPreMF::DecompressSequence(ifstream &input, ofstream &output){
	cds_word length_block = 0, last_read_buffer = 0, pos_buffer = 0, actualPos = 0;
	cds_word countBlocks = 0, reads_to_copy = 0, posPre = 0, rel_pos = 0;
	string preseq, reads;
	bool first = false;
	cds_word last_length = 0;
	size_t pos_input = init_seq;
	input.seekg(pos_input);
	if(buffer != NULL)
		delete [] buffer;
	tmp_file_string = "";
	length_block = (end_seq - pos_input);  //get block of data
	if(length_block > MAXRAMBYTES)                
		length_block = MAXRAMBYTES;           
	buffer = (cds_word *)(LoadValue<char>(input, length_block));
	pos_buffer = 0;
	while(countBlocks < ClusterIndex.size()){
		if(ClusterIndex[countBlocks] != (cds_word)-1){
			ClusterReads = GetVarField(buffer, pos_buffer,  pos_buffer + kWordSize - 1);
			pos_buffer += kWordSize;
			posPre = GetVarField(buffer, pos_buffer,  pos_buffer + kWordSize - 1);
			pos_buffer += kWordSize;
			first = true;
			if(variableSize)      
				last_length = ClusterOverlap[countBlocks];
			else
				last_length = 0;
			preseq = GetPressumeSeq(input, posPre, posPre + IndexRate + 350); //get pressume sequence of the block, assuming 350 maximum read length
			actualPos = 0;
			reads_to_copy = ClusterReads;
			while(reads_to_copy > 0){
				reads = DecompressLine(buffer, &pos_buffer, preseq, &actualPos, &first, &last_length, &rel_pos) + "\n";
				tmp_file_string += reads;
				check_tmp_string(&tmp_file_string, output);
				reads_to_copy --;
				if((8 * length_block - pos_buffer) < (kWordSize * 350) && (pos_input + length_block) != end_seq){	//check if other block of data is neccesary
					pos_input += (pos_buffer / 8);
					input.seekg(pos_input);
					last_read_buffer = pos_buffer - 8 * (pos_buffer / 8);
					length_block = (end_seq - pos_input);  //get block of data
					if(length_block > MAXRAMBYTES)                      
						length_block = MAXRAMBYTES;    
					buffer = (cds_word *)(LoadValue<char>(input, length_block)); 
					pos_buffer = last_read_buffer;
				}	
			}
		}
		countBlocks ++;
	}
	delete [] buffer;
	buffer = NULL;
	if(tmp_file_string.length() > 0){
		output << tmp_file_string;            
		tmp_file_string = "";           
	}
}


cds_word CRPSPreMF::GetInterval(ifstream &input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, size_t MaxBytes){
	cds_word total_reads = 0;
	cds_word cluster_ref = (cds_word)-1, next = 0; 
	cds_word POS = 0,  rel_pos = 0, line_count = 0; 
	cds_word cluster_x, cluster_y, max_cluster, Ini_cluster;
	cds_word pos_buffRPS = 0, reads_in_cluster = 0, posPre = 0;
	cds_word last_length = 0, actualPos = 0;
	string datablock = "", preseq = "", seq = "";
	size_t pos_input = 0, sizeAux = 0;
	bool first = false; 
	cluster_ref = getFirstCluster(rname, &next); //look if rname is valid
	if(cluster_ref == (cds_word)-1)
		return total_reads;    //empty interval, rname not valid
	if(next == RNamesCluster.size())
		max_cluster = ClusterIndex.size();
	else
		max_cluster =  RNamesCluster[next] - 1;
	cluster_x = cluster_ref + pos_x / IndexRate;
	Ini_cluster =  IndexRate * (cds_word) (pos_x / IndexRate);
	cluster_y = cluster_ref + pos_y / IndexRate;
	if(cluster_x > max_cluster)
		return total_reads;    //empty interval, interval not valid
	if(cluster_y >= max_cluster)
		cluster_y = max_cluster - 1;
	while(cluster_x  <= cluster_y){
		POS = Ini_cluster;
		if(POS > pos_y)
			break;
		if(ClusterIndex[cluster_x] != (cds_word)-1){
			pos_input = init_seq + (ClusterIndex[cluster_x] / 8);
			if((buffRPS == NULL) || (pos_input < init_buffRPS) || ((init_buffRPS + MaxBytes) <= (pos_input + kBytesPerWord * 350))){
				input.seekg(pos_input);
				sizeAux = (end_seq - pos_input); //get block of data
				if(sizeAux > MaxBytes)
					sizeAux = MaxBytes; 
				if(buffRPS != NULL)
					delete [] buffRPS;
				buffRPS = (cds_word *)(LoadValue<char>(input, sizeAux));
				pos_buffRPS = ClusterIndex[cluster_x] - 8 * (ClusterIndex[cluster_x] / 8);
				init_buffRPS = pos_input;
			}
			else{
				sizeAux = MaxBytes;
				if((init_buffRPS + sizeAux) > end_seq)
					sizeAux = end_seq - init_buffRPS;
				pos_buffRPS = ClusterIndex[cluster_x] - 8 * (init_buffRPS - init_seq);
			}
			ClusterReads = GetVarField(buffRPS, pos_buffRPS,  pos_buffRPS + kWordSize - 1);
			pos_buffRPS += kWordSize;
			posPre = GetVarField(buffRPS, pos_buffRPS,  pos_buffRPS + kWordSize - 1);
			pos_buffRPS += kWordSize;
			first = true;
			if(variableSize)
				last_length = ClusterOverlap[cluster_x];
			else
				last_length = 0;
			preseq = GetPressumeSeq(input, posPre, posPre + IndexRate + 350); //get pressume sequence of the block, assuming 350 maximum read length
			actualPos = 0;
			reads_in_cluster = ClusterReads;			
			while(reads_in_cluster > 0){				
				seq = DecompressLine(buffRPS, &pos_buffRPS, preseq, &actualPos, &first, &last_length, &rel_pos);
				POS += rel_pos;
				if(POS > pos_y)
					break;
				if(POS >= pos_x){
					datablock += rname + "\t" +  NumtoString(POS) + "\t" + seq + "\n";
					line_count ++;
					total_reads ++;
					if((line_count % IndexRate) == 0){
						check_tmp_string(&datablock, output);
						if(datablock == "")
							line_count = 0;
					}
				}
				reads_in_cluster --;
				if(((8 * sizeAux - pos_buffRPS) < (kWordSize * 350)) && ((init_buffRPS + sizeAux) != end_seq)){  //check if other block of data is neccesary
					pos_input = init_buffRPS + pos_buffRPS / 8;
					sizeAux = (end_seq - pos_input); //get block of data
					if(sizeAux > MaxBytes)
						sizeAux = MaxBytes;              
					if(buffRPS != NULL)
						delete [] buffRPS;
					input.seekg(pos_input);
					buffRPS = (cds_word *)(LoadValue<char>(input, sizeAux));
					init_buffRPS = pos_input;
					pos_buffRPS = pos_buffRPS - 8 * (pos_buffRPS / 8);
				}
			}
		}
		cluster_x ++;
		Ini_cluster += IndexRate;
	}	
	if(datablock.length() != 0)
		output << datablock;
	return total_reads;
}


CRPSPreMF *CRPSPreMF::Load(ifstream &fp){
	cds_word r = LoadValue<cds_word>(fp);
	if (r != ROW) {
		assert(false);
		return NULL; 
	}
	cds_word sizeblock = 0, rnames_size = 0, aux_length = 0;
	cds_word size_readlength = 0;
	string ref = "";
	char *aux_text;
	string s;
	cds_word *index, *readlengths;
	CRPSPreMF *crps = new CRPSPreMF(0, false);
	sizeblock = LoadValue<cds_word>(fp);
	index = LoadValue<cds_word>(fp, sizeblock);
	crps->ClusterIndex.assign(index, index + sizeblock);
	delete [] index; 
	crps->variableSize = LoadValue<int>(fp);
	//if(crps->variableSize){
		index = LoadValue<cds_word>(fp, sizeblock);
		crps->ClusterOverlap.assign(index, index + sizeblock);
		delete [] index;
	//}
	cds_word *Variables = LoadValue<cds_word>(fp, 6);	
	crps->numberOfLinesNR = Variables[0];
	crps->numberOfLines = Variables[1];
	size_readlength = Variables[2];
	crps->preSeqLength = Variables[3];
	rnames_size = Variables[4];
	crps->IndexRate = Variables[5];
	readlengths = LoadValue<cds_word>(fp, size_readlength);
	crps->sizeLine.assign(readlengths, readlengths + size_readlength);
	crps->bitsLength =  (cds_word)(ceil(log2(size_readlength)));
	for(cds_word i = 0; i < rnames_size; i++){
		aux_length = LoadValue<cds_word>(fp);
		aux_text = LoadValue<char>(fp, aux_length);
		ref.assign(aux_text, aux_length);
		crps->RNames.push_back(ref);
		delete [] aux_text;
		crps->RNamesCluster.push_back(LoadValue<cds_word>(fp));
	}	
	crps->huffA = Huffman::Load(fp); 
	crps->huffC = Huffman::Load(fp);
	crps->huffG = Huffman::Load(fp);  
	crps->huffT = Huffman::Load(fp);
	crps->init_preseq = LoadValue<size_t>(fp);
	crps->init_seq = LoadValue<size_t>(fp);
	crps->end_seq = LoadValue<size_t>(fp);
	delete [] Variables;	
	fp.seekg(crps->end_seq);
	return crps;
}

