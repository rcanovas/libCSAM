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

void CRPSPreMF::DecompressAll(ifstream **input, ofstream &output, CQualBlock *cqual, CFields *cfields){
	vector<string> vseq, vref;
	string datablock = "", preseq = "", qnafla = "", mcrpt = "", other= "";
	string ref = RNames[0], seq = "";
	cds_word NumBlock = 0, linesDeco = 0, rname_pos = 0, POS = 0, Ini_Block = 0, rel_pos = 0;
	cds_word length_block = 0, last_read_buffer = 0, pos_buffRPS = 0, actualPos = 0;
	cds_word countBlocks = 0, reads_to_copy = 0, posPre = 0;
	bool first = false;
	cds_word last_length = 0;
	size_t pos_input = init_seq;
	((*input)[0]).seekg(pos_input);
	if(buffRPS != NULL)
		delete [] buffRPS;
	length_block = (end_seq - pos_input);  //get block of data
	if(length_block > MAXRAMBYTES)                
		length_block = MAXRAMBYTES;           
	buffRPS = (cds_word *)(LoadValue<char>((*input)[0], length_block));
	pos_buffRPS = 0;
	while(countBlocks < ClusterIndex.size()){
		if( ((rname_pos + 1) < RNames.size()) && (RNamesCluster[rname_pos + 1] == countBlocks)){
			rname_pos ++;
			ref = RNames[rname_pos];
			Ini_Block = 0;
		}
		POS = Ini_Block;
		if(ClusterIndex[countBlocks] != (cds_word)-1){
			ClusterReads = GetVarField(buffRPS, pos_buffRPS,  pos_buffRPS + kWordSize - 1);
			pos_buffRPS += kWordSize;
			posPre = GetVarField(buffRPS, pos_buffRPS,  pos_buffRPS + kWordSize - 1);
			pos_buffRPS += kWordSize;
			first = true;			
			if(variableSize)
				last_length = ClusterOverlap[countBlocks];
			else
				last_length = 0;			
			preseq = GetPressumeSeq((*input)[0], posPre, posPre + IndexRate + 350); //get pressume sequence of the block, assuming 350 maximum read length
			actualPos = 0;
			reads_to_copy = ClusterReads;
			while(reads_to_copy > 0){
				seq = DecompressLine(buffRPS, &pos_buffRPS, preseq, &actualPos, &first, &last_length, &rel_pos);
				vseq.push_back(seq);
				POS += rel_pos;
				vref.push_back(ref + "\t" + NumtoString(POS));
				if(vseq.size() == cfields->rate){
					vqual = cqual->GetInterval((*input)[2], linesDeco, linesDeco + cfields->rate - 1, &buffQual, &init_buffQual, MAXRAMBYTES);
					linesDeco += cfields->rate;
					cfields->ExtractBlock((*input)[1], NumBlock, &fieldblock, &buffField, &init_buffField, MAXRAMBYTES);				
					for(cds_word i = 0; i < vseq.size(); i++)
						datablock += (fieldblock.QNAME)[i] + "\t" + (fieldblock.FLAG)[i] + "\t" + vref[i] + "\t" + 
												 (fieldblock.MAPQ)[i] + "\t" + (fieldblock.CIGAR)[i] + "\t" + (fieldblock.RNEXT)[i] + "\t" + 
												 (fieldblock.PNEXT)[i] + "\t" + (fieldblock.TLEN)[i] + "\t" + vseq[i] + "\t" + 
												 vqual[i] + "\t" + (fieldblock.OTHER)[i] + "\n";  
					check_tmp_string(&datablock, output);
					vseq.clear(); vqual.clear(); vref.clear();
					cfields->ClearFieldBlock(&fieldblock);
					NumBlock ++;
				}
				reads_to_copy --;
				if((8 * length_block - pos_buffRPS) < (kWordSize * 350) && (pos_input + length_block) != end_seq){	//check if other block of data is neccesary
					pos_input += (pos_buffRPS / 8);
					((*input)[0]).seekg(pos_input);
					last_read_buffer = pos_buffRPS - 8 * (pos_buffRPS / 8);
					length_block = (end_seq - pos_input);  //get block of data
					if(length_block > MAXRAMBYTES)                      
						length_block = MAXRAMBYTES;    
					buffRPS = (cds_word *)(LoadValue<char>((*input)[0], length_block)); 
					pos_buffRPS = last_read_buffer;
				}	
			}
		}
		countBlocks ++;
		Ini_Block += IndexRate;
	}	
	if(vseq.size() > 0){
		vqual = cqual->GetInterval((*input)[2], linesDeco, linesDeco + cfields->rate - 1, &buffQual, &init_buffQual, MAXRAMBYTES);
		cfields->ExtractBlock((*input)[1], NumBlock, &fieldblock, &buffField, &init_buffField, MAXRAMBYTES);                    
		for(cds_word i = 0; i < vseq.size(); i++)
			datablock += (fieldblock.QNAME)[i] + "\t" + (fieldblock.FLAG)[i] + "\t" + vref[i] + "\t" +  
									 (fieldblock.MAPQ)[i] + "\t" + (fieldblock.CIGAR)[i] + "\t" + (fieldblock.RNEXT)[i] + "\t" +   
									 (fieldblock.PNEXT)[i] + "\t" + (fieldblock.TLEN)[i] + "\t" + vseq[i] + "\t" +  
									 vqual[i] + "\t" + (fieldblock.OTHER)[i] + "\n"; 
		vseq.clear(); vqual.clear(); vref.clear();        
		cfields->ClearFieldBlock(&fieldblock);
	}
	if(datablock.length() != 0)
		output << datablock;
}


cds_word CRPSPreMF::GetIntervalAll(ifstream **input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, CQualBlock *cqual, CFields *cfields, vector<cds_word>* ClusToLine, int selection, size_t MaxBytes){
	cds_word total_reads = 0, line_count = 0, linesDeco = 0;
	cds_word cluster_ref = (cds_word)-1, next = 0; 
	cds_word POS = 0,  rel_pos = 0, Field_count = 0, NumBlock = 0; 
	cds_word cluster_x, cluster_y, max_cluster, Ini_cluster;
	cds_word pos_buffRPS = 0, reads_in_cluster = 0, posPre = 0;
	cds_word last_length = 0, actualPos = 0;
	string preseq = "", seq = "", QnaFla = "";
	string datablock = "", mcrpt = "", other= "";
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
	if( (Ini_cluster + ClusterOverlap[cluster_x]) > pos_x){ //check if clusteroverlap affect the count
		cluster_x --;
		Ini_cluster -= IndexRate;
	}
	while(cluster_x  <= cluster_y){
		POS = Ini_cluster;
		if(POS > pos_y)
			break;
		if(ClusterIndex[cluster_x] != (cds_word)-1){
			pos_input = init_seq + (ClusterIndex[cluster_x] / 8);
			if((buffRPS == NULL) || (pos_input < init_buffRPS) || ((init_buffRPS + MaxBytes) <= (pos_input + kBytesPerWord * 350))){
				((*input)[0]).seekg(pos_input);
				sizeAux = (end_seq - pos_input); //get block of
				if(sizeAux > MaxBytes)
					sizeAux = MaxBytes; 
				if(buffRPS != NULL)
					delete [] buffRPS;
				buffRPS = (cds_word *)(LoadValue<char>((*input)[0], sizeAux));
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
			preseq = GetPressumeSeq((*input)[0], posPre, posPre + IndexRate + 350); //get pressume sequence of the block, assuming 350 maximum read length
			actualPos = 0;
			reads_in_cluster = ClusterReads;			
			NumBlock = (*ClusToLine)[cluster_x] / cfields->rate;
			Field_count = (*ClusToLine)[cluster_x] - NumBlock * cfields->rate;
			while(reads_in_cluster > 0){				
				seq = DecompressLine(buffRPS, &pos_buffRPS, preseq, &actualPos, &first, &last_length, &rel_pos);
				POS += rel_pos;
				if(Field_count == cfields->rate){ //get FLAG block
					Field_count = 0;
					NumBlock ++;
				}
				if(POS > pos_y)
					break;
				if(POS >= pos_x){
					if(QnaFla == "" || Field_count == 0){						
						if(lastFieldBlock != NumBlock){
							cfields->ClearFieldBlock(&fieldblock);
							cfields->ExtractSelection((*input)[1], NumBlock, &fieldblock, &buffField, &init_buffField, selection);
							vqual.clear();
							if(selection & 1024){
								linesDeco = NumBlock * cfields->rate;
								vqual = cqual->GetInterval((*input)[2], linesDeco, linesDeco + cfields->rate - 1, &buffQual, &init_buffQual);
							}
							else
								vqual.assign (cfields->rate, "*");
						}
						lastFieldBlock = NumBlock;
					}
					total_reads ++;
					datablock += (fieldblock.QNAME)[Field_count] + "\t" + (fieldblock.FLAG)[Field_count] + "\t" + rname + "\t" +
									 		 NumtoString(POS) + "\t" + (fieldblock.MAPQ)[Field_count] + "\t" + (fieldblock.CIGAR)[Field_count] + "\t" + 
											 (fieldblock.RNEXT)[Field_count] + "\t" + (fieldblock.PNEXT)[Field_count] + "\t" + 
											 (fieldblock.TLEN)[Field_count] + "\t" + seq + "\t" + vqual[Field_count] + "\t" + (fieldblock.OTHER)[Field_count] + "\n";
					line_count ++;
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
					((*input)[0]).seekg(pos_input);
					buffRPS = (cds_word *)(LoadValue<char>((*input)[0], sizeAux));
					init_buffRPS = pos_input;
					pos_buffRPS = pos_buffRPS - 8 * (pos_buffRPS / 8);
				}
				Field_count ++;
			}
		}
		cluster_x ++;
		Ini_cluster += IndexRate;
	}
	if(datablock.length() != 0)
		output << datablock;
	return total_reads;

}


cds_word CRPSPreMF::CountReads(ifstream **input, string rname, cds_word pos_x, cds_word pos_y, int reverse, CFields *cfields, vector<cds_word>* ClusToLine, size_t MaxBytes){
	cds_word total_reads = 0;
	vector<string> vqf, small_record;
	cds_word cluster_ref = (cds_word)-1, next = 0; 
	cds_word POS = 0,  rel_pos = 0, Flag_count = 0, NumBlock = 0; 
	cds_word cluster_x, cluster_y, max_cluster, Ini_cluster;
	cds_word pos_buffRPS = 0, reads_in_cluster = 0, posPre = 0;
	cds_word last_length = 0, actualPos = 0;
	string preseq = "", seq = "", QnaFla = "";
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
	if( (Ini_cluster + ClusterOverlap[cluster_x]) > pos_x){ //check if clusteroverlap affect the count
		cluster_x --;
		Ini_cluster -= IndexRate;
	}
	while(cluster_x  <= cluster_y){
		POS = Ini_cluster;
		if(POS > pos_y)
			break;
		if(ClusterIndex[cluster_x] != (cds_word)-1){
			pos_input = init_seq + (ClusterIndex[cluster_x] / 8);
			if((buffRPS == NULL) || (pos_input < init_buffRPS) || ((init_buffRPS + MaxBytes) <= (pos_input + kBytesPerWord * 350))){
				((*input)[0]).seekg(pos_input);
				sizeAux = (end_seq - pos_input); //get block of data
				if(sizeAux > MaxBytes)
					sizeAux = MaxBytes; 
				if(buffRPS != NULL)
					delete [] buffRPS;
				buffRPS = (cds_word *)(LoadValue<char>((*input)[0], sizeAux));
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
			preseq = GetPressumeSeq((*input)[0], posPre, posPre + IndexRate + 350); //get pressume sequence of the block, assuming 350 maximum read length
			actualPos = 0;
			reads_in_cluster = ClusterReads;			
			NumBlock = (*ClusToLine)[cluster_x] / cfields->rate;
			Flag_count = (*ClusToLine)[cluster_x] - NumBlock * cfields->rate;
			while(reads_in_cluster > 0){				
				seq = DecompressLine(buffRPS, &pos_buffRPS, preseq, &actualPos, &first, &last_length, &rel_pos);
				POS += rel_pos;
				if(Flag_count == cfields->rate){ //get FLAG block
					Flag_count = 0;
					NumBlock ++;
				}
				if(POS > pos_y)
					break;
				if((POS  + seq.length()) > pos_x){
					if(QnaFla == "" || Flag_count == 0){
						vqf.clear();
				//		cfields->ExtractQnaFla((*input)[1], NumBlock, &QnaFla, &buffField, &init_buffField);
						Tokenize(QnaFla, vqf, "\n");
					}
					Tokenize(vqf[Flag_count], small_record, "\t");
					if(reverse == 0)
						total_reads ++;
					else{ 
						if( (((atoi((small_record[1]).c_str())) & 0x10) == 16) && (reverse == -1))
							total_reads ++;
						else if( (((atoi((small_record[1]).c_str())) & 0x10) == 0) && (reverse == 1))
							total_reads ++;
					}
					small_record.clear();
				}
				reads_in_cluster --;
				if(((8 * sizeAux - pos_buffRPS) < (kWordSize * 350)) && ((init_buffRPS + sizeAux) != end_seq)){  //check if other block of data is neccesary
					pos_input = init_buffRPS + pos_buffRPS / 8;
					sizeAux = (end_seq - pos_input); //get block of data
					if(sizeAux > MaxBytes)
						sizeAux = MaxBytes;              
					if(buffRPS != NULL)
						delete [] buffRPS;
					((*input)[0]).seekg(pos_input);
					buffRPS = (cds_word *)(LoadValue<char>((*input)[0], sizeAux));
					init_buffRPS = pos_input;
					pos_buffRPS = pos_buffRPS - 8 * (pos_buffRPS / 8);
				}
				Flag_count ++;
			}
		}
		cluster_x ++;
		Ini_cluster += IndexRate;
	}	
	return total_reads;
}

