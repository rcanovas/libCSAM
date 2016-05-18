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
#include "./CRePoSe.h"

cds_word CRePoSe::GetNumberOfLines(){
	return numberOfLines;
}

CRePoSe * CRePoSe::Load(ifstream &fp){
	cds_word r = LoadValue<cds_word>(fp);
	size_t pos = fp.tellg();
	fp.seekg(pos - sizeof(cds_word));
	switch (r) {
		case ROW:
			return CRPSPreMF::Load(fp);
		default:
			throw CDSException("Unknown type");
	}
	return NULL;
}

void CRePoSe::initData(cds_word sample, bool create){
	tmp_file_string = "";
	buffer = NULL;
	buffer_use = bitsLength = 0; 
	IndexRate = sample;  //must be greater than 0 always.
	variableSize = 0;
	init_seq = init_preseq = end_seq = 0; 
	last_ref = "*";
	ClusterReads = ClusterMaxPos = MaxClusterReads = MaxDepthCluster = 0;
	huffA = huffC = huffG = huffT = NULL;
	last_pos = farther_pos = overlapcluster = preSeqLength = 0;
	posWin = 0;
	if(create){
		WinSeq = new int*[WINSIZE]; //window approach of most frequent base 
		for(cds_word j = 0; j < WINSIZE; j++){
			WinSeq[j] = new int[6];
			for(int u = 0; u < 4; u++)
				WinSeq[j][u] = 0;
			WinSeq[j][4] = -1;
			WinSeq[j][5] = 0; //number of Ns
		}
	}
	//Global variables used to get intervals of the other fields of the SAM file
	buffQual = buffRPS = NULL;
	buffField = NULL;
	init_buffQual = init_buffField = init_buffRPS = 0;
}

void CRePoSe::restartWinSeq(string seq){
	char letter;
	for(cds_word j = 0; j < WINSIZE; j++){
		for(cds_word u = 0; u < 4; u++)
			WinSeq[j][u] = 0;
		WinSeq[j][4] = -1;
		WinSeq[j][5] = 0;
	}
	for(cds_word j = 0; j < seq.length(); j++){
		letter = toupper(seq[j]);
		if(letter == 'A'){      WinSeq[j][0] = 1; WinSeq[j][4] = 0; }
		else if(letter == 'C'){ WinSeq[j][1] = 1; WinSeq[j][4] = 1; }
		else if(letter == 'G'){ WinSeq[j][2] = 1; WinSeq[j][4] = 2; }
		else if(letter == 'T'){ WinSeq[j][3] = 1; WinSeq[j][4] = 3; }
		else{ WinSeq[j][4] = 0; WinSeq[j][5] = 1;}
	}
}

void CRePoSe::getValuesPreSeq(int since, cds_word n){
	int aux = 0, count = 0;
	for(cds_word j = 0; j < n; j++){
		aux = (since + j) % WINSIZE;
		switch(WinSeq[aux][4]){
			case -1: break;
			case 0: tmp_file_string += "A";  preSeqLength ++;
							occ_notA[cds_word('C')] += WinSeq[aux][1];
							occ_notA[cds_word('G')] += WinSeq[aux][2];
							occ_notA[cds_word('T')] += WinSeq[aux][3];
							break;
			case 1: tmp_file_string += "C";  preSeqLength ++;
							occ_notC[cds_word('A')] += WinSeq[aux][0];
							occ_notC[cds_word('G')] += WinSeq[aux][2];
							occ_notC[cds_word('T')] += WinSeq[aux][3];
							break;
			case 2: tmp_file_string += "G";  preSeqLength ++;
							occ_notG[cds_word('A')] += WinSeq[aux][0];
							occ_notG[cds_word('C')] += WinSeq[aux][1];
							occ_notG[cds_word('T')] += WinSeq[aux][3];
							break;
			case 3: tmp_file_string += "T";  preSeqLength ++;
							occ_notT[cds_word('A')] += WinSeq[aux][0];
							occ_notT[cds_word('C')] += WinSeq[aux][1];
							occ_notT[cds_word('G')] += WinSeq[aux][2];
							break;
		}
		count = WinSeq[aux][0] + WinSeq[aux][1] + WinSeq[aux][2] + WinSeq[aux][3] + WinSeq[aux][5];
		if((cds_word)count > MaxDepthCluster)
			MaxDepthCluster = count;
		WinSeq[aux][0] = WinSeq[aux][1] = WinSeq[aux][2] = WinSeq[aux][3] = 0;
		WinSeq[aux][4] = -1;
		WinSeq[aux][5] = 0;
	}
}

void CRePoSe::ProcessLine(string Ref, cds_word Pos, string Seq, vector<cds_word>* IndexClusterToLines){
	cds_word readLength = Seq.length(), l = 0, aux = 0;
	char letter;
	if(sizeLine.size() == 0){  //it means that is the first line to be process
		sizeLine.push_back(readLength);
		//tmp_name = tmpnam (NULL);
		size_t pi = getpid();  //id of the process
		tmp_name = "tmp_crepose_" + to_string(pi); //should work
		fileSeq.open(tmp_name.c_str()); //create temporal file
	}
	else{
		while(l < sizeLine.size()){  //check the line length
			if(readLength == sizeLine[l])
				break;
			l ++;
		}
		if(l == sizeLine.size()){
			sizeLine.push_back(readLength);
			variableSize = 1;
			l = 0;
		}
	}
	if(!Ref.compare("*")){  //we assume that all the lines without ref are at the start
		if(RNames.size() == 0)
			RNames.push_back("*");
		numberOfLinesNR ++;
		cout << "Not referenced reads are ignored so far" << endl;
	}
	else{ //SEQUENCE WITH REFERENCE
		if(!last_ref.compare(Ref)){
			if(farther_pos  > Pos){
				l = Pos - last_pos;
				preSeqLength += l;
				getValuesPreSeq(posWin, l); //copy the part of the Win that do not overlap
				posWin = (posWin +  l) % WINSIZE;
				for(cds_word i = 0; i < readLength; i++){ //copy Seq to the window
					aux = (posWin + i) % WINSIZE;
					letter = toupper(Seq[i]);
					switch(letter){
						case 'A':  WinSeq[aux][0] += 1;
											 if(WinSeq[aux][4] == -1)  WinSeq[aux][4] = 0;
											 else if(WinSeq[aux][0] > WinSeq[aux][WinSeq[aux][4]])  WinSeq[aux][4] = 0; break;
						case 'C':  WinSeq[aux][1] += 1;
											 if(WinSeq[aux][4] == -1)   WinSeq[aux][4] = 1;
											 else if(WinSeq[aux][1] > WinSeq[aux][WinSeq[aux][4]])  WinSeq[aux][4] = 1; break;
						case 'G':  WinSeq[aux][2] += 1;
											 if(WinSeq[aux][4] == -1)  WinSeq[aux][4] = 2;
											 else if(WinSeq[aux][2] > WinSeq[aux][WinSeq[aux][4]])  WinSeq[aux][4] = 2; break;
						case 'T':  WinSeq[aux][3] += 1;
											 if(WinSeq[aux][4] == -1)  WinSeq[aux][4] = 3;
											 else if(WinSeq[aux][3] > WinSeq[aux][WinSeq[aux][4]])  WinSeq[aux][4] = 3; break;
						default :  //N case 
											 if(WinSeq[aux][4] == -1)  WinSeq[aux][4] = 0;  WinSeq[aux][5] += 1; break;
					}
				}
			}
			else{ //NO OVERLAP READ
				getValuesPreSeq(posWin, WINSIZE);
				posWin = 0;
				restartWinSeq(Seq);
			}
			if(Pos >= ClusterMaxPos){  //CHECK Cluster info
				ClusterIndex[ClusterIndex.size() - 1] = ClusterReads;
				while(ClusterMaxPos + IndexRate <= Pos){
					ClusterIndex.push_back((cds_word)-1);
					ClusterOverlap.push_back(overlapcluster);
					overlapcluster = 0;
					ClusterMaxPos += IndexRate;
					if(IndexClusterToLines != NULL)                       
						IndexClusterToLines->push_back(numberOfLines);
				}
				if(IndexClusterToLines != NULL)
					IndexClusterToLines->push_back(numberOfLines);
				ClusterIndex.push_back((cds_word)-1); //change this later for pointer to cluster
				ClusterOverlap.push_back(overlapcluster);
				overlapcluster = 0;
				ClusterMaxPos += IndexRate;
				if(MaxClusterReads < ClusterReads)
					MaxClusterReads = ClusterReads;
				ClusterReads = 1;
			}
			else
				ClusterReads++;
			if(Pos + readLength >= ClusterMaxPos) //check if some of the reads extend to the next
				overlapcluster = max(overlapcluster, Pos + readLength - ClusterMaxPos);
		}
		else{ //NEW REFERENCE
			getValuesPreSeq(posWin, WINSIZE); //get values from WinSeq
			posWin = 0;
			restartWinSeq(Seq);
			RNames.push_back(Ref);
			ClusterMaxPos = 0;
			if(ClusterIndex.size() > 0){
				ClusterIndex[ClusterIndex.size() - 1] = ClusterReads;
				if(overlapcluster){ //add extra block if there is at least one read overlap with the next block
					ClusterIndex.push_back((cds_word)-1);
					ClusterOverlap.push_back(overlapcluster);
					if(IndexClusterToLines != NULL)                       
						IndexClusterToLines->push_back(numberOfLines);
				}
			}
			overlapcluster = 0;
			RNamesCluster.push_back(ClusterIndex.size());
			while(ClusterMaxPos + IndexRate <= Pos){
				ClusterIndex.push_back((cds_word)-1);
				ClusterOverlap.push_back(0);
				ClusterMaxPos += IndexRate;
				if(IndexClusterToLines != NULL)
					IndexClusterToLines->push_back(numberOfLines);
			}
			if(IndexClusterToLines != NULL)
				IndexClusterToLines->push_back(numberOfLines);
			ClusterIndex.push_back((cds_word)-1); //change this later for pointer to block
			ClusterOverlap.push_back(0);
			if(MaxClusterReads < ClusterReads)
				MaxClusterReads = ClusterReads;
			ClusterReads = 1;
			ClusterMaxPos += IndexRate;
			last_ref = Ref;
			farther_pos = 0;
		}
		last_pos = Pos;
		farther_pos = max(farther_pos, Pos + readLength);
	}
	check_tmp_string(&tmp_file_string, fileSeq);
	numberOfLines ++;
}

cds_word CRePoSe::getFirstCluster(string Ref, cds_word *next){
	for(cds_word i = 0; i < RNames.size(); i++){
		if(!Ref.compare(RNames[i])){
			*next = i + 1;
			return RNamesCluster[i];
		}
	}
	return (cds_word)-1;
}


/*we assumed that init_buff is zero*/
void CRePoSe::CompressPressumeSeq(ofstream &fp, cds_word **buff, cds_word *init_buff){
	ifstream fileTemp;
	cds_word letter_value = 0, reading = 0, write_buffer = 0, keep = 0;
	size_t pos_input = 0, end_input = 0;
	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	tmp_file_string = "";
	while(pos_input < end_input){
		reading = (end_input - pos_input);  //get block of data
		if(reading > MAXRAMBYTES)
			reading = MAXRAMBYTES;
		tmp_file_string = LoadValue<char>(fileTemp, reading);
		/*write all letter to the buffer and copy to the file*/
		for(cds_word w = 0; w < tmp_file_string.length(); w++){
			letter_value = (cds_word)((tmp_file_string[w] >> 1) & 0x03);
			SetVarField(*buff, *init_buff,  *init_buff + 1, letter_value);
			*init_buff += 2;
		}
		write_buffer = *init_buff / kWordSize;
		SaveValue(fp, *buff, write_buffer);
		keep = 0;
		while((write_buffer + keep) * kWordSize < buffer_use){
			(*buff)[keep] = (*buff)[write_buffer + keep];
			keep ++;
		}
		*init_buff = *init_buff - write_buffer * kWordSize;
		pos_input += reading;
	}
	if(*init_buff != 0){
		write_buffer = (*init_buff + kWordSize - 1) / kWordSize;
		SaveValue(fp, *buff, write_buffer);
	}
	fileTemp.close();
	fileSeq.close();
}


cds_word CRePoSe::CompressLine(string Ref, cds_word Pos, string Seq, string PSeq, cds_word **buff, cds_word *init_buff){
	cds_word runCo = 0, runRe = 0, pos_enc = 0, Ns = 0;
	cds_word location = 0, aux = 0, relative_pos = 0, readLength = 0;
	vector<cds_word> posN;
	posN.push_back(0);
	char letter;
	bool copy_all = true;
	cds_word len_aux = (Seq.size() * 2 + kWordSize - 1) / kWordSize ;
	cds_word *replace_letters = new cds_word[len_aux];
	cds_word parameter_m = (cds_word)(ceil(0.69 * (IndexRate + ClusterReads) / ClusterReads));
	for(cds_word i = 0; i < len_aux; i ++)
		replace_letters[i] = 0;
	if(!(Ref.compare("*")))
		cout << "Not referenced reads are not consider" << endl;
	else{
		if(!last_ref.compare(Ref)){
			if((cds_word)(Pos / IndexRate) == (cds_word)(last_pos / IndexRate)) //Storing relative positions: if both are in the same Cluster
				relative_pos = Pos - last_pos;
			else
				relative_pos = Pos - ((cds_word)(Pos / IndexRate)) * IndexRate;
		}
		else{
			relative_pos = Pos - ((cds_word)(Pos / IndexRate)) * IndexRate;
			last_ref = Ref;
			farther_pos = 0;
		}
		aux = encodeGolomb(relative_pos, *buff, *init_buff, parameter_m);
		location += aux - *init_buff;
		*init_buff = aux;
		readLength = Seq.length();
		if(variableSize){
			for(cds_word i = 0; i < sizeLine.size(); i++){
				if(readLength == sizeLine[i]){
					SetVarField(*buff, *init_buff, *init_buff + bitsLength - 1, i); //saving read length
					*init_buff += bitsLength; location += bitsLength;
					break;
				}
			}
		}
		if(Seq[0] != PSeq[0] && Seq[0] != 'N'){   //check if we start Coping or Replacing
			BitSet(*buff, *init_buff, 0);
			copy_all = false;
		}
		else
			BitSet(*buff, *init_buff, 1);
		*init_buff = *init_buff + 1; location ++;
		for(cds_word i = 0; i < readLength; i++){
			letter = toupper(Seq[i]);
			if(letter != PSeq[i] && letter != 'N'){
				if(runCo != 0){
					if(copy_all){
						BitSet(*buff, *init_buff, 0);
						copy_all = false;
						*init_buff = *init_buff + 1; location ++;
					}
					aux = encodeGamma(runCo, *buff, *init_buff);
					location += aux - *init_buff;
					*init_buff = aux;
				}
				runCo = 0;
				runRe ++;
				switch(PSeq[i]){
					case 'A':  aux = huffA->encode((int)letter,  replace_letters, pos_enc);
										 break;
					case 'C':  aux = huffC->encode((int)letter,  replace_letters, pos_enc); 
										 break;
					case 'G':  aux = huffG->encode((int)letter,  replace_letters, pos_enc);
										 break;
					default: 	 aux = huffT->encode((int)letter,  replace_letters, pos_enc);
										 break;
				}
				pos_enc = aux;
			}
			else{
				runCo ++;
				if(runRe != 0){
					aux = encodeGamma(runRe, *buff, *init_buff);
					location += aux - *init_buff;
					*init_buff = aux;
					SetVarSeq(*buff, *init_buff,  *init_buff + pos_enc - 1, replace_letters);
					*init_buff += pos_enc;
					location += pos_enc;
					pos_enc = 0;
				}
				runRe = 0;
				if(letter == 'N'){
					Ns ++;
					posN.push_back(i + 1);
				}
			}
		}
		if(runCo != 0){
			if(runCo == readLength){
				BitSet(*buff, *init_buff, 1);
				*init_buff = *init_buff + 1; location ++;
			}
			else{
				aux = encodeGamma(runCo, *buff, *init_buff);
				location += aux - *init_buff;
				*init_buff = aux;
			}
		}
		if(runRe != 0){
			aux = encodeGamma(runRe, *buff, *init_buff);
			location += aux - *init_buff;
			*init_buff = aux;
			SetVarSeq(*buff, *init_buff, *init_buff + pos_enc - 1, replace_letters);
			*init_buff += pos_enc; location += pos_enc;
			pos_enc = 0;
		}
		runCo = runRe = 0;
		last_pos = Pos;
		farther_pos = max(farther_pos, last_pos + readLength);
	}
	if(Ns > 0){
		BitSet(*buff, *init_buff, 1);
		*init_buff = *init_buff + 1;
		aux = encodeGamma(Ns, *buff, *init_buff);
		location += aux -  *init_buff;
		*init_buff = aux;
		if(Ns != Seq.length()){
			for(cds_word i = 1; i <= Ns; i++){
				aux = encodeGamma(posN[i] - posN[i - 1], *buff, *init_buff);
				location += aux - *init_buff;
				*init_buff = aux;
			}
		}
	}
	else{
		BitSet(*buff, *init_buff, 0);
		*init_buff = *init_buff + 1;
	}
	location += 1;
	delete [] replace_letters;
	return location;
}

string CRePoSe::GetPressumeSeq(ifstream &input, cds_word ini, cds_word end){
	string preseq = "";
	cds_word aux = 0, pos = 0, l = 0;
	size_t original_pos = input.tellg();
	size_t move_to =  init_preseq + (cds_word)(ini / 4);
	input.seekg(move_to);
	cds_word *chunck;
	if(end > preSeqLength)
		end = preSeqLength;
	aux = end - 4 * (ini / 4); //2 bits per value
	chunck = LoadValue<cds_word>(input, (aux * 2 + kWordSize - 1) / kWordSize );
	//read block of data
	pos = 2 * (ini - 4 * (ini / 4));
	while(l < (end - ini)){
		aux = GetVarField(chunck, pos,  pos + 1);
		pos += 2;
		if(aux == 0) preseq += 'A';
		else if(aux == 1) preseq += 'C';
		else if(aux == 2) preseq += 'T';
		else preseq += 'G';
		l ++;
	}
	delete [] chunck;
	input.seekg(original_pos);
	return preseq;
}


string CRePoSe::DecompressLine(cds_word *data, cds_word *PosData, string preseq, cds_word *actualPos, bool *first, cds_word *last_length, cds_word *rel_pos){
	string seq = "";
	bool copy = true;
	cds_uint symb;
	cds_word aux = 0, posN = 0;
	cds_word Pos = 0, location = *PosData, sizeRead = sizeLine[0], decoded = 0, cr = 0;
	cds_word parameter_m = (cds_word)(ceil(0.69 * (IndexRate + ClusterReads) / ClusterReads));
	cds_word numN = 0;
	aux = decodeGolomb(&Pos, data, location, parameter_m);
	*rel_pos = Pos;
	location = aux;
	if(variableSize){
		sizeRead = sizeLine[GetVarField(data, location, location + bitsLength - 1)];
		location += bitsLength;
	}
	//get position in the reference sequence
	if(*first){
		if(Pos > *last_length)
			*last_length = sizeRead;
		else
			*last_length = max(sizeRead, *last_length - Pos);
		*first = false;
	}
	else if(Pos > *last_length){
		*actualPos = *actualPos + *last_length;
		*last_length = sizeRead;
	}
	else{
		*actualPos = *actualPos + Pos;
		*last_length = max(sizeRead, *last_length - Pos);
	}
	copy = BitGet(data, location);
	location ++;
	seq = preseq.substr(*actualPos, sizeRead);
	if(copy){
		if(BitGet(data, location))
			decoded = sizeRead;
		location ++;
	}
	while(decoded < sizeRead){
		location = decodeGamma(&cr, data, location);
		if(copy)
			copy = false;
		else{
			copy = true;
			for(cds_word j = 0; j < cr; j ++){
				if(seq[decoded + j] == 'A')
					location = huffA->decode(&symb, data, location);
				else if(seq[decoded + j] == 'C')
					location = huffC->decode(&symb, data, location);
				else if(seq[decoded + j] == 'G')
					location = huffG->decode(&symb, data, location);
				else
					location = huffT->decode(&symb, data, location);
				seq[decoded + j] = (char)symb;
			}
		}
		decoded += cr;
	}
	copy = BitGet(data, location);
	location ++;
	if(copy){ //Ns
		location = decodeGamma(&numN, data, location);
		if(numN == sizeRead){
			seq = "";
			seq.append(sizeRead, 'N');
		}
		else{
			location = decodeGamma(&posN, data, location);
			posN --;
			seq[posN] = 'N';
			for(cds_word i = 1; i < numN; i++){
				location = decodeGamma(&aux, data, location);
				posN += aux;
				seq[posN] = 'N';
			}
		}
	}
	*PosData = location;
	return seq;
}


