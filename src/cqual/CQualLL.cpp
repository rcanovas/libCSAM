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
#include "./CQualLL.h"


CQualLL::CQualLL(cds_word sample){
	IndexRate = sample;
	init_qual = end_qual = 0;
	CQual = "";
	tmp_file_string = "";
	buffer = NULL;
	buffer_use = 0;
	numberOfLines = 0;
	IndexLine = NULL;
}

CQualLL::~CQualLL(){
	if(IndexLine != NULL)
		delete [] IndexLine;
}

void CQualLL::ProcessLine(string qual){
	if(numberOfLines == 0){
		//tmp_name = tmpnam (NULL);
		size_t pi = getpid();  //id of the process
		tmp_name = "tmp_llqual_" + to_string(pi); //should work
		fileQual.open(tmp_name.c_str()); //create temporal file
	}
	CQual += qual + "\n";
	if(((numberOfLines + 1) % IndexRate) == 0){
		CQual = GZcompressString(CQual);
		SaveValue(fileQual, CQual.length());
		SaveValue(fileQual, (char *)CQual.c_str(), CQual.length());
		CQual = "";
	}
	numberOfLines ++;
}

void CQualLL::CreateFile(ofstream &fp){
	char *block;
  ifstream fileTemp;
	size_t IndexPos = 0, pos_input = 0, end_input = 0, aux = 0, inpos = 1;
	if(CQual.length() != 0){
		CQual = GZcompressString(CQual);
		SaveValue(fileQual, CQual.length());
		SaveValue(fileQual, (char *)CQual.c_str(), CQual.length());
		CQual = "";
	}
	fileQual.close();
	SaveValue(fp, LOSSLESS);
	cds_word Variables[2];
	Variables[0] = IndexRate;
	Variables[1] = numberOfLines;
	SaveValue(fp, Variables, 2);
	buffer = init_buffer();
	buffer_use = 0;
	if(IndexRate > 0){//this will be overwrite at the end
		IndexPos = fp.tellp();
		IndexLine = new cds_word[1 + (numberOfLines + IndexRate - 1)/IndexRate];
		for(cds_word w = 0; w < (1 + (numberOfLines + IndexRate - 1)/IndexRate); w++)
			IndexLine[w] = 0;
		SaveValue(fp, IndexLine, 1 + ((numberOfLines + IndexRate - 1)/IndexRate));
	}
	init_qual = fp.tellp();
	init_qual += 2 * sizeof(size_t);
	SaveValue(fp, init_qual);
	SaveValue(fp, end_qual);

	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	if(IndexRate > 0)	
		IndexLine[0] = fp.tellp();
	while(pos_input < end_input){
		aux = LoadValue<size_t>(fileTemp);
		block = LoadValue<char>(fileTemp, aux);
		pos_input += aux + sizeof(size_t);
		SaveValue(fp, block, aux);
		if(IndexRate > 0){
			IndexLine[inpos] = fp.tellp();
			inpos ++;
		}
		delete [] block;
	}
	fileTemp.close();
	end_qual = fp.tellp();
	fp.seekp(IndexPos);
	if(IndexRate > 0)
		SaveValue(fp, IndexLine, (1 + (numberOfLines + IndexRate - 1)/IndexRate));
	delete [] buffer;
	fp.seekp(init_qual - sizeof(size_t));
	SaveValue(fp, end_qual);
	fp.seekp(end_qual);
	remove(tmp_name.c_str());
}

CQualLL * CQualLL::Load(ifstream &fp){
	cds_word sizeAux = 0;
	cds_word r = LoadValue<cds_word>(fp);
	if (r != LOSSLESS) {
		assert(false);
		return NULL;
	}
	CQualLL *cqual = new CQualLL(0);
	cds_word *Variables = LoadValue<cds_word>(fp, 2);
	cqual->IndexRate = Variables[0];
	cqual->numberOfLines = Variables[1];
	if(cqual->IndexRate > 0){
		sizeAux = 1 + ((cqual->numberOfLines + cqual->IndexRate - 1) / cqual->IndexRate);
		cqual->IndexLine = LoadValue<cds_word>(fp, sizeAux);
	}
	cqual->init_qual = LoadValue<size_t>(fp);
	cqual->end_qual = LoadValue<size_t>(fp);
	delete [] Variables;
	fp.seekg(cqual->end_qual);
	return cqual;
}


void CQualLL::DecompressQual(ifstream &input, ofstream &output){
	string qual_block = "";
	char *buff = NULL;
	size_t posBlock = IndexLine[0], posBuff = 0;
	cds_word blockNum = 1, sizeIndex = 1 + ((numberOfLines + IndexRate - 1) / IndexRate);
	cds_word sizeAux = 0, init_buff = 0;  //bloack size;
	while(blockNum < sizeIndex){
		if((buff == NULL)  || ((init_buff + MAXRAMBYTES) <= IndexLine[blockNum])){
			init_buff = posBlock; 
			posBuff = 0;
			input.seekg(posBlock);
			sizeAux = (end_qual  - posBlock); //get block of data
			if(sizeAux > MAXRAMBYTES)                                   
				sizeAux = MAXRAMBYTES;                                
			if(buff != NULL)   
				delete [] buff;                                    
			buff = LoadValue<char>(input, sizeAux);                                     
		}
		sizeAux = IndexLine[blockNum] - posBlock;
		qual_block.assign(&buff[posBuff], sizeAux);
		qual_block = GZdecompressString(qual_block);
		output << qual_block;
		qual_block = "";
		posBuff += sizeAux;
		posBlock = IndexLine[blockNum];
		blockNum ++;
	}
}
 

vector<string> CQualLL::GetInterval(ifstream &input, cds_word x, cds_word y, cds_word **buff, cds_word *init_buff, size_t MaxBytes){
	vector<string> block, auxBlock;
	string qual = "";
	cds_word lastLine = y, lastRead = x;
	char *tmp_buff;
	size_t posBlock = 0, numBlock = 0, diff = 0, posBuff = 0, sizeAux = 0;
	if(y >= numberOfLines)
		lastLine = numberOfLines - 1;
	//go to the closer sample point to x
	numBlock = x / IndexRate;
	diff = x - numBlock * IndexRate;
	while(lastRead <= lastLine){
		posBlock = IndexLine[numBlock];
		if((*buff == NULL) || (posBlock < *init_buff)  || ((*init_buff + MaxBytes) <= IndexLine[numBlock + 1])){ 
			*init_buff = posBlock;
			input.seekg(posBlock);
			sizeAux = (end_qual - posBlock); //get block of data
			if(sizeAux > MaxBytes)
				sizeAux = MaxBytes;
			if(*buff != NULL)
				delete [] (*buff);
			*buff = (cds_word *)(LoadValue<char>(input, sizeAux));
		}
		tmp_buff = (char *)(*buff);
		posBuff = posBlock - *init_buff;
		sizeAux = IndexLine[numBlock + 1] - posBlock;
		qual.assign(&tmp_buff[posBuff], sizeAux);
		posBuff += sizeAux;
		Tokenize(GZdecompressString(qual), auxBlock, "\n");
		for(cds_word j = diff; j < IndexRate && (lastRead + j - diff) <= lastLine; j++)
			block.push_back (auxBlock[j]);
		lastRead += IndexRate - diff;
		auxBlock.clear();
		numBlock ++;
		diff = 0;
	}
	return block;
}


