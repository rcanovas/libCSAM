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
#include "CQualBlock.h"

CQualBlock::~CQualBlock(){

}

void CQualBlock::storeOnlyQuality(ofstream &fp){
	ifstream fileTemp;
	char *tmp_string = NULL;
	cds_word lineCount = 0, aux = 0, length_block = 0, IndexPos = 1, location = 0, occ = 0;
	size_t pos_input = 0, end_input = 0, postemp = 0;
	int local_mode = (mode == 2);
	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);	
	while(pos_input < end_input){
		if((length_block - postemp) < (size_t)(4 + 2 * local_mode)){ //value + occ(2 bytes) + endline
			fileTemp.seekg(pos_input);
			length_block = (end_input - pos_input);
			if(length_block > MAXRAMBYTES)
				length_block = MAXRAMBYTES;
			if(tmp_string != NULL) 
				delete [] tmp_string;
			tmp_string = LoadValue<char>(fileTemp, length_block);
			postemp = 0;
		}		
		if(local_mode){
			postemp += 2; pos_input += 2;
		}
		SetVarField(buffer, buffer_use, buffer_use + 7, (cds_word)(tmp_string[postemp]));
		buffer_use += 8; location += 8;
		postemp ++;
		if(variableSize){
			occ = (cds_word)((unsigned char)tmp_string[postemp]);
			if(occ >= 128){
				occ = ((occ - 128 ) | (((unsigned char)tmp_string[postemp + 1]) << 7));
				postemp ++; pos_input ++;
			}
			aux = encodeGamma(occ, buffer, buffer_use);//run
			location += aux - buffer_use;
			buffer_use = aux;
		}
		postemp ++;		
		if(tmp_string[postemp] == '\n'){
			check_buffer(buffer, &buffer_use, fp);   //check the buffer
			lineCount ++;
			if(IndexRate > 0 && variableSize){
				if( (lineCount % IndexRate) == 0){		
					IndexLine[IndexPos] = location;
					IndexPos ++;
				}
			}
			postemp ++; pos_input ++;
		}
		pos_input += 2;
	}
	//write whatever is left into buffer
	if(tmp_string != NULL) 
		delete [] tmp_string;
	if(buffer_use != 0){
		postemp = (buffer_use + kWordSize - 1) / kWordSize;
		SaveValue(fp, buffer, postemp);
	}	
}						


void CQualBlock::storeQualityASCII(ofstream &fp){
	ifstream fileTemp;
	char *tmp_string = NULL;
	cds_word lineCount = 0, aux = 0, length_block = 0, IndexPos = 1, location = 0, occ = 0;
	size_t pos_input = 0, end_input = 0, postemp = 0;
	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	while(pos_input < end_input){
		if((length_block - postemp) < 4){ //value + occ(2 bytes) + endline
			fileTemp.seekg(pos_input);
			length_block = (end_input - pos_input);
			if(length_block > MAXRAMBYTES)
				length_block = MAXRAMBYTES;
			if(tmp_string != NULL) 
				delete [] tmp_string;
			tmp_string = LoadValue<char>(fileTemp, length_block);
			postemp = 0;
		}		
		SetVarField(buffer, buffer_use, buffer_use + 7, (cds_word)(tmp_string[postemp]));
		buffer_use += 8; location += 8;
		mem_rep += 8;
		postemp ++;
		occ = (cds_word)((unsigned char)tmp_string[postemp]);
		if(occ >= 128){
			occ = ((occ - 128 ) | (((unsigned char)tmp_string[postemp + 1]) << 7));
			postemp ++; pos_input ++;
		}
		aux = encodeGamma(occ, buffer, buffer_use);//run
		location += aux - buffer_use;
		mem_run += aux - buffer_use;
		buffer_use = aux;
		postemp ++;		
		if(tmp_string[postemp] == '\n'){
			if(variableSize){
				SetVarField(buffer, buffer_use, buffer_use + 7, (cds_word)('\n'));
				buffer_use += 8; location += 8;
			}
			check_buffer(buffer, &buffer_use, fp);   //check the buffer
			lineCount ++;
			if(IndexRate > 0){
				if( (lineCount % IndexRate) == 0){
					IndexLine[IndexPos] = location;
					IndexPos ++;
				}
			}
			postemp ++; pos_input ++;
		}
		pos_input += 2;
	}
	//write whatever is left into buffer
	if(tmp_string != NULL) 
		delete [] tmp_string;
	if(buffer_use != 0){
		postemp = (buffer_use + kWordSize - 1) / kWordSize;
		SaveValue(fp, buffer, postemp);
	}
}


void CQualBlock::storeQualityGlobalmm(ofstream &fp){
	ifstream fileTemp;
	char *tmp_string = NULL;
	cds_word lineCount = 0, aux = 0, length_block = 0, bitsPerValue = 0, shift = 0;
	cds_word IndexPos = 1, location = 0, occ = 0;
	size_t pos_input = 0, end_input = 0, postemp = 0;
	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	if(variableSize){
		bitsPerValue = (cds_word)ceil(log2(1.0 * (globalMax - globalMin + 2)));
		shift = 1;
	}
	else
		bitsPerValue = (cds_word)ceil(log2(1.0 * (globalMax - globalMin + 1)));
	while(pos_input < end_input){
		if((length_block - postemp) < 4){ //value + occ(2 bytes) + endline
			fileTemp.seekg(pos_input);
			length_block = (end_input - pos_input);
			if(length_block > MAXRAMBYTES)
				length_block = MAXRAMBYTES;
			if(tmp_string != NULL)
				delete [] tmp_string;
			tmp_string = LoadValue<char>(fileTemp, length_block);
			postemp = 0;
		}
		SetVarField(buffer, buffer_use, buffer_use + bitsPerValue - 1, (cds_word)(tmp_string[postemp]) - globalMin + shift); //value
		buffer_use += bitsPerValue; location += bitsPerValue;
		mem_rep += bitsPerValue;
		postemp ++;
		occ = (cds_word)((unsigned char)tmp_string[postemp]);
		if(occ >= 128){
			occ = ((occ - 128 ) | (((unsigned char)tmp_string[postemp + 1]) << 7));
			postemp ++;	pos_input ++;
		}
		aux = encodeGamma(occ, buffer, buffer_use);//run
		location += aux - buffer_use;
		mem_run += aux - buffer_use;
		buffer_use = aux;
		postemp ++;
		if(tmp_string[postemp] == '\n'){
			if(variableSize){
				SetVarField(buffer, buffer_use, buffer_use + bitsPerValue - 1, 0);
				buffer_use += bitsPerValue; location += bitsPerValue;
			}
			check_buffer(buffer, &buffer_use, fp);   //check the buffer
			lineCount ++;
			if(IndexRate > 0){
				if( (lineCount % IndexRate) == 0){
					IndexLine[IndexPos] = location;
					IndexPos ++;
				}
			}
			postemp ++; pos_input ++;
		}
		pos_input += 2;
	}
	//write whatever is left into buffer
	if(tmp_string != NULL)
		delete [] tmp_string;
	if(buffer_use != 0){
		postemp = (buffer_use + kWordSize - 1) / kWordSize;
		SaveValue(fp, buffer, postemp);
	}
}


void CQualBlock::storeQualityLocalmm(ofstream &fp){
	ifstream fileTemp;
	int start_line = 1;
	char *tmp_string = NULL;
	cds_word lineCount = 0, aux = 0, length_block = 0, bitsPerValue = 0, shift = 0;
	cds_word IndexPos = 1, location = 0, occ = 0, min = 0, diff = 0;
	size_t pos_input = 0, end_input = 0, postemp = 0;
	fileTemp.open(tmp_name.c_str());
	pos_input = fileTemp.tellg();
	fileTemp.seekg(0, std::ifstream::end);
	end_input = fileTemp.tellg();
	fileTemp.seekg(pos_input);
	while(pos_input < end_input){
		if((length_block - postemp) < (size_t)(4 + 2 * start_line)){ //value + occ(2 bytes) + endline
			fileTemp.seekg(pos_input);
			length_block = (end_input - pos_input);
			if(length_block > MAXRAMBYTES)
				length_block = MAXRAMBYTES;
			if(tmp_string != NULL)
				delete [] tmp_string;
			tmp_string = LoadValue<char>(fileTemp, length_block);
			postemp = 0;
		}
		if(start_line){
			min = (cds_word)(tmp_string[postemp]);
			SetVarField(buffer, buffer_use, buffer_use + 7, min); //local min
			buffer_use += 8; location += 8;
			mem_rep += 8;
			diff = (cds_word)(tmp_string[postemp + 1]) - min;	
			if(variableSize){
				bitsPerValue = (cds_word)ceil(log2(1.0 * (diff + 2)));
				shift = 1;
			}
			else
				bitsPerValue = (cds_word)ceil(log2(1.0 * (diff + 1)));		
			SetVarField(buffer, buffer_use, buffer_use + 2, bitsPerValue); //number of bits per value
			mem_rep += 3;
			buffer_use += 3; location += 3;
			postemp += 2; pos_input += 2;
			start_line = 0;
		}
		SetVarField(buffer, buffer_use, buffer_use + bitsPerValue - 1, (cds_word)(tmp_string[postemp]) - min + shift); //value
		buffer_use += bitsPerValue; location += bitsPerValue;
		mem_rep += bitsPerValue;
		postemp ++;
		occ = (cds_word)((unsigned char)tmp_string[postemp]);
		if(occ >= 128){
			occ = ((occ - 128 ) | (((unsigned char)tmp_string[postemp + 1]) << 7));
			postemp ++;	pos_input ++;
		}
		aux = encodeGamma(occ, buffer, buffer_use);//run
		location += aux - buffer_use;
		mem_run += aux - buffer_use;
		buffer_use = aux;
		postemp ++;
		if(tmp_string[postemp] == '\n'){
			if(variableSize){
				SetVarField(buffer, buffer_use, buffer_use + bitsPerValue - 1, 0);
				buffer_use += bitsPerValue; location += bitsPerValue;
			}
			check_buffer(buffer, &buffer_use, fp);   //check the buffer
			lineCount ++;
			if(IndexRate > 0){
				if( (lineCount % IndexRate) == 0){
					IndexLine[IndexPos] = location;
					IndexPos ++;
				}
			}
			start_line = 1;
			postemp ++; pos_input ++;
		}
		pos_input += 2;
	}
	//write whatever is left into buffer
	if(tmp_string != NULL)
		delete [] tmp_string;
	if(buffer_use != 0){
		postemp = (buffer_use + kWordSize - 1) / kWordSize;
		SaveValue(fp, buffer, postemp);
	}
}


void CQualBlock::DecompressQual(ifstream &input, ofstream &output){	
	cds_word length_block = 0, last_read = 0, pos = 0;
	cds_word countLines = 0;
	size_t pos_input = init_qual;
	if(buffer != NULL)
		delete [] buffer;
	tmp_file_string = "";
	while(countLines < numberOfLines){
		length_block = (end_qual - pos_input);  //get block of data
		if(length_block > MAXRAMBYTES)
			length_block = MAXRAMBYTES;
		input.seekg(pos_input);
		buffer = (cds_word *)(LoadValue<char>(input, length_block));
		pos = last_read;
		while(pos < (length_block * 8) && countLines < numberOfLines){
			tmp_file_string += DecompressLine(buffer, &pos) + "\n";	//recover line by line
			countLines++;
			check_tmp_string(&tmp_file_string, output);
			if((8 * length_block - pos) < MaxLineSize && (pos_input + length_block) != end_qual){	//check if other block of data is neccesary
				pos_input += (pos / 8);
				last_read = pos - 8 * (pos / 8);
				break;
			}
		}
		delete [] buffer;
		buffer = NULL;
	}
	if(tmp_file_string.length() > 0){
		output << tmp_file_string;            
		tmp_file_string = "";           
	}
}

string CQualBlock::DecompressLine(cds_word *data, cds_word *PosData){
	string qual = "";
	cds_word value = 0, occ = 0, len_line = 0, pos = *PosData;
	cds_word jump_line = 10, bits_per_qual = 8, base = 0, shift = 0, min = 0;
	if(numberOfValues == numberOfLines){
		value = GetVarField(data, pos, pos + bits_per_qual -1);             
		pos += bits_per_qual;                     
		if(variableSize){
			pos = decodeGamma(&occ, data, pos);
			qual.append(occ, (char)(value));
		}
		else
			 qual.append(sizeQualLine, (char)(value));
	}
	else{
		if(mode != 2){
			if(mode == 1){
				jump_line = base = globalMin;
				if(variableSize){
					bits_per_qual = (cds_word)ceil(log2(1.0 * (globalMax - globalMin + 2)));
					shift = 1;
				}
				else
					 bits_per_qual = (cds_word)ceil(log2(1.0 * (globalMax - globalMin + 1)));
			}
			if(variableSize){
				value = base + GetVarField(data, pos, pos + bits_per_qual -1);          
				pos += bits_per_qual;
				while(value != jump_line){ 
					pos = decodeGamma(&occ, data, pos);
					qual.append(occ, (char)(value - shift));
					value = base + GetVarField(data, pos, pos + bits_per_qual - 1);          				
					pos += bits_per_qual;
				}
			}
			else{
				while(len_line < sizeQualLine){
					value = base + GetVarField(data, pos, pos + bits_per_qual - 1);	
					pos += bits_per_qual;
					pos = decodeGamma(&occ, data, pos);
					qual.append(occ, (char)value);
					len_line += occ;
				}
			}
		}
		else{ //mode == 2
			min = GetVarField(data, pos, pos + 7);
			pos += 8;
			bits_per_qual = GetVarField(data, pos, pos + 2);        
			pos += 3;
			jump_line = base = min;
			if(variableSize){
				shift = 1;
				value = base + GetVarField(data, pos, pos + bits_per_qual -1);
				pos += bits_per_qual;
				while(value != jump_line){
					pos = decodeGamma(&occ, data, pos);
					qual.append(occ, (char)(value - shift));
					value = base + GetVarField(data, pos, pos + bits_per_qual - 1);
					pos += bits_per_qual;
				}
			}
			else{
				while(len_line < sizeQualLine){
					value = base + GetVarField(data, pos, pos + bits_per_qual - 1);
					pos += bits_per_qual;
					pos = decodeGamma(&occ, data, pos);
					qual.append(occ, (char)value);
					len_line += occ;
				}
			}
		}
	}
	*PosData = pos;
	return qual;
}


vector<string> CQualBlock::GetInterval(ifstream &input, cds_word x, cds_word y, cds_word **buff, cds_word *init_buff, size_t MaxBytes){
	vector<string> block;
	string qual;
	cds_word *tmp_buff = NULL;
	cds_word length_block = 0, last_read = 0, pos = 0;
	cds_word line = 0;
	cds_word LastLine = y;
	size_t pos_input = init_qual;
	bool needData = false;
	if(y >= numberOfLines)
		LastLine = numberOfLines - 1;
	if(IndexLine != NULL){       //go to the closer sample point to x
		 line = x / IndexRate;
		 pos = IndexLine[line];
		 pos_input += (pos / 8);
		 last_read = pos - 8 * (pos / 8);
		 line = line * IndexRate;
	}
	if(numberOfValues == numberOfLines && !variableSize){
		pos_input += x;
		line = x;
	}

	length_block = *init_buff + MaxBytes;
	if((*buff == NULL) || (pos_input < *init_buff) || (length_block <= (pos_input + MaxLineSizeByte)))
		needData = true;
	else{
		if(length_block > end_qual)
			length_block = (end_qual - *init_buff);
		else
			length_block = MaxBytes;
		pos = 8 * (pos_input - *init_buff)  + last_read;
	}
	
	while(line <= LastLine){
		if(needData){
			length_block = (end_qual - pos_input);  //get block of data
			if(length_block > MaxBytes)
				length_block = MaxBytes;
			input.seekg(pos_input);
			*buff = (cds_word *)(LoadValue<char>(input, length_block));
			*init_buff = pos_input;
			pos = last_read;
			needData = false;
		}
		tmp_buff = *buff;
		while(pos < (length_block * 8) && line < x){  //read lines without writing until line = x
			qual = DecompressLine(tmp_buff, &pos);	//recover line by line
			line++;
			if((8 * length_block - pos) < (kWordSize * 350) && (pos_input + length_block) != end_qual){	//check if other block of data is neccesary
				pos_input = *init_buff + (pos / 8);
				last_read = pos - 8 * (pos / 8);
				needData = true;
				break;
			}
		}
		if(line >= x){
			while(pos < (length_block * 8) && line <= LastLine){ 
				qual = DecompressLine(tmp_buff, &pos);  //recover line by line
				block.push_back(qual);
				line++;       
				if((8 * length_block - pos) < (kWordSize * 350) && (pos_input + length_block) != end_qual){  //check if other block of data is neccesary
					pos_input = *init_buff + (pos / 8);              
					last_read = pos - 8 * (pos / 8);
					needData = true;
					break;        
				}
			}
		}
	}
	return block;
}


cds_word CQualBlock::GetNumberOfLines(){
	return numberOfLines;
}


CQualBlock * CQualBlock::Load(ifstream &fp){
	cds_word r = LoadValue<cds_word>(fp);
	size_t pos = fp.tellg();
	fp.seekg(pos - sizeof(cds_word));
	switch (r) {
		case LOSSLESS:
			return CQualLL::Load(fp);
		case PBLOCK: 
			return CQualPBlock::Load(fp);
		case RBLOCK:
			return CQualRBlock::Load(fp);
		default:
			throw CDSException("Unknown type");
	}
	return NULL;
}
