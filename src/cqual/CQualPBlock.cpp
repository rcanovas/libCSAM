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

#include "CQualPBlock.h"


/*Create a CQualPblock object
 * m: 0- Quality values stored as plain ASCII bytes
 *  	1- Quality values stored as binary value using range between min and 	
 * 		   max values over the whole file.
 * 		2- Quality values stored as binary value using range between local min 
 *			 and max values over each quality line
 */
CQualPBlock::CQualPBlock(int lossyparam, int m, cds_word sample){
	mode = m;
	lossy = lossyparam;
	tmp_file_string = "";
	buffer = NULL;
	buffer_use = 0;
	IndexRate = sample;
	numberOfValues = 0;  
	numberOfLines = 0;
	variableSize = 0;
	sizeQualLine = 0;
	globalMax = 0; 
	globalMin = 256;
	init_qual = end_qual = 0;
	IndexLine = NULL;
}
		
CQualPBlock::~CQualPBlock(){
	if(IndexLine != NULL)
		delete [] IndexLine;
}

void CQualPBlock::ProcessLine(string qual){
	cds_word localmin = 256, localmax = 0, occ = 0;
	cds_word min = 0, max = 0, valueL = 0, storequal = 0;
	cds_word pos_local = tmp_file_string.length();
	if(mode == 2)
		tmp_file_string += "  "; //add two spaces for the local min and max
	if(variableSize == 0){
		if(sizeQualLine == 0){ //first line to be read
			//tmp_name = tmpnam (NULL); //generates warnings 
			size_t pi = getpid();  //id of the process
			tmp_name = "tmp_pqual_" + to_string(pi); //should work
			fileQual.open(tmp_name.c_str()); //create temporal file
			sizeQualLine = qual.length();	
		}
		else
			if(qual.length() != sizeQualLine)
				variableSize = 1;
	}
	for(cds_word w = 0; w < qual.length(); w++){
		valueL = (cds_word)(qual[w]);
		if(min == 0){
			min = max = valueL;
			occ = 1;
		}
		else{
			if(valueL >= min && valueL <= max)
				occ++;
			else{
				if(valueL < min  &&  (max - valueL) <= (cds_word)(2 * lossy)){
					occ++;
					min = valueL;
				}
				else if(valueL > max  &&  (valueL - min) <= (cds_word)(2 * lossy)){
					occ++;
					max = valueL;
				}
				else{
					numberOfValues++;
					storequal = (max + min) / 2;
					tmp_file_string += (char)storequal;
					if(occ >= 128){
						tmp_file_string +=  (unsigned char)(occ | (1 << 7));
						tmp_file_string +=  (unsigned char)(occ >> 7);
					}
					else
						tmp_file_string +=  (unsigned char)(occ);  
					if(mode == 1){
						if(storequal < globalMin)
							globalMin = storequal;
						if(storequal > globalMax)
							globalMax = storequal;
					}
					else{ 
						if(mode == 2){
							if(storequal < localmin)
								localmin = storequal;
							if(storequal > localmax)
								localmax = storequal;
						}
					}
					min = max = valueL;
					occ = 1;
				}
			}

		}
	}
	numberOfValues++;
	storequal = (max + min) / 2;
	tmp_file_string += (char)storequal;
	if(occ >= 128){
		tmp_file_string +=  (unsigned char)(occ | (1 << 7));
		tmp_file_string +=  (unsigned char)(occ >> 7);
	}
	else
		tmp_file_string +=  (unsigned char)(occ);
	tmp_file_string += '\n';  // add an end line
	if(mode == 1){            
		if(storequal < globalMin)                         
			globalMin = storequal;                      
		if(storequal > globalMax)                                   
			globalMax = storequal;
	}
	else{ 
		if(mode == 2){     
			if(storequal < localmin)            
				localmin = storequal;         
			if(storequal > localmax)            
				localmax = storequal;
			tmp_file_string[pos_local] = (char)localmin;
			tmp_file_string[pos_local + 1] = (char)localmax;
		}
	}
	check_tmp_string(&tmp_file_string, fileQual);
	numberOfLines ++;
}

void CQualPBlock::CreateFile(ofstream &fp){
	size_t IndexPos = 0;
	if(tmp_file_string.length() > 0){
		fileQual << tmp_file_string;  //write to fileQual (it must exist)
		tmp_file_string = "";
	}
	fileQual.close(); 
	SaveValue(fp, PBLOCK);
	SaveValue(fp, mode);
	SaveValue(fp, variableSize);
	SaveValue(fp, lossy);
	cds_word Variables[6];
	Variables[0] = sizeQualLine;
	Variables[1] = IndexRate;
	Variables[2] = globalMin;
	Variables[3] = globalMax;
	Variables[4] = numberOfValues;
	Variables[5] = numberOfLines;
	SaveValue(fp, Variables, 6);
	//compress the information into fp and erase temporal files
	//use buffer to write every MAXRAM bits. We alloc more memory to avoid overflow writting in the buffer

	mem_run = mem_rep = 0;

	buffer = init_buffer(); 
	buffer_use = 0;
	if(numberOfValues == numberOfLines && !variableSize){
		init_qual = fp.tellp();
		init_qual += 2 * sizeof(size_t);
		SaveValue(fp, init_qual);
		SaveValue(fp, end_qual);
		storeOnlyQuality(fp);
		end_qual = fp.tellp();
	}
	else{
		if(IndexRate > 0){//this will be overwrite at the end
			IndexPos = fp.tellp();
			IndexLine = new cds_word[1 + (numberOfLines + IndexRate - 1)/IndexRate];
			for(cds_word w = 0; w < (1 + (numberOfLines + IndexRate - 1)/IndexRate); w++)
				IndexLine[w] = 0;
			IndexLine[0] = 0;
			SaveValue(fp, IndexLine, 1 + ((numberOfLines + IndexRate - 1)/IndexRate));
		}
		init_qual = fp.tellp();
		init_qual += 2 * sizeof(size_t);
		SaveValue(fp, init_qual);
		SaveValue(fp, end_qual);
		if(numberOfValues == numberOfLines)
			storeOnlyQuality(fp);
		else{
			switch(mode){
				case 0: storeQualityASCII(fp); break;
				case 1: storeQualityGlobalmm(fp); break;
				case 2: storeQualityLocalmm(fp); break;
				default: throw CDSException("Unknown QualPBLock mode type");
			}
		}
		end_qual = fp.tellp();
		fp.seekp(IndexPos);
		if(IndexRate > 0)
			SaveValue(fp, IndexLine, (1 + (numberOfLines + IndexRate - 1)/IndexRate));
	}
	if(buffer != NULL){
		delete [] buffer;
		buffer = NULL;
	}
	fp.seekp(init_qual - sizeof(size_t));
	SaveValue(fp, end_qual);
	fp.seekp(end_qual);
	remove(tmp_name.c_str());
	cout << "Run: " << mem_run*1.0/(8*1048576) << endl;
	cout << "Rep: " << mem_rep*1.0/(8*1048576) << endl;
}

CQualPBlock * CQualPBlock::Load(ifstream &fp){
	cds_word r = LoadValue<cds_word>(fp);
	cds_word sizeAux = 0;
	if (r != PBLOCK) {
		assert(false);
		return NULL;
	}
	CQualPBlock *cqual = new CQualPBlock(0, 0, 0);
	cqual->mode = LoadValue<int>(fp);
	cqual->variableSize = LoadValue<int>(fp);
	cqual->lossy = LoadValue<int>(fp);
	cds_word *Variables = LoadValue<cds_word>(fp, 6);
	cqual->sizeQualLine = Variables[0];
	cqual->IndexRate = Variables[1];
	cqual->globalMin = Variables[2];
	cqual->globalMax = Variables[3];
	cqual->numberOfValues = Variables[4];
	cqual->numberOfLines = Variables[5];
	if(cqual->numberOfValues != cqual->numberOfLines || cqual->variableSize){
		if(cqual->IndexRate > 0){
			sizeAux = 1 + ((cqual->numberOfLines + cqual->IndexRate - 1)/cqual->IndexRate);
			cqual->IndexLine = LoadValue<cds_word>(fp, sizeAux);
		}
	}
	cqual->init_qual = LoadValue<size_t>(fp);
	cqual->end_qual = LoadValue<size_t>(fp);
	delete [] Variables;
	fp.seekg(cqual->end_qual);
	return cqual;
}


