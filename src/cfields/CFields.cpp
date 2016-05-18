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

#include "./CFields.h"


CFields::CFields(cds_word sample){
	rate = sample;
	init_index = end_field = 0;
	no_other = 0;
	Cqname = Cflag = Cmapq = Ccigar = Crnext = Cpnext = Ctlen = Cother = "";
	memCigar = memFlag = memMapq = memQname = memRnext = memPnext = memTlen = memOther = 0;
}

CFields::~CFields(){

}

void CFields::ProcessLine(string qname, string flag, string mapq, string cigar, string rnext, string pnext, 
		                      string tlen, string others, ofstream &fp, cds_word lineNumber){
	if(lineNumber == 0){
		SaveValue(fp, rate);
		init_index = fp.tellp();
		SaveValue(fp, init_index);
		SaveValue(fp, end_field);
		IndexBlock.push_back(fp.tellp());
	}
	Cqname += qname + "\n"; Cflag += flag + "\n"; Cmapq += mapq + "\n"; Ccigar += cigar + "\n";
	Crnext += rnext + "\n"; Cpnext += pnext + "\n"; Ctlen += tlen + "\n"; Cother += others + "\n";
	if(((lineNumber + 1) % rate) == 0)
		SaveBlock(&Cqname, &Cflag, &Cmapq, &Ccigar, &Crnext, &Cpnext, &Ctlen, &Cother, fp);
}


void CFields::SaveBlock(string *qname, string *flag, string *mapq, string *cigar, string *rnext, 
											  string *pnext, string *tlen, string *others, ofstream &fp){
	*cigar = GZcompressString(*cigar);
	SaveValue(fp, (char *)((*cigar).c_str()), (*cigar).length());
	memCigar +=  (*cigar).length(); 
	*cigar = "";
	IndexFlag.push_back(fp.tellp());
	*flag = GZcompressString(*flag);
	SaveValue(fp, (char *)((*flag).c_str()), (*flag).length());
	memFlag +=  (*flag).length();  
	*flag = "";
	IndexMapq.push_back(fp.tellp());
	*mapq = GZcompressString(*mapq);
	SaveValue(fp, (char *)((*mapq).c_str()), (*mapq).length()); //we do not need to store its length
	memMapq +=  (*mapq).length();  
	*mapq = "";
	IndexQname.push_back(fp.tellp());
	*qname = GZcompressString(*qname);
	SaveValue(fp, (char *)((*qname).c_str()), (*qname).length());
	memQname +=  (*qname).length();  
	*qname = "";
	IndexRnext.push_back(fp.tellp());
	*rnext = GZcompressString(*rnext);
	SaveValue(fp, (char *)((*rnext).c_str()), (*rnext).length());   
	memRnext +=  (*rnext).length();
	*rnext = "";
	IndexPnext.push_back(fp.tellp());
	*pnext = GZcompressString(*pnext);
	SaveValue(fp, (char *)((*pnext).c_str()), (*pnext).length());   
	memPnext +=  (*pnext).length();
	*pnext = "";
	IndexTlen.push_back(fp.tellp());
	*tlen = GZcompressString(*tlen);                  
	SaveValue(fp, (char *)((*tlen).c_str()), (*tlen).length());  //we do not need to store its length
	memTlen +=  (*tlen).length();     
	*tlen = "";
	if((*others).length() != 0){
		IndexOthers.push_back(fp.tellp());
		*others = GZcompressString(*others);
		SaveValue(fp, (char *)((*others).c_str()), (*others).length()); //we do not need to store its length
		memOther += (*others).length();
		*others = "";
	}
	IndexBlock.push_back(fp.tellp());
}


void CFields::CreateFile(ofstream &fp){
	size_t *index_aux;
	size_t aux = init_index;
	cds_word sizeIndex = 0;
	cds_word *index_field;
	if(Ccigar.length() != 0) 
		SaveBlock(&Cqname, &Cflag, &Cmapq, &Ccigar, &Crnext, &Cpnext, &Ctlen, &Cother, fp);
	if(IndexOthers.size() == 0)
		no_other = 1;
	cout << "Cigar: " << (memCigar * 1.0 / 1048576) << "  MB" << endl;
	cout << "Flag: " << (memFlag * 1.0 / 1048576) << "  MB" << endl;
	cout << "Mapq: " << (memMapq * 1.0 / 1048576) << "  MB" << endl;
	cout << "Tlen: " << (memTlen * 1.0 / 1048576) << "  MB" << endl;
	cout << "Qname: " << (memQname * 1.0 / 1048576) << "  MB" << endl;
	cout << "Rnext: " << (memRnext * 1.0 / 1048576) << "  MB" << endl;
	cout << "Pnext: " << (memPnext * 1.0 / 1048576) << "  MB" << endl;
	cout << "OTHER: " << (memOther * 1.0 / 1048576) << "  MB" << endl;

	cds_word mem_index = 0;
	
	sizeIndex = IndexBlock.size();
	//save indexes and pointer to them
	init_index = fp.tellp();
	index_aux = &IndexBlock[0];
	SaveValue(fp,  sizeIndex);
	SaveValue(fp, index_aux, sizeIndex);
	
	index_field = &IndexFlag[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	index_field = &IndexMapq[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	index_field = &IndexQname[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	index_field = &IndexRnext[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	index_field = &IndexPnext[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	index_field = &IndexTlen[0];
	SaveValue(fp, index_field, sizeIndex - 1);
	
	SaveValue(fp, no_other);
	mem_index += (7 * sizeIndex - 6) * sizeof(cds_word);
	if(!no_other){  // IndexOthers.size() == sizeIndex -1
		index_field = &IndexOthers[0];
		SaveValue(fp, index_field, sizeIndex - 1);
		mem_index += (sizeIndex - 1) * sizeof(cds_word);
	}

	cout << "IndexFields: " << (mem_index * 1.0 / 1048576) << "  MB" << endl;
	
	end_field = fp.tellp();
	fp.seekp(aux);
	SaveValue(fp, init_index);
	SaveValue(fp, end_field);
	fp.seekp(end_field);
}

CFields * CFields::Load(ifstream &fp){
	cds_word sizeIndex = 0; 
	cds_word *index;
	CFields *cf = new CFields(0);
	cf->rate = LoadValue<cds_word>(fp);
	cf->init_index = LoadValue<size_t>(fp);
	cf->end_field = LoadValue<size_t>(fp);
	fp.seekg(cf->init_index);
	sizeIndex = LoadValue<cds_word>(fp);
	index = LoadValue<cds_word>(fp, sizeIndex);
	cf->IndexBlock.assign(index, index + sizeIndex);
	delete [] index;

	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexFlag.assign(index, index + sizeIndex - 1);
	delete [] index;
	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexMapq.assign(index, index + sizeIndex - 1);
	delete [] index;
	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexQname.assign(index, index + sizeIndex - 1);
	delete [] index;
	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexRnext.assign(index, index + sizeIndex - 1);
	delete [] index;
	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexPnext.assign(index, index + sizeIndex - 1);
	delete [] index;
	index = LoadValue<cds_word>(fp, sizeIndex - 1);
	cf->IndexTlen.assign(index, index + sizeIndex - 1);
	delete [] index;	
	
	cf->no_other = LoadValue<int>(fp);
	if(!(cf->no_other)){
		index = LoadValue<cds_word>(fp, sizeIndex - 1);
		cf->IndexOthers.assign(index, index + sizeIndex - 1);
		delete [] index;
	}
	fp.seekg(cf->end_field);
	return cf;
}


void CFields::ExtractBlock(ifstream &fp, cds_word blockNum, FieldBlock *fields, char **buff, size_t *init_buff, size_t MaxBytes){
	string data = "";
	char *tmp_buff;
	size_t posBlock = IndexBlock[blockNum], posBuff = 0;
	cds_word sizeAux = 0;  //block size
	if((*buff == NULL) || (posBlock < *init_buff)  || ((*init_buff + MaxBytes) <= IndexBlock[blockNum + 1])){
		*init_buff = posBlock;
		fp.seekg(posBlock);
		sizeAux = (end_field  - posBlock); //get block of data
		if(sizeAux > MaxBytes)
			sizeAux = MaxBytes;
		if(*buff != NULL)
			delete [] (*buff);
		*buff = LoadValue<char>(fp, sizeAux);
	}
	tmp_buff = *buff;
	posBuff = posBlock - *init_buff;
	
	sizeAux = IndexFlag[blockNum] - posBlock;		
	data.assign(&tmp_buff[posBuff], sizeAux);
	Tokenize(GZdecompressString(data), fields->CIGAR, "\n");
	posBuff += sizeAux;
	sizeAux = IndexMapq[blockNum] - IndexFlag[blockNum];
	data.assign(&tmp_buff[posBuff], sizeAux);
	Tokenize(GZdecompressString(data), fields->FLAG, "\n");
	posBuff += sizeAux;
	sizeAux = IndexQname[blockNum] - IndexMapq[blockNum];  
	data.assign(&tmp_buff[posBuff], sizeAux);     
	Tokenize(GZdecompressString(data), fields->MAPQ, "\n");   
	posBuff += sizeAux;
	sizeAux = IndexRnext[blockNum] - IndexQname[blockNum]; 
	data.assign(&tmp_buff[posBuff], sizeAux);
	Tokenize(GZdecompressString(data), fields->QNAME, "\n"); 
	posBuff += sizeAux;
	sizeAux = IndexPnext[blockNum] - IndexRnext[blockNum];  
	data.assign(&tmp_buff[posBuff], sizeAux); 
	Tokenize(GZdecompressString(data), fields->RNEXT, "\n");    
	posBuff += sizeAux;
	sizeAux = IndexTlen[blockNum] - IndexPnext[blockNum];  
	data.assign(&tmp_buff[posBuff], sizeAux);     
	Tokenize(GZdecompressString(data), fields->PNEXT, "\n");   
	posBuff += sizeAux;

	if(no_other)
		sizeAux = IndexBlock[blockNum + 1] - IndexTlen[blockNum]; 
	else
		sizeAux = IndexOthers[blockNum] - IndexTlen[blockNum];
	data.assign(&tmp_buff[posBuff], sizeAux);	
	Tokenize(GZdecompressString(data), fields->TLEN, "\n");   
	posBuff += sizeAux;
	
	if(!no_other){
		sizeAux = IndexBlock[blockNum + 1] - IndexOthers[blockNum];
		data.assign(&tmp_buff[posBuff], sizeAux); 
		Tokenize(GZdecompressString(data), fields->OTHER, "\n");
	}
}

void CFields::ExtractSelection(ifstream &fp, cds_word blockNum, FieldBlock *fields, char **buff, size_t *init_buff, int selection, size_t MaxBytes){
	string data = "";
	char *tmp_buff;
	size_t posBlock = IndexBlock[blockNum], posBuff = 0;
	cds_word sizeAux = 0;  //block size
	if((*buff == NULL) || (posBlock < *init_buff)  || ((*init_buff + MaxBytes) <= IndexBlock[blockNum + 1])){
		*init_buff = posBlock;
		fp.seekg(posBlock);
		sizeAux = (end_field  - posBlock); //get block of data
		if(sizeAux > MaxBytes)
			sizeAux = MaxBytes;
		if(*buff != NULL)
			delete [] (*buff);
		*buff = LoadValue<char>(fp, sizeAux);
	}
	tmp_buff = *buff;
	posBuff = posBlock - *init_buff;
	
	sizeAux = IndexFlag[blockNum] - posBlock;		
	if(selection & 32){
		data.assign(&tmp_buff[posBuff], sizeAux);
		Tokenize(GZdecompressString(data), fields->CIGAR, "\n");
	}
	else
		(fields->CIGAR).assign (rate, "*");
	posBuff += sizeAux;
	sizeAux = IndexMapq[blockNum] - IndexFlag[blockNum];
	if(selection & 2){
		data.assign(&tmp_buff[posBuff], sizeAux);
		Tokenize(GZdecompressString(data), fields->FLAG, "\n");
	}
	else
		(fields->FLAG).assign (rate, "0");
	posBuff += sizeAux;
	sizeAux = IndexQname[blockNum] - IndexMapq[blockNum];  
	if(selection & 16){
		data.assign(&tmp_buff[posBuff], sizeAux);     
		Tokenize(GZdecompressString(data), fields->MAPQ, "\n");   
	}
	else
		(fields->MAPQ).assign (rate, "0");
	posBuff += sizeAux;
	sizeAux = IndexRnext[blockNum] - IndexQname[blockNum]; 
	if(selection & 1){
		data.assign(&tmp_buff[posBuff], sizeAux);
		Tokenize(GZdecompressString(data), fields->QNAME, "\n"); 
	}
	else
		(fields->QNAME).assign (rate, "*");
	posBuff += sizeAux;
	sizeAux = IndexPnext[blockNum] - IndexRnext[blockNum];  
	if(selection & 64){
		data.assign(&tmp_buff[posBuff], sizeAux); 
		Tokenize(GZdecompressString(data), fields->RNEXT, "\n");    
	}
	else
		(fields->RNEXT).assign (rate, "*");
	posBuff += sizeAux;
	sizeAux = IndexTlen[blockNum] - IndexPnext[blockNum];  
	if(selection & 128){
		data.assign(&tmp_buff[posBuff], sizeAux);     
		Tokenize(GZdecompressString(data), fields->PNEXT, "\n");   
	}
	else
		(fields->PNEXT).assign (rate, "0");
	posBuff += sizeAux;
	if(no_other)
		sizeAux = IndexBlock[blockNum + 1] - IndexTlen[blockNum]; 
	else
		sizeAux = IndexOthers[blockNum] - IndexTlen[blockNum];
	if(selection & 256){
		data.assign(&tmp_buff[posBuff], sizeAux);	
		Tokenize(GZdecompressString(data), fields->TLEN, "\n");   
	}
	else
		(fields->TLEN).assign (rate, "0");
	posBuff += sizeAux;
	if(!no_other){
		sizeAux = IndexBlock[blockNum + 1] - IndexOthers[blockNum];
		if(selection & 2048){
			data.assign(&tmp_buff[posBuff], sizeAux); 
			Tokenize(GZdecompressString(data), fields->OTHER, "\n");
		}
		else
			(fields->OTHER).assign (rate, "*");
	}
}

void CFields::ClearFieldBlock(FieldBlock *fields){
	(fields->QNAME).clear();
	(fields->FLAG).clear();
	(fields->MAPQ).clear();
	(fields->CIGAR).clear();
	(fields->RNEXT).clear();
	(fields->PNEXT).clear();
	(fields->TLEN).clear();
	(fields->OTHER).clear();
}


vector<CigOp> CFields::CigarToArray(string cigar){
	vector<CigOp> operations;
	CigOp aux;
	int val = 0;
	for(cds_word i = 0; i < cigar.length(); i++){
		if(cigar[i] >= '0' && cigar[i] <= '9')
			val = val * 10 +  (cigar[i] - '0');
		else{
			aux.op = cigar[i];
			aux.value = val;
			operations.push_back(aux);
			val = 0;
		}
	}
	return operations;
}
