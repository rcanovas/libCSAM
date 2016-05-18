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
#include "./CSAM.h"

CSAM::CSAM(int empty, cds_word qualmode, cds_word encodemode, int lossyp, cds_word sampleRPS, cds_word sampleQF){
	if(!empty){
		sampleRate = sampleQF;
		posSample = sampleRPS;
		crps = new CRPSPreMF(posSample);
		switch(qualmode){  
			case 0: cqual = new CQualLL(sampleRate); break;
			case 1: cqual = new CQualPBlock(lossyp, encodemode, sampleRate); break;
			case 2: cqual = new CQualRBlock(lossyp, encodemode, sampleRate); break;				
			default: cout << "Error: Qual mode selected do not exist" << endl; exit(0);
		}
		cfields = new CFields(sampleRate);
	}
	max_DN = 0;
	fileCSAM = NULL;
}
		
CSAM::~CSAM(){
	delete crps;
	delete cqual;
	delete cfields;
	if(fileCSAM != NULL){
		for(int i = 0; i < 3; i++)
			fileCSAM[i].close();
		delete [] fileCSAM;
	}
}


void CSAM::CompressSAM(string filename){
	cds_word lineNumber = 0;
	cds_word *index;
	string line, other;
	string header = "";
	size_t bef = 0, aft = 0;
	size_t aux = 0, initIndex = 0, aux_end = 0;
	vector<string> tokens;
	ifstream SamFile(filename.c_str());
	ofstream CSAMFile;
	string compressFileName = filename + ".csam";	
	if (!SamFile.is_open()){
		cout << "Unable to open " << filename << " file" << endl;
		return;
	}
	CSAMFile.open(compressFileName.c_str());
	SaveValue(CSAMFile, sampleRate);
	SaveValue(CSAMFile, posSample);
	aux = CSAMFile.tellp();
	SaveValue(CSAMFile, NumberofLines);
	SaveValue(CSAMFile, max_DN);
	SaveValue(CSAMFile, aux); //pointer to start of IndexClusterToLines
	getline(SamFile, line);
	while (SamFile.good() ){
		if(line[0] == '@'){ //part of the HEADER
			header += line + "\n";
		}
		else{
			if(header.compare("") != 0){
				cheader = GZcompressString(header);
				SaveValue(CSAMFile, cheader.length());
				SaveValue(CSAMFile, (char *)cheader.c_str(), cheader.length());
				cheader = header = "";
			}
			Tokenize(line, tokens, "\t");
		  crps->ProcessLine(tokens[2], atoi(tokens[3].c_str()), tokens[9], &IndexClusterToLines);
			cqual->ProcessLine(tokens[10]);
			if(tokens.size() > 11){
				other = tokens[11];
				for(cds_word i = 12; i < tokens.size(); i++)
					other += "\t" + tokens[i];
			}
			else
				other = "";
			cfields->ProcessLine(tokens[0], tokens[1], tokens[4], tokens[5], tokens[6], 
													 tokens[7], tokens[8], other, CSAMFile, lineNumber);

			CheckCIGAR(tokens[5], tokens[9].length());
			tokens.clear();
			lineNumber ++;
		}
		getline(SamFile,line);
	}
	NumberofLines = lineNumber;
	cfields->CreateFile(CSAMFile);
	bef = CSAMFile.tellp();
	cqual->CreateFile(CSAMFile);
	aft =  CSAMFile.tellp();
	cout << "Max_DN: " << max_DN << endl;
	cout << "Qual: " << (aft - bef) * 1.0 / 1048576 << "  MB" << endl;
	
	SamFile.clear();
	SamFile.seekg(0, ios_base::beg);

	bef = CSAMFile.tellp();
	crps->CreateFile(SamFile, CSAMFile, 0);
	aft =  CSAMFile.tellp();
	cout << "RPS: " << (aft - bef) * 1.0 / 1048576 << "  MB" << endl;
	initIndex = CSAMFile.tellp();
	//add info of IndexClusterToLines and pointer to recovered when load
	index = &IndexClusterToLines[0]; //Index data
	SaveValue(CSAMFile, (cds_word)IndexClusterToLines.size());
	SaveValue(CSAMFile, index, IndexClusterToLines.size());
	aux_end = CSAMFile.tellp();
	CSAMFile.seekp(aux);
	SaveValue(CSAMFile, NumberofLines);
	SaveValue(CSAMFile, max_DN);
	SaveValue(CSAMFile, initIndex);
	CSAMFile.seekp(aux_end);
	SamFile.close();
	CSAMFile.close();
}

void CSAM::CheckCIGAR(string cigar, cds_word len){
	vector<CigOp> info;
	cds_word l = len;
	info = cfields->CigarToArray(cigar);
	for(cds_word i = 0; i < info.size(); i ++){
		if(info[i].op == 'D' || info[i].op == 'N')
			l += info[i].value;
	}
	if(l > max_DN)
		max_DN = l;
}


void CSAM::DecompressSAM(ofstream &output){
	output << GZdecompressString(cheader);
	crps->DecompressAll(&fileCSAM, output, cqual, cfields);
}


cds_word CSAM::GetInterval(ofstream &output, string rname, cds_word pos_x, cds_word pos_y, bool header, int selection, size_t MaxBytes){
	if(header)
		output << GZdecompressString(cheader);
	return crps->GetIntervalAll(&fileCSAM, output, rname, pos_x, pos_y, cqual, cfields, &IndexClusterToLines, selection, MaxBytes);
}

cds_word CSAM::CountReads(string rname, cds_word pos_x, cds_word pos_y, int strand){
	return crps->CountReads(&fileCSAM, rname, pos_x, pos_y, strand, cfields, &IndexClusterToLines);
}

CSAM * CSAM::Load(string filename){
	CSAM *csam = new CSAM(1);
	csam->fileCSAM = new ifstream[3];
	for(int i = 0; i < 3; i++)
		 (csam->fileCSAM)[i].open(filename.c_str());
	cds_word sizeIndex;
	cds_word *index;
	size_t aux, init_index;
	char *tmpString;
	csam->sampleRate = LoadValue<cds_word>((csam->fileCSAM)[0]);
	csam->posSample = LoadValue<cds_word>((csam->fileCSAM)[0]);
	csam->NumberofLines = LoadValue<cds_word>((csam->fileCSAM)[0]);
	csam->max_DN = LoadValue<cds_word>((csam->fileCSAM)[0]);
	init_index = LoadValue<size_t>((csam->fileCSAM)[0]);
	aux = LoadValue<size_t>((csam->fileCSAM)[0]); 	//we are assuming that all SAM files have a header
	tmpString = LoadValue<char>((csam->fileCSAM)[0], aux); (csam->cheader).assign(tmpString, aux); delete [] tmpString;
	csam->cfields = CFields::Load((csam->fileCSAM)[0]);
	csam->cqual = CQualBlock::Load((csam->fileCSAM)[0]);
	csam->crps = CRePoSe::Load((csam->fileCSAM)[0]);
	//move to init_index and load IndexClusterToLines
	(csam->fileCSAM)[0].seekg(init_index);
	sizeIndex = LoadValue<cds_word>((csam->fileCSAM)[0]);
	index = LoadValue<cds_word>((csam->fileCSAM)[0], sizeIndex);
	csam->IndexClusterToLines.assign(index, index + sizeIndex);
	return csam;
}
																		
cds_word CSAM::GetNumberOfLines(){
	return cqual->GetNumberOfLines();
}
																						

