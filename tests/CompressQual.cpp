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

#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "./../src/cqual/CQualBlock.h"

int main (int argc, char *argv[]){
	CQualBlock *cqual;
	int lossyparameter = 0;
	cds_word sample = 0;
	cds_word qualmode = 0;
	cds_word mode = 1;
	
	if(argc ==1){
		cout << "Use: ./CompressQual <arch> <opt>" << endl;
		cout << "opt: " << endl;
		cout << "-q qm: How the Quality values are stored. qm=0 gzip, qm=1 pblock, qm=2 rblock. Default: qm 0" << endl;
		cout << "-m mode: Mode to store the Representative Array. mode = 0 ASCII, mode = 1 Binary Global, mode = 2 Binary Local. Default mode = 1" << endl;
		cout << "-s sample:  size of the sample rate that will be use. Default: no sample" << endl;
		cout << "-l lossy: lossy parameter use to compress the quality score depending on the mode use. Default: 0" << endl; 
		return 0;
	}

	string filename = argv[1];
	ofstream fileCQual;
	string file = ".cqual";
	file = filename + file;

	int c;
	while((c = getopt (argc, argv, "q:l:s:m:")) != -1){
		switch (c){
			case 'q':	qualmode = atoi(optarg); 	break;
			case 's':	sample = atoi(optarg);	break;
			case 'm': mode = atoi(optarg); break;
			case 'l': lossyparameter = atoi(optarg); break;
			case '?': if(optopt == 's') fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else if(optopt == 'q') fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else if(optopt == 'm') fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else if(optopt == 'l') fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else if(isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
								else fprintf(stderr,"Unknown option character `\\x%x'.\n",	optopt);
								return 1;
			default:	abort ();
		}
	}

	switch(qualmode){
		case 0: cqual = new CQualLL(sample); break;
		case 1: cqual = new CQualPBlock(lossyparameter, mode, sample); break;
		case 2: cqual = new CQualRBlock(lossyparameter, mode, sample); break;
		default: cout << "Error: Qual mode selected do not exist" << endl; exit(0);
	}

	/*get the Quality fields*/
	string line;
	vector<string> tokens; 
	ifstream SamFile(filename.c_str());
	if (!SamFile.is_open()){
		cout << "Unable to open file" << endl;
		return 1;
	}
	getline(SamFile, line);
	while (SamFile.good() ){
		if(line[0] != '@'){
			Tokenize(line, tokens, "\t");	
			cqual->ProcessLine(tokens[10]);
		}
		tokens.clear();
		getline(SamFile,line);
	}

	fileCQual.open(file.c_str());
	cqual->CreateFile(fileCQual);
	fileCQual.close();

	delete cqual;
	return 0;
}

