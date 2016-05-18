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

#include "./../src/crps/CRePoSe.h"

int main (int argc, char *argv[]){
	CRePoSe *cseq;
	cds_word sample = 0;
	
	if(argc ==1){
		cout << "Use: ./CompressSeq <arch> <opt>" << endl;
		cout << "<arch> must be a .sam or .rps (rname pos seq) file" << endl;  
		cout << "opt: " << endl;
		cout << "-s sample:  size of the sample rate that will be use. Default: no sample" << endl;
		return 0;
	}

	string filename = argv[1];
	string format = filename.substr(filename.size() - 3);
	ofstream fileCSeq;
	string file = ".cseq";
	file = filename + file;

	int c;
	while((c = getopt (argc, argv, "s:")) != -1){
		switch (c){
			case 's':	sample = atoi(optarg);	break;
			case '?': if(optopt == 's') fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else fprintf(stderr,"Unknown option character `\\x%x'.\n",	optopt);
								return 1;
			default:	abort ();
		}
	}
	
	cseq = new CRPSPreMF(sample); 

	/*get the SEQ fields*/
	string line;
	vector<string> tokens; 
	ifstream SamFile(filename.c_str());
	if (!SamFile.is_open()){
		cout << "Unable to open file" << endl;
		return 1;
	}
	if(format == "sam"){
		cout << "format: sam" << endl;
		getline(SamFile, line);
		while (SamFile.good() ){
			if(line[0] != '@'){		
				Tokenize(line, tokens, "\t");	
				cseq->ProcessLine(tokens[2], atoi(tokens[3].c_str()), tokens[9]);
			}
			tokens.clear();
			getline(SamFile,line);
		}
		fileCSeq.open(file.c_str());

		SamFile.clear();
		SamFile.seekg(0, ios_base::beg);

		cseq->CreateFile(SamFile, fileCSeq, 0);
		fileCSeq.close();
	}
	else if(format == "rps"){
		 cout << "format: rps" << endl;
		getline(SamFile, line);
		while (SamFile.good() ){
			if(line[0] != '@'){		
				Tokenize(line, tokens, "\t");	
				cseq->ProcessLine(tokens[0], atoi(tokens[1].c_str()), tokens[2]);
			}
			tokens.clear();
			getline(SamFile,line);
		}
		fileCSeq.open(file.c_str());

		SamFile.clear();
		SamFile.seekg(0, ios_base::beg);

		cseq->CreateFile(SamFile, fileCSeq, 1);
		fileCSeq.close();
	}
	else{
		fprintf(stderr,"Unknown input format file.\n");
		return 1;
	}

	delete cseq;
	return 0;
}

