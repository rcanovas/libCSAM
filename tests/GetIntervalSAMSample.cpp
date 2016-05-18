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

#include "./../src/basics/time.h"
#include "./../src/csam/CSAM.h"

int main (int argc, char *argv[]){
	Timer *t;
	CSAM *csam;
	string filename, file = "_inte.sam";
	string line="", rname="", last_rname = "";
	vector<string> tokens;
	cds_word pos_x = 0, pos_y = 0;
	cds_word last_y = 0;
	string auxSeq = "";
	size_t sizeBuf = MAXRAMBYTESDECO;
	bool header = false, query = true;
	cds_word count = 0, cont = 0;
	if(argc!=4 && argc!=3){
		cout << "Use: ./GetIntervalSAMSample <arch>.csam sample_interval_file" << endl;
		cout << "Use: ./GetIntervalSAMSample <arch>.csam sample_interval_file BuffSizeInBytes" << endl;
		return 0;
	}

	t = new Timer();
	/*load info*/
	ifstream fileSample;
	filename = argv[1];	
	fileSample.open(argv[2]);
	if(argc==4)
		sizeBuf = atoi(argv[3]);
	ofstream fileDecom;
	file = filename + file;
	fileDecom.open(file.c_str());
	csam = CSAM::Load(argv[1]);
	//we will assume that the interval are in order by rname, first position and end postion, 
	//if not just add a sorting function here
	getline(fileSample, line);
	while (fileSample.good()){
		Tokenize(line, tokens, "\t");
		rname = tokens[0];
		pos_x = atoi(tokens[1].c_str());
		pos_y = atoi(tokens[2].c_str());
		if(rname == last_rname){
			if(pos_x <= last_y){
				if(pos_y <= last_y){
					pos_y = last_y;
					query = false;
				}
				else
					pos_x = last_y + 1;
			}
		}
		if(query)
			count += csam->GetInterval(fileDecom, rname, pos_x, pos_y, header, 4095, sizeBuf);
		header = false;
		query = true;
		getline(fileSample,line);
		tokens.clear();
		cont ++;
		last_y = pos_y;
		last_rname = rname;
	}
	fileDecom.close();
	delete csam;
	t->Stop();
//	cout << "Numer of lines: " << count << endl;
	cout << /*"Time per line: " <<*/ t->ElapsedTime() / count << /* " microseconds" <<*/ endl;
	delete t;
	return 0;
}

