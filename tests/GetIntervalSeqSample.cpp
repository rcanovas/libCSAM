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
#include "./../src/crps/CRePoSe.h"

int main (int argc, char *argv[]){
	Timer *t;
	CRePoSe *crps;
	string filename, file = ".rps";
	string line="", rname="";
	vector<string> tokens;
	cds_word pos_x = 0, pos_y = 0;
	cds_word inte_len = 0, max_inte = 0, min_inte = (cds_word)-1, avg_inte = 0;
	string auxSeq = "";
	cds_word count = 0, cont = 0;
	if(argc!=3){
		cout << "Use: ./GetIntervalSeqSample <arch>.cseq sample_interval_file" << endl;
		return 0;
	}
	t = new Timer();
	/*load info*/
	ifstream fileSeq, fileSample;
	filename = argv[1];
	fileSeq.open(argv[1]);
	fileSample.open(argv[2]);
	ofstream fileDecomSeq;
	file = filename + file;
	fileDecomSeq.open(file.c_str());
	crps = CRePoSe::Load(fileSeq);
	
	getline(fileSample, line);
	while (fileSample.good()){
		Tokenize(line, tokens, "\t");
		rname = tokens[0];
		pos_x = atoi(tokens[1].c_str());
		pos_y = atoi(tokens[2].c_str());
		inte_len = pos_y - pos_x;
		avg_inte += inte_len;
		if(inte_len < min_inte)
			min_inte = inte_len;
		if(inte_len > max_inte)
			max_inte = inte_len;
		count += crps->GetInterval(fileSeq, fileDecomSeq, rname, pos_x, pos_y);
		getline(fileSample,line);
		tokens.clear();
		cont ++;
	}
	fileSeq.close();
	fileDecomSeq.close();
	delete crps;
	t->Stop();
	cout << "Numer of lines: " << count << endl;
	cout << "Average Interval: " << (1.0 * avg_inte)/ cont << endl;
	cout << "Min. Interval: " << min_inte << endl;
	cout << "Max. Interval: " << max_inte << endl;
	cout << "Wall clock: " << t->ElapsedTime() << " microseconds" << endl;
	cout << "Time per line: " << t->ElapsedTime() / count <<  " microseconds" << endl;
	cout << "CPU process: " << t->ElapsedTimeCPU() << " microseconds" <<  endl;
	delete t;
	return 0;
}

