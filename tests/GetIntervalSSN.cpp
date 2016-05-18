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

struct Query{
	string rname;
	int pos_x;
	int pos_y;
};

bool compareSSN(Query first, Query second){
	cds_word i=0;
	string f = first.rname;
	string s = second.rname;
	while( (i < f.length()) && (i < s.length()) ){
		if (f[i] < s[i])
			return true;
		else if(f[i] > s[i])
			return false;
		++i;
	}
	if(f.length() < s.length())
		return true;
	else if(f.length() > s.length())
		return false;
	else{
		if(first.pos_x < second.pos_x)
			return true;
		else
			return false;
	}
}



int main (int argc, char *argv[]){
	Timer *t;
	CSAM *csam;
	vector<Query> queries;
	Query aux_query;
	string filename, file = "_inte.sam";
	string prev_rname = "";
	int prev_x = 0, prev_y = 0;
	string line="";
	vector<string> tokens;
	bool header = true;
	cds_word cont = 0;
	if(argc!=3){
		cout << "Use: ./GetIntervalSSN <arch>.csam sample_interval_file" << endl;
		return 0;
	}

	t = new Timer();
	/*load info*/
	ifstream fileSample;
	filename = argv[1];
	fileSample.open(argv[2]);
	ofstream fileDecom;
	file = filename + file;
	fileDecom.open(file.c_str());
	csam = CSAM::Load(argv[1]);
	
	getline(fileSample, line);
	while (fileSample.good()){
		Tokenize(line, tokens, "\t");
		aux_query.rname = tokens[0];
		aux_query.pos_x = atoi(tokens[1].c_str());
		if((cds_word)aux_query.pos_x < csam->max_DN)
			aux_query.pos_x = 0;
		else
			aux_query.pos_x -= csam->max_DN;
		aux_query.pos_y = atoi(tokens[2].c_str());
		queries.push_back(aux_query);
		getline(fileSample, line);
		tokens.clear();
		cont ++;
	}
	cout << "Cont queries: " << cont << endl;
	sort(queries.begin(), queries.end(), compareSSN);
	prev_rname = queries[0].rname;
	prev_x = queries[0].pos_x;
	prev_y = queries[0].pos_y;
	cont = 0;
	for(cds_word i = 1; i < queries.size(); i++){
		if(queries[i].rname != prev_rname){
			//4095 == all
			//1023 all less other and qual
			//575 minimal (QN FL RN PO MQ SEQ)
			csam->GetInterval(fileDecom, prev_rname, prev_x, prev_y, header, 575);
			cont ++;
			prev_rname = queries[i].rname;
			prev_x = queries[i].pos_x;
			prev_y = queries[i].pos_y;
		}
		else{
			if(queries[i].pos_x > prev_y){
				csam->GetInterval(fileDecom, prev_rname, prev_x, prev_y, header, 575);
				cont ++;
				prev_x = queries[i].pos_x;
				prev_y = queries[i].pos_y;
			}
			else
				prev_y = max(queries[i].pos_y, prev_y);
		}
		header = false;
	}
	csam->GetInterval(fileDecom, prev_rname, prev_x, prev_y, header, 575);
	cout << "Cont after queries: " << cont + 1 << endl;

	fileDecom.close();
	delete csam;
	t->Stop();
	delete t;
	return 0;
}

