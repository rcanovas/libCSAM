#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "./../src/csam/CSAM.h"

int main (int argc, char *argv[]){
	CSAM *csam;
	string filename, file = "_inter.sam";
	string rname = "";
	cds_word pos_x = 0, pos_y = 0;
	string auxSeq = "";
	if(argc!=5){
		cout << "Use: ./GetIntervalSAM <arch> ref_name pos_x pos_y" << endl;
		return 0;
	}
	/*load info*/
	filename = argv[1];
	rname = argv[2];
	pos_x = atoi(argv[3]);
	pos_y = atoi(argv[4]);
	ofstream fileDecomSAM;
	file = filename + file;
	fileDecomSAM.open(file.c_str());
	csam = CSAM::Load(argv[1]);
	cout << "Get Interval: " << rname << " [" << pos_x << ", " << pos_y << "]" << endl;
	csam->GetInterval(fileDecomSAM, rname, pos_x, pos_y);
	fileDecomSAM.close();
	delete csam;
	return 0;
}

