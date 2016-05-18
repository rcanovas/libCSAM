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
	if(argc ==1){
		cout << "Use: ./CompressSAM <arch> <opt>" << endl;
		cout << "opt: " << endl;
		 cout << "-q qm: How the Quality values are stored. q = 0 lossles q = 1 pblock, 2 rblock. Default: q = 1" << endl;
		 cout << "-m mode: Mode use to store the Representative Array. mode = 0 ASCII, 1 Binary Global, 2 Binary Local. Default mode = 1" << endl;
		 cout << "-l lossy: lossy parameter use to compress the quality score depending on the mode use. Default: 0" << endl;
		cout << "-s sample: Sample rate used for Fields and Quality structure. Default: s = 1000" << endl;
		cout << "-p position: Sample position rate used for Seq, Rname, and Pos. Default: p = 1000" << endl;
		return 0;
	}

	string fileName = argv[1];
	cds_word sampleQF = 1000, sampleRPS = 1000;
	int lossyparameter = 0;
	cds_word qualmode = 1, encodemode = 1;

	CSAM *csam;

	int c;
	while((c = getopt (argc, argv, "q:m:l:s:p:")) != -1){
		switch (c){
			case 'q': qualmode = atoi(optarg);  break;
			case 'm': encodemode = atoi(optarg); break;
			case 'l': lossyparameter = atoi(optarg); break;
			case 's':	sampleQF = atoi(optarg);	break;
			case 'p': sampleRPS = atoi(optarg);  break;
			case '?': if(optopt == 's' || optopt == 'p' || optopt == 'l' || optopt == 'm' || optopt == 'q') 
									fprintf (stderr, "Option -%c requires an argument.\n", optopt);
								else 
									fprintf(stderr,"Unknown option character `\\x%x'.\n",	optopt);
								return 1;
			default:	abort ();
		}
	}

	csam = new CSAM(0, qualmode, encodemode, lossyparameter, sampleRPS, sampleQF);
	csam->CompressSAM(fileName);

//	delete csam;
	return 0;
}

