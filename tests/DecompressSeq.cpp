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
	CRePoSe *crps;
	
	string filename, file = ".seq";
	string auxSeq = "";
	if(argc!=2){
		cout << "Use: ./DecompressSeq <arch>" << endl;
		return 0;
	}
	/*load info*/
	ifstream fileSeq;
	filename = argv[1];
	fileSeq.open(argv[1]);
	ofstream fileDecomSeq;
	file = filename + file;
	fileDecomSeq.open(file.c_str());
	
	crps = CRePoSe::Load(fileSeq);

  crps->DecompressSequence(fileSeq, fileDecomSeq);
	
	fileSeq.close();
	fileDecomSeq.close();
	delete crps;
	return 0;
}

