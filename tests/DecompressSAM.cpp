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
#include "./../src/csam/CSAM.h"

int main (int argc, char *argv[]){
	CSAM *csam;
	string filename, file = ".sam";
	
	if(argc!=2){
		cout << "Use: ./DecompressSAM <arch>" << endl;
		return 0;
	}
	/*load info*/
	filename = argv[1];
	ofstream fileDecomSAM;
	file = filename + file;
	fileDecomSAM.open(file.c_str());	
	csam = CSAM::Load(argv[1]);
	csam->DecompressSAM(fileDecomSAM);
	fileDecomSAM.close();
	//delete csam;
	return 0;
}

