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

using namespace std;

int main (int argc, char *argv[]){
	if(argc != 2){
		cout << "Use: ./ComputeEntropyCOL <arch>" << endl;
		return 0;
	}

	string filename = argv[1];

	/*get the Quality fields*/
	string line, qual;
	unsigned long cont = 0, line_count = 0;
	unsigned long letters[256];
	unsigned long posLetter[200][257];

	for(unsigned int i = 0; i < 256; i++)
		letters[i] = 0;

	for(unsigned int i = 0; i < 200; i++)
		for(unsigned int j = 0; j <= 256; j++)
			posLetter[i][j] = 0;


	
	double entropy = 0, avg_en = 0;

	ifstream RefFile(filename.c_str());
	if (!RefFile.is_open()){
		cout << "Unable to open file" << endl;
		return 1;
	}
	
	getline(RefFile, line);
	while (RefFile.good() ){
		for(unsigned int i = 0; i < line.length(); i++){
			letters[(int)(line[i])]++;
			posLetter[i][(int)(line[i])]++;
			cont ++;
			posLetter[i][256]++;
		}
		getline(RefFile, line);
		line_count ++;
	}
	
	RefFile.close();
	
	for(unsigned int i = 0; i < 256; i++)
		if(letters[i] > 0)		
			entropy += ((letters[i]*1.0)/(cont*1.0))*log2((cont*1.0)/(letters[i]*1.0));
	cout << "Number of lines: " << line_count << endl;
	cout << "Total number of quality scores: " << cont << endl;
	cout << "Entropy File: " << entropy << endl;
	cout << "Theoretical Space: " << entropy * cont / 8388608.0 << endl;

	cout << endl;

	cont = 0;
	for(unsigned int i = 0; i < 200; i++){
		entropy = 0;
		for(unsigned int j = 0; j < 256; j++)
			if(posLetter[i][j] > 0)
				entropy += ((posLetter[i][j]*1.0)/(posLetter[i][256]*1.0))*log2((posLetter[i][256]*1.0)/(posLetter[i][j]*1.0));
		if(entropy > 0){
			cont ++;
			cout << "Entropy Line " << i << ": " << entropy << endl;
			avg_en +=  entropy;
		}
	}

	cout << endl << "Average Entropy: " <<  avg_en/cont << endl;
	return 0;
}

