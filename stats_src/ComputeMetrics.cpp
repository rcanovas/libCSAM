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
	if(argc ==1){
		cout << "Use: ./ComputeMetrcis <ref_arch> <arch>" << endl;
		return 0;
	}

	string ref_filename = argv[1];
	string filename = argv[2];


	/*get the Quality fields*/
	string line, qual;
	unsigned long cont = 0;
	
	double manhattan = 0, max_min = 0, mse =0;
	double chebyshev = 0, soergel = 0, lorentzian = 0;


	double local_max_min, aux_max_min, local_mse;
	double sum_x_y, aux_x_y, max_x_y, local_lor, sum_max;

	ifstream RefFile(ref_filename.c_str());
	if (!RefFile.is_open()){
		cout << "Unable to open file" << endl;
		return 1;
	}
	ifstream LossyFile(filename.c_str());
	if (!LossyFile.is_open()){
		cout << "Unable to open file" << endl;
		return 1;
	}

	getline(RefFile, line);
	getline(LossyFile, qual);
	cout << "length reads: " << line.length() << endl;
	while (RefFile.good() ){
		sum_x_y = 0, aux_max_min = 0, local_max_min = 0, local_mse = 0;
		max_x_y = 0, local_lor = 0, sum_max = 0;
		for(unsigned int i = 0; i < line.length(); i++){
			aux_x_y = abs((int)(line[i]) - (int)(qual[i]));
			sum_x_y += aux_x_y;
			aux_max_min = (max((int)(line[i]), (int)(qual[i])) - 33)* 1.0  / (min((int)(line[i]), (int)(qual[i])) - 33.0);
			local_max_min = max(local_max_min, aux_max_min);
			local_mse += (aux_x_y * aux_x_y);
			max_x_y = max(max_x_y, aux_x_y);
			sum_max += max((int)(line[i]), (int)(qual[i]));
			local_lor += log2(1 + aux_x_y);
		}
		manhattan += (sum_x_y / line.length()); 
		max_min += local_max_min;
		mse += (local_mse * 1.0 / line.length());
		chebyshev += max_x_y;
		soergel += (sum_x_y * 1.0 / sum_max);
		lorentzian += local_lor;
		cont ++;
		getline(RefFile, line);
		getline(LossyFile, qual);
	}
	
	while(LossyFile.good() ){
		cout << qual << endl;
		getline(LossyFile, qual);
	}
		
	RefFile.close();
	LossyFile.close();
	cout << "number of lines: " << cont << endl;
	cout << endl;
	cout << "Manhattan:  " << (manhattan / cont) << endl;
	cout << "Max:Min:  " << (max_min / cont) << endl;
	cout << "MSE:  " << (mse / cont) << endl;
	cout << "Chebyshev:  " << (chebyshev / cont) << endl;
	cout << "Soergel:  " << (soergel / cont) << endl;
	cout << "Lorentzian:  " << (lorentzian / cont) << endl;

	return 0;
}

