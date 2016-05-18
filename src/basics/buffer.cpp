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

#include "./buffer.h"

cds_word * init_buffer(){
	cds_word *buffer = new cds_word[(cds_word)( (1.2 * MAXRAM + kWordSize)  / kWordSize)];  //ask for 20% in case of overflow
	for(cds_word w = 0; w < (cds_word)((1.2 * MAXRAM + kWordSize)  / kWordSize); w++)
		buffer[w] = 0;
	return buffer;
}


void check_tmp_string(string *tmp, ofstream &fp){
	  if((8 * (*tmp).length()) >= MAXRAM){ //aprox 250 MB      
			fp << (*tmp);  //write to fileQual (it must exist)
			*tmp = "";
		}
}

void check_buffer(cds_word *buffer, cds_word *buffer_use, ofstream &fp){
	cds_word copy = 0, pos = 0;
	if(*buffer_use >= MAXRAM){
		copy = (MAXRAM + kWordSize - 1) / kWordSize;
		SaveValue(fp, buffer, copy);
		while((copy + pos) * kWordSize < *buffer_use){
			buffer[pos] = buffer[copy + pos];
			pos ++;
		}
		*buffer_use = *buffer_use - copy * kWordSize;
	}
}


double ComputeEntropy(cds_word *array, cds_word length){
	double entropy = 0;
	cds_word total = 0;
	for(cds_word i = 0; i < length; i++)
		total += array[i];
	cout << "Total: " << total << endl;
	for(cds_word i = 0; i < length; i++)
		if(array[i] > 0){	
	//		cout << "Value: " << i << "  Occ: " << array[i] << endl;
			entropy += ((array[i]*1.0)/(total*1.0))*log2((total*1.0)/(array[i]*1.0));
		}
	double bits =  total * entropy;
	cout << "Theoretical Space in MB: " << bits / 8388608 << endl; 
	//cout << endl;
	return entropy;
}


double ComputeEntropy(cds_uint *array, cds_word length){
	double entropy = 0;
	cds_word total = 0, count = 0;
	for(cds_word i = 0; i < length; i++)
		total += array[i];
	cout << "Total: " << total << endl;
	for(cds_word i = 0; i < length; i++)
		if(array[i] > 0){	
			count ++;
	//		cout << "Value: " << i << "  Occ: " << array[i] << endl;
			entropy += ((array[i]*1.0)/(total*1.0))*log2((total*1.0)/(array[i]*1.0));
		}
	cout << "count: " << count << endl;
	double bits =  total * entropy;
	cout << "Theoretical Space in MB: " << bits / 8388608 << endl; 
	bits = total * ceil(log2(count));
	cout << "Worse case Space in MB: " << bits / 8388608 << endl;
	//cout << endl;
	return entropy;
}



void Tokenize(const string& str, vector<string>& tokens, const string& delimiters){
	string::size_type last_pos = 0;
	string::size_type pos = str.find_first_of(delimiters);
	while(pos!=string::npos) {
		tokens.push_back(str.substr(last_pos,pos-last_pos));
		last_pos = pos+1;
		if(last_pos >= str.length())
			break;
		pos = str.find_first_of(delimiters,pos+1);
	}
	if(last_pos<str.length())
		tokens.push_back(str.substr(last_pos));
}


