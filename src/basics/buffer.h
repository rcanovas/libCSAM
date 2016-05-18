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

#ifndef SRC_BASIC_BUFFER_H_
#define SRC_BASIC_BUFFER_H_

#include <cmath>
#include <algorithm>
#include "./libcds.h"
#include "./io.h"



using namespace cds;
using namespace cds::basic;
using namespace std;


const size_t MAXRAM = 2147483648;    //in bits
//2147483648; // aprox 250 MB     
//83886080; //10 mb
//8388608;  //1 mb 
const size_t MAXRAMBYTES = 268435456;
//268435456; 
//10485760; //10 mb
//1048576; //1 mb 

const size_t MAXRAMBYTESDECO = 1048576; 
//209715200; //200 mb
//104857600; //100 mb
//10485760; //10 mb
//1048576; //1 mb
//200000; //200kb

/*Check if tmp variable is using more than MAXRAM.
 * If that is the case, write the string to fp and reset 
 * the string to empty*/
void check_tmp_string(string *tmp, ofstream &fp);

/*Check if buffer_use is greater or equal than MAXRAM. If that 
 *is the case, write MAXRAM bits from the buffer to fp and 
 *modify the buffer to only keep only the extra data and update 
 *the information on buffer_use*/
void check_buffer(cds_word *buffer, cds_word *buffer_use, ofstream &fp);


cds_word * init_buffer();

/*compute entropy of an array*/
double ComputeEntropy(cds_word *array, cds_word length);

double ComputeEntropy(cds_uint *array, cds_word length);

/*separate string by the delimiter storing the resulting string into tokens*/
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);

#endif




