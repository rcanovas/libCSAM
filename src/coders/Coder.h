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
#ifndef CODER_H
#define CODER_H


#include "./../basics/libcds.h"
#include "./../basics/io.h"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace cds::basic;


#define BASE_BITS 4
#define BASE 16

/** Vbyte(4) encode symb into stream at bit-position, 
 *  returns the ending positon (bits)*/
cds_word encodeVbyte(cds_word symb, cds_word *stream, cds_word pos);

/** Decode into symb from stream at bit-position
 *  pos, returns the new position. Each symbol was Vbyte(4) encoded*/
cds_word decodeVbyte(cds_word *symb, cds_word *stream, cds_word pos);

/** Unary encode symb into stream at bit-position pos,
 *  returns the ending position (bits) */
cds_word encodeUnary(cds_word symb, cds_word *stream, cds_word pos);

/** Decode into symb from stream at bit-position
 * pos, returns the new position. Each symbol was Unary encoded*/
cds_word decodeUnary(cds_word *symb, cds_word *stream, cds_word pos);


/** Elias Gamma encode symb into stream at bit-position pos,
 *  returns the ending position (bits) */
cds_word encodeGamma(cds_word symb, cds_word *stream, cds_word pos);

/** Decode into symb from stream at bit-position
 * pos, returns the new position. Each symbol was Elias Gamma encoded*/
cds_word decodeGamma(cds_word *symb, cds_word *stream, cds_word pos);


/** Golomb encode symb into stream at bit-position pos,
 *  returns the ending position (bits). Golomb is parameterized by m */       
cds_word encodeGolomb(cds_word symb, cds_word *stream, cds_word pos, cds_word m);
								

/** Decode into symb from stream at bit-position
 * pos, returns the new position. Each symbol was Golomb encode*/     
cds_word decodeGolomb(cds_word *symb, cds_word *stream, cds_word pos, cds_word m);

/*Save in the stream from pos, the first b bits of symb in inverse order*/
cds_word saveInverseValue(cds_word symb, cds_word *stream, cds_word pos, cds_word b);

/*Read from the stream from pos, the first b bits of symb in inverse order*/
cds_word readInverseValue(cds_word *symb, cds_word *stream, cds_word pos, cds_word b);


void computeEntropy(string array, cds_word size);


#endif
