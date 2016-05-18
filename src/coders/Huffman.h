/* Huffman.h
   Copyright (C) 2008, Francisco Claude, all rights reserved.
	 							 2012, Rodrigo Canovas
   Wrapper for huff written by Gonzalo Navarro

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef HUFFMAN_H
#define HUFFMAN_H

#include <cmath>
#include <algorithm>
#include "THuff.h"

using namespace std;


class Huffman{

	public:
		/** Creates the codes for the sequence seq of length n.
		 * Note that this only create the codification tree, this
		 * program do not compress the array.
		 * In case that frec is true it means that seq contain the
		 * frequency of each value*/
		Huffman(cds_uint * seq, cds_word n, bool frec = false);

		 /** Creates the codes for the vector seq.
			* Note that this only create the codification tree, 
			* this program do not compress the array. In case 
			* that frec is true it means that the vector seq 
			* contain the frequency of each value*/
		Huffman(vector<int> seq, bool frec = false);



		
		/** Creates the codes for the chars in the string seq**/
		Huffman(string seq);

		virtual ~Huffman();

		/** Encodes symb into stream at bit-position pos,
		 * returns the ending position (bits) */
		virtual cds_word encode(cds_uint symb, cds_word *stream, cds_word pos);

		/** decodes into symb from stream at bit-position
		 * pos, returns the new position.*/
		virtual cds_word decode(cds_uint * symb, cds_word *stream, cds_word pos);

		/** Returns the maximum length of a code */
		virtual cds_word maxLength();

		/** Returns the size of the table */
		virtual cds_word getSize();

		 /*Save only the data needed to decode information previously Huffman coded*/
		virtual void Save(ofstream & f);

		/*Load only the data needed to decode information previously Huffman coded*/
		static Huffman * Load(ifstream &f);
	
	protected:
		Huffman();
		THuff *huff_table;
};

#endif
