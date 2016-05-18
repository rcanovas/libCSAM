/* Huffman.cpp
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

#include "Huffman.h"

Huffman::Huffman(cds_uint * symb, cds_word n, bool frec){
	cds_uint max_v = 0;
	if(!frec){
		for(cds_word i = 0; i < n; i++)
			max_v = max(max_v, symb[i]);
		cds_uint *occ = new cds_uint[max_v + 1];
		for(cds_word i = 0; i < max_v + 1; i++)
			occ[i] = 0;
		for(cds_word i = 0; i < n; i++)
			occ[symb[i]]++;
		huff_table = new THuff(occ, max_v);
		delete [] occ;
	}
	else{
		max_v = n;
		huff_table = new THuff(symb, max_v);
	}
}


Huffman::Huffman(vector<int> seq, bool frec){
	vector<int>::iterator it;
	cds_uint max_v = 0;
	cds_uint *occ;
	if(!frec){
		for(it = seq.begin(); it != seq.end(); ++it)
			max_v = max(max_v, (cds_uint)(*it));
		occ = new cds_uint[max_v + 1];
		for(cds_word i = 0; i < max_v + 1; i++)
			occ[i] = 0;
		for(it = seq.begin(); it != seq.end(); ++it)
			occ[(*it)]++;
	}
	else{
		max_v = seq.size() - 1;
		occ = new cds_uint[max_v + 1];
		for(cds_word w = 0; w <= max_v; w++)
			occ[w] = seq[w];
	}
	huff_table = new THuff(occ, max_v);
	delete [] occ;
}


Huffman::Huffman(string seq){
	cds_uint max_v = 0;
	for(cds_word i = 0; i < seq.size(); i++)
		max_v = max(max_v, (cds_uint)seq[i]);
	cds_uint *occ = new cds_uint[max_v + 1];
	for(cds_word i = 0; i < max_v + 1; i++)
		occ[i] = 0;
	for(cds_word i = 0; i < seq.size(); i++)
		occ[(cds_uint)seq[i]]++;
	huff_table = new THuff(occ, max_v);
	delete [] occ;
}


Huffman::Huffman(){
	huff_table = NULL;
}

Huffman::~Huffman() {
	delete huff_table;
}

cds_word Huffman::maxLength(){
	return huff_table->getDepth();
}

cds_word Huffman::getSize(){
	return sizeof(Huffman) + huff_table->getSize();
}

cds_word Huffman::encode(cds_uint symb, cds_word *stream, cds_word pos){
	return huff_table->encodeHuff(symb, stream, pos);
}

cds_word Huffman::decode(cds_uint * symb, cds_word *stream, cds_word pos){
	return huff_table->decodeHuff(symb, stream, pos);
}

void Huffman::Save(ofstream & f){
	huff_table->Save(f);
}
            
    /*Load only the data needed to decode information previously Huffman coded*/            
Huffman * Huffman::Load(ifstream &f){
	Huffman *h = new Huffman();
	h->huff_table = THuff::Load(f);
	return h;
}



