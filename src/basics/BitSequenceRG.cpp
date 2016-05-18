/* BitSequenceRG.cpp
   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   New RANK, SELECT, SELECT-NEXT and SPARSE RANK implementations.
   Addaptation to libcds by Francisco Claude

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

#include "./BitSequenceRG.h"
#include <cassert>
#include <cmath>

const cds_uint W = sizeof(cds_uint) * 8;
const uint mask31 = 0x0000001F; // mask for obtaining the first 5 bits 


/////////////
//Rank(B,i)//
/////////////
//_factor = 0  => s=W*lgn
//_factor = P  => s=W*P
//Is interesting to notice
//factor=2 => overhead 50%
//factor=3 => overhead 33%
//factor=4 => overhead 25%
//factor=20=> overhead 5%

BitSequenceRG::BitSequenceRG() {
	data = NULL;
	n = 0;
	factor = 0;
}


BitSequenceRG::BitSequenceRG(cds_uint *bitarray, cds_word _n, cds_word _factor) {
	if(_factor == 0) 
		exit(-1);
	data=new cds_uint[_n / W + 1];
	for(cds_word i = 0; i < uint_len(_n, 1); i++) 
		data[i] = bitarray[i];
	for(cds_word i = uint_len(_n, 1); i < _n / W + 1; i++)
		data[i] = 0;
	n = _n;
	cds_uint lgn = bits(n - 1);
	if(_factor == 0) 
		factor = lgn;
	else 
		factor = _factor;

	b = 32;
	s = b * factor;
	integers = n / W + 1;
	BuildRank();
	length = n;
	ones = rank1(n-1);
}

BitSequenceRG::~BitSequenceRG() {
	delete [] Rs;
	delete [] data;
}

//Metodo que realiza la busqueda d
void BitSequenceRG::BuildRank() {
	cds_word num_sblock = n / s;
	// +1 pues sumo la pos cero
	Rs = new cds_uint[num_sblock + 5];
	for(cds_uint i = 0; i < num_sblock + 5; i++)
		Rs[i] = 0;
	cds_word j;
	Rs[0] = 0;
	for (j = 1; j <= num_sblock; j++){
		Rs[j] = Rs[j - 1];
		Rs[j] += BuildRankSub((j - 1) * factor, factor);
	}
}

cds_word BitSequenceRG::BuildRankSub(cds_word ini, cds_word bloques){
	cds_word rank = 0, aux;
	for(cds_uint i = ini; i < ini + bloques; i++){
		if(i < integers) {
			aux = data[i];
			rank += popcount((cds_word)aux);
		}
	}
	return rank;             //retorna el numero de 1's del intervalo
}

cds_word BitSequenceRG::rank1(const cds_word i1) const{
	cds_uint i = (cds_uint)i1;
	++i;
	cds_word resp = Rs[i / s];
	cds_uint aux = (i / s) * factor;
	for(cds_uint a = aux; a < i / W; a++)
		resp += popcount((cds_word)(data[a]));
	resp += popcount((cds_word)(data[i / W] & ((1 << (i & mask31)) - 1)));
	return resp;
}

bool BitSequenceRG::access(const cds_word i) const{
	return (1u << (i % W)) & data[i / W];
}

void BitSequenceRG::Save(ofstream & fp){
	  cds_uint var[2];
		var[0] = n;
		var[1] = factor;
		SaveValue(fp, var, 2);
		SaveValue(fp, data, integers);
		SaveValue(fp, Rs, n / s + 1);
}

BitSequenceRG * BitSequenceRG::Load(ifstream & fp) {
	assert(fp.good());
	BitSequenceRG * ret = new BitSequenceRG();
	cds_uint *var = LoadValue<cds_uint>(fp, 2);
	ret->n = var[0];
	ret->b = 32;
	ret->factor = var[1];
	ret->s = ret->b * ret->factor;
	ret->integers = (ret->n + 1) / W + ((ret->n + 1) % W != 0? 1:0);
	ret->data = LoadValue<cds_uint>(fp, ret->integers);
	ret->Rs = LoadValue<cds_uint>(fp, ret->n / ret->s + 1);
	ret->length = ret->n;
	ret->ones = ret->rank1(ret->n - 1);
	delete [] var;
	return ret;
}

cds_word BitSequenceRG::SpaceRequirementInBits() const{
	return uint_len(n, 1) * W + (n / s) * W + sizeof(this) * 8;
}

cds_word BitSequenceRG::getSize() const{
	return SpaceRequirementInBits()/8;
}

cds_word BitSequenceRG::SpaceRequirement() const{
	return n / 8 + (n / s) * sizeof(cds_uint) + sizeof(BitSequenceRG);
}

cds_word BitSequenceRG::select1(const cds_word x1) const{
	cds_uint x = x1;
	// returns i such that x=rank(i) && rank(i-1)<x or n if that i not exist
  // first binary search over first level rank structure
  // then sequential search using popcount over a int
  // then sequential search using popcount over a char
  // then sequential search bit a bit
	if(x > ones) 
		return (cds_word)(-1);
	//binary search over first level rank structure
	cds_uint l = 0, r = n / s;
	cds_uint mid = (l + r) / 2;
	cds_uint rankmid = Rs[mid];
	while(l <= r) {
		if (rankmid < x)
			l = mid + 1;
		else
			r = mid - 1;
		mid = (l + r) / 2;
		rankmid = Rs[mid];
	}
	//sequential search using popcount over a int
	cds_uint left = mid * factor;
	x -= rankmid;
	cds_uint j = data[left];
	cds_uint _ones = (cds_uint)popcount(j);
	while (_ones < x) {
		x -= _ones;
		left++;
		if(left > integers) 
			return n;
		j = data[left];
		_ones = (cds_uint)popcount(j);
	}
	//sequential search using popcount over a char
	left = left * b;
	rankmid = popcount8(j);
	if(rankmid < x) {
		j = j >> 8;
		x -= rankmid;
		left += 8;
		rankmid = popcount8(j);
		if(rankmid < x){
			j = j >> 8;
			x -= rankmid;
			left += 8;
			rankmid = popcount8(j);
			if (rankmid < x){
				j = j >> 8;
				x -= rankmid;
				left += 8;
			}
		}
	}
	// then sequential search bit a bit
	while(x > 0){
		if(j & 1) 
			x--;
		j = j >> 1;
		left++;
	}
	return left - 1;
}

cds_word BitSequenceRG::select0(const cds_word x1) const{
	uint x = (cds_uint)x1;
	// returns i such that x=rank_0(i) && rank_0(i-1)<x or n if that i not exist
	// first binary search over first level rank structure
	// then sequential search using popcount over a int
	// then sequential search using popcount over a char
	// then sequential search bit a bit
	if( x> n - ones) 
		return (cds_word)(-1);

	//binary search over first level rank structure
	if(x == 0) 
		return 0;
	cds_uint l = 0, r = n / s;
	cds_uint mid = (l + r) / 2;
	cds_uint rankmid = mid * factor * W - Rs[mid];
	while(l <= r){
		if(rankmid < x)
			l = mid + 1;
		else
			r = mid - 1;
		mid = (l + r) / 2;
		rankmid = mid * factor * W - Rs[mid];
	}
	//sequential search using popcount over a int
	cds_uint left = mid * factor;
	x -= rankmid;
	uint j = data[left];
	uint zeros = W - (cds_uint)popcount(j);
	while(zeros < x) {
		x -= zeros;
		left++;
		if(left > integers) 
			return n;
		j = data[left];
		zeros = W - (cds_uint)popcount(j);
	}
	//sequential search using popcount over a char
	left = left * b;
	rankmid = 8 - popcount8(j);
	if(rankmid < x) {
		j = j >> 8;
		x -= rankmid;
		left += 8;
		rankmid = 8 - popcount8(j);
		if(rankmid < x){
			j = j >> 8;
			x -= rankmid;
			left += 8;
			rankmid = 8 - popcount8(j);
			if(rankmid < x){
				j = j >> 8;
				x -= rankmid;
				left += 8;
			}
		}
	}
	// then sequential search bit a bit
	while(x > 0){
		if(j % 2 == 0 ) 
			x--;
		j = j >> 1;
		left++;
	}
	left--;
	if(left > n)  
		return n;
	else 
		return left;
}

