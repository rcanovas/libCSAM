/* THuff.cpp
   Copyright (C) 2008, Gonzalo Navarro, all rights reserved.
	 							 2012, Rodrigo Canovas
   Canonical Huffman

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
// implements canonical Huffman

#include "THuff.h"

typedef struct{ 
	cds_uint freq;
  cds_uint symb;
	union{
		int prev;
		cds_uint depth;
	} h;
	int ch1,ch2;
} Ttree;

static void sort(Ttree *tree, int lo, int up) {
	cds_uint i, j;
	Ttree temp;
	while(up > lo){
		i = lo;
		j = up;
		temp = tree[lo];
		while(i < j){
			while(tree[j].freq > temp.freq) 
				j --;
			tree[i] = tree[j];
			while(i < j && tree[i].freq <= temp.freq) 
				i++;
			tree[j] = tree[i];
		}
		tree[i] = temp;
		if(i - lo < up - i){ 
			sort(tree, lo, i - 1); 
			lo = i + 1; 
		}
		else{ 
			sort(tree, i + 1, up); 
			up = i - 1; 
		}
	}
}

static void setdepths(Ttree *tree, cds_uint node, int depth){
	// leaf
	if(tree[node].ch1 == -1){
		tree[node].h.depth = depth;
		return;
	}
	setdepths(tree, tree[node].ch1, depth + 1);
	setdepths(tree, tree[node].ch2, depth + 1);
}

THuff::THuff(cds_uint *freq, cds_uint _lim){
	int i = 0, j, d;
	Ttree *tree;
	uint ptr,last,fre;
	// remove zero frequencies
	max = _lim;
	tree = new Ttree[2 * (_lim + 1) - 1];
	j = 0;
	for(i = 0; i <= (int)_lim; i++){
		if(freq[i] > 0){
			tree[j].freq = freq[i];
			tree[j].symb = i;
			j++;
		}
	}
	lim = _lim = j - 1;
	// now run Huffman algorithm
	sort (tree, 0, _lim);
	for(i = 0; i <= (int)_lim; i++){
			tree[i].h.prev = i + 1;
			tree[i].ch1 = tree[i].ch2 = -1;
	}
	tree[_lim].h.prev = -1;
	// last = next node to process, ptr = search point, fre = next free cell
	// leaves are in 0..lim in decreasing freq order
	// internal nodes are in lim+1.. 2*lim, created in incr. fre order
	last=0; ptr = 0; fre = _lim + 1;
	for(i = 0; i < (int)_lim; i++){
		tree[fre].ch1 = last;
		last = tree[last].h.prev;
		tree[fre].ch2 = last;
		tree[fre].freq = tree[tree[fre].ch1].freq+tree[tree[fre].ch2].freq;
		while((tree[ptr].h.prev != -1) && (tree[tree[ptr].h.prev].freq <= tree[fre].freq))
			ptr = tree[ptr].h.prev;
		tree[fre].h.prev = tree[ptr].h.prev;
		tree[ptr].h.prev = fre;
		last = tree[last].h.prev;
		fre ++;
	}
	// now assign depths recursively
	setdepths(tree, 2 * _lim, 0);
	spos = new cds_uint[max + 1];
	for(i = 0; i <= (int)max; i++) 
		spos[i] = ~0;
	num_enc = new cds_uint[_lim+ 1 ]; // max possible depth
	d=0;
	for(i = _lim; i >= 0; i--){
		spos[tree[i].symb] = (cds_uint)i;
		while((int)tree[i].h.depth > d){
			num_enc[d] = i + 1; 
			d++; 
		}
	}
	num_enc[d] = 0;
	depth = d;
	for(d = depth; d > 0; d--) 
		num_enc[d] = num_enc[d - 1] - num_enc[d];
	num_enc[0] = (_lim == 0);
	cds_uint * Htmp = new cds_uint[depth + 1];
	for(cds_uint i = 0; i < depth + 1; i++)
		Htmp[i] = 0;
	for(cds_uint i = 0; i < depth + 1; i++)
		Htmp[i] = num_enc[i];
	delete [] num_enc;
	num_enc = Htmp;
	total = 0;
	for(i = 0; i <= (int)_lim; i++)
		total += freq[tree[i].symb] * tree[i].h.depth;
	delete [] tree;

	/**for decode**/
	symb = new uint[lim + 1];
	cds_uint aux = 0;
	for(i = 0; i < (int)(lim + 1); i++)
		symb[i] = 0;
	for(i = 0; i <= (int)max; i++){
		aux = spos[i];
		if(aux != (cds_uint)~0) 
			symb[aux] = i;
	}
	cds_uint dold, dact;
	num_dec = new cds_uint[depth + 1];
	fst = new cds_uint[depth + 1];
	fst[depth] = 0; 
	dold = 0;
	for(d = depth - 1; d >= 0; d--){
		dact = num_enc[d + 1];
		fst[d] = (fst[d + 1] + dact) >> 1;
		num_dec[d + 1] = dold;
		dold += dact;
		if(d==0)
			break;
	}
	num_dec[0] = dold;	
}

void bitzero(register cds_word *e, register cds_word p, register cds_word len){
	e += p/kWordSize; 
	p %= kWordSize;
	if(p+len >= kWordSize){
		*e &= ~((1 << p) - 1);
		len -= p;
		e++; 
		p = 0;
	}
	while(len >= kWordSize){
		*e++ = 0;
		len -= kWordSize;
	}
	if(len > 0)
		*e &= ~(((1 << len) - 1) << p);
}

cds_word THuff::encodeHuff(cds_uint sym, cds_word *stream, cds_word ptr){
	cds_word pos, code, d;
	pos = spos[sym];
	code = 0;
	d = depth;
	while(pos >= num_enc[d]){
		code = (code + num_enc[d]) >> 1;
		pos -= num_enc[d--];
	}
	code += pos;
	if(d > kWordSize){
		bitzero(stream, ptr, d - kWordSize); 
		ptr += d - kWordSize;
		d = kWordSize; 
	}
	while(d--){
		if((code >> d) & 1)
			BitOne(stream,ptr);
		else 
			BitZero(stream,ptr);
		ptr++;
	}
	return ptr;
}

cds_word THuff::decodeHuff(cds_uint *sym, cds_word *stream, cds_word ptr){
	cds_word pos, d;
	pos = 0;
	d = 0;
	while(pos < fst[d]){
		pos = (pos << 1) | BitGet(stream,ptr);
		ptr++; 
		d++;
	}
	*sym = symb[num_dec[d] + pos - fst[d]];
	return ptr;
}

cds_word THuff::getSize(){
	uint mem = sizeof(THuff);
	mem += 2 *(depth + 1) * sizeof(cds_uint); //num_dec + fst 
	mem += (lim + 1) * sizeof(cds_uint);    //symb
	return mem;
}

cds_uint THuff::getDepth(){
	return depth;
}

THuff::~THuff(){
	delete [] spos;
	delete [] symb;
	delete [] fst;
	delete [] num_enc;
	delete [] num_dec;
}

THuff::THuff(){
	spos = NULL;
	symb = NULL;
	fst = NULL;
	num_enc = NULL;
	num_dec = NULL;
}

void THuff::Save(ofstream & fp){
	cds_uint var[3];
	var[0] = max;
	var[1] = lim;
	var[2] = depth;
	SaveValue(fp, var, 3);
	SaveValue(fp, symb, lim + 1); 
	SaveValue(fp, num_dec, depth + 1);
//	SaveValue(fp, num_enc, depth + 1);
	SaveValue(fp, fst, depth + 1);
	SaveValue(fp, total);
}

THuff * THuff::Load(ifstream &fp){
	THuff *th = new THuff();
	cds_uint *var = LoadValue<cds_uint>(fp, 3);
	th->max = var[0];
	th->lim = var[1];
	th->depth = var[2];
	th->symb = LoadValue<cds_uint>(fp, th->lim + 1);
	/*th->spos = new uint[th->max+1];
	for(cds_uint i=0;i<=H.max;i++) 
		H.s.spos[i] = (uint)~0;
	for (i=0;i<=H.lim;i++) 
		H.s.spos[symb[i]] = i;*/
	th->num_dec = LoadValue<cds_uint>(fp, th->depth + 1);
	//th->num_enc = LoadValue<cds_uint>(fp, th->depth + 1);
	th->fst = LoadValue<cds_uint>(fp, th->depth + 1);
	th->total = LoadValue<cds_ulong>(fp);
	delete [] var;
	return th;
}




