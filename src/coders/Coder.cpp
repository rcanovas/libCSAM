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

#include "Coder.h"

cds_word encodeVbyte(cds_word symb, cds_word *stream, cds_word pos){ 
	cds_word last_pos = pos, newvalue = 0, x = 0;
	if(symb == 0){
		SetVarField(stream, last_pos, last_pos + BASE_BITS, newvalue);
		last_pos += BASE_BITS + 1;
	}
	else{
		newvalue = symb;
		while(newvalue > 0){
				x = newvalue;
				newvalue = newvalue >> BASE_BITS;
				if(newvalue > 0)
					BitSet(&x, BASE_BITS, true);
				else
					BitSet(&x, BASE_BITS, false);
				SetVarField(stream, last_pos, last_pos + BASE_BITS, x);
				last_pos += BASE_BITS + 1;
		}
	}
	return last_pos;
}

cds_word decodeVbyte(cds_word *symb, cds_word *stream, cds_word pos){
	cds_word partialSum = 0, readblock = 0, mult = BASE_BITS;
	cds_word last_pos = pos, nextblock = 0;
	readblock = GetVarField(stream, last_pos, last_pos + BASE_BITS);
	partialSum = readblock % BASE;
	last_pos += BASE_BITS + 1;
	nextblock = BitGet(&readblock, BASE_BITS);
	while(nextblock){
		readblock = GetVarField(stream, last_pos, last_pos + BASE_BITS);
		partialSum += (readblock % BASE) << mult;
		last_pos += BASE_BITS + 1;
		mult += BASE_BITS;
		nextblock = BitGet(&readblock, BASE_BITS);
	}
	*symb = partialSum;
	return last_pos;
}

//symb > 0 always and number lower than 64!!
cds_word encodeUnary(cds_word symb, cds_word *stream, cds_word pos){
	cds_word zeros = 0;  //set all bits in 0
	cds_word newPos = 0;
	if(symb == 1){
		BitOne(stream, pos);
		newPos = pos + 1;
	}
	else{
		SetVarField(stream, pos, (pos + symb - 2), zeros);
		BitOne(stream, pos + symb - 1);
		newPos = pos + symb;
	}
	return newPos;
}

cds_word decodeUnary(cds_word *symb, cds_word *stream, cds_word pos){
	//this should be improve using pop-count or just implementing a 
	//method to fast look for the next 1
	cds_word newPos = pos;
	cds_word x = 1;
	while(!BitGet(stream, newPos)){
		x++;
		newPos++;
	}
	*symb = x;
	return newPos + 1;
}

//symb > 0 always!!
cds_word encodeGamma(cds_word symb, cds_word *stream, cds_word pos){
	cds_word lenbit, x;
	cds_word newPos = pos;
	if(symb == 1){
		BitOne(stream, pos);
		newPos++;
	}
	else{
		lenbit = (cds_word)(floor(log2(symb)));
		newPos = encodeUnary(lenbit + 1, stream, pos);
		x = GetVarField(&symb, 0, lenbit-1);
		SetVarField(stream, newPos, (newPos + lenbit-1), x);
		newPos += lenbit;
	}
	return newPos;
}

cds_word decodeGamma(cds_word *symb, cds_word *stream, cds_word pos){
	cds_word lenbit = 0;
	cds_word newPos = 0;
	cds_word x;
	if(BitGet(stream, pos)){
		*symb = 1;
		newPos = pos + 1;
	}
	else{
		newPos = decodeUnary(&lenbit, stream, pos); 
		x =  GetVarField(stream, newPos, newPos + lenbit - 2);
		*symb = (1<<(lenbit-1)) + x ;
		newPos +=  (lenbit - 1);
	}
	return newPos;
}

//symb >= 0
cds_word encodeGolomb(cds_word symb, cds_word *stream, cds_word pos, cds_word m){
	cds_word newPos = pos;
	cds_word q = symb/ m, lim = 0;
	//write q in unary
	for(cds_word i =  0; i < q; i++)
		BitZero(stream, newPos + i);
	BitOne(stream, newPos + q);
	newPos += q + 1;
	cds_word r = (symb % m);
	//write r in binary
	cds_word b = (cds_word)(ceil(log2(m)));
	if((m & (m - 1)) == 0){ //m is a power of two
		SetVarField(stream, newPos, (newPos + b - 1), r);
		newPos += b;
	}
	else{
		lim = (1 << b) - m; // 2^b - m
		if(r < lim) 
			newPos = saveInverseValue(r, stream, newPos, b -1);   
		else
			newPos = saveInverseValue(r + lim, stream, newPos, b);
	}
	return newPos;
}


cds_word decodeGolomb(cds_word *symb, cds_word *stream, cds_word pos, cds_word m){
	cds_word newPos = pos;
	cds_word q = 0, r = 0, lim = 0;
	cds_word b = (cds_word)(ceil(log2(m)));
	//get q
	while(!BitGet(stream, newPos + q))
		q += 1;
	newPos += q + 1;
	//get r
	if((m & (m - 1)) == 0){ //m is a power of two
		r = GetVarField(stream, newPos, newPos + b - 1);
		newPos += b;
	}
	else{ 
		lim = (1 << b) - m; // 2^b - m
		newPos = readInverseValue(&r, stream, newPos, b - 1);
		if((r + 1) > lim){
			r = (r << 1) + BitGet(stream, newPos) - lim;
			newPos++;
		}
		else
			b = b - 1;
	}
	*symb = q  * m + r;
	return newPos;
}

cds_word saveInverseValue(cds_word symb, cds_word *stream, cds_word pos, cds_word b){
	cds_word newPos = pos;
	for(cds_word i = b; i > 0; i --){
		if(BitGet(&symb, i - 1))
			BitOne(stream, newPos);
		else
			BitZero(stream, newPos);
		newPos ++;
	}
	return newPos;
}

cds_word readInverseValue(cds_word *symb, cds_word *stream, cds_word pos, cds_word b){
	cds_word newPos = pos, value = 0;
	for(cds_word i = 0; i < b; i ++){
		value = value << 1;
		value += BitGet(stream, newPos);
		newPos ++;
	}
	*symb = value;
	return newPos;
}

void computeEntropy(string array, cds_word size){
	cds_word alphabet[256];
	double entropy =0.0;
	for(int w = 0; w < 256; w++)
		alphabet[w]=0;
	for(cds_word i = 0; i < size; i++)
		alphabet[(int)(array[i])] += 1;
	for(int w = 0; w < 256; w++){
		if(alphabet[w] > 0)
			entropy += (alphabet[w]*1.0/size)*log2(size*1.0/alphabet[w]);
	}
	cout << "Length: " << size << "  Entropy: " <<  entropy <<  "  Theoretical Space use in MB: "  <<  (entropy * (size) / 8388608) << endl;
}




