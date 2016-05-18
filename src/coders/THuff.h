/* THuff.h
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

#ifndef HUFF_H
#define HUFF_H

#include "./../basics/libcds.h"
#include "./../basics/io.h"

using namespace cds;
using namespace cds::basic;


class THuff{
	public:
		THuff();
		/** Creates Huffman encoder given symbols 0..lim with frequencies
		 * freq[i], ready for compression
		 *
		 * @author Gonzalo Navarro
		 */
		THuff(cds_uint *freq, cds_uint _lim);
		 
		/** Encodes symb, over stream[ptr...lim] (ptr and lim are
		 *  bit positions of stream). Returns the new ptr.
		 *
		 *  @author Gonzalo Navarro
		 */
		virtual cds_word encodeHuff(cds_uint symb, cds_word *stream, cds_word ptr);

		/** Decodes *symb, over stream[ptr...lim] (ptr and lim are
		 * bit positions of stream). Returns the new ptr.
		 *
		 * @author Gonzalo Navarro
		 */
		virtual cds_word decodeHuff(cds_uint *symb, cds_word *stream, cds_word ptr);


		/** Size of Huffman encoder written on file
		 * *
		 * *  @author Gonzalo Navarro
		 * */
		virtual cds_word getSize();

		/** Frees Huffman encoder
		 * *
		 * *  @author Gonzalo Navarro
		 * */
		virtual ~THuff();

		virtual cds_uint getDepth();

		/*Save only the data needed to decode information previously Huffman coded*/
		virtual void Save(ofstream & f);

		/*Load only the data needed to decode information previously Huffman coded*/
		static THuff * Load(ifstream &f);


	protected:
		cds_uint max,lim;
		cds_uint depth;        // max symbol length
		cds_uint *spos;            //symbol positions after sorting by decr freq (enc)
		cds_uint *symb;           // symbols sorted by freq (dec)
		cds_uint *num_dec;          // first pos of each length (dec)
		cds_uint *num_enc;          // number of each length (enc)
		cds_uint *fst;              // first code (numeric) of each length (dec)
		cds_ulong total;            // total length to achieve, in bits
};

#endif

