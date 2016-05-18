/* BitSequenceRG.h
   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   New RANK, SELECT, SELECT-NEXT and SPARSE RANK implementations.
   Adaptation to libcds by Francisco Claude

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

#ifndef _STATIC_BITSEQUENCE_BRW32_H
#define _STATIC_BITSEQUENCE_BRW32_H

#include "./libcds.h"
#include "./io.h"

using namespace cds;
using namespace cds::basic;



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

/** Implementation of Rodrigo Gonzalez et al. practical rank/select solution [1].
*  The interface was adapted. Work well when cds_uint use 32 bits
*
*  [1] Rodrigo Gonzalez, Szymon Grabowski, Veli Makinen, and Gonzalo Navarro.
*      Practical Implementation of Rank and Select Queries. WEA05.
*
*  @author Rodrigo Gonzalez
*/
class BitSequenceRG{
	private:
		/** Length of the bitstring */
		cds_word length;
		/** Number of ones in the bitstring */
		cds_word ones;
		cds_word n,integers;

		cds_word factor,b,s;
		cds_uint *Rs;            //superblock array

		//uso interno para contruir el indice rank
		cds_word BuildRankSub(cds_word ini, cds_word fin);
		void BuildRank();    //crea indice para rank
		BitSequenceRG();
		cds_word SpaceRequirementInBits() const;
		cds_word SpaceRequirement() const;

	public:
		cds_uint *data;
		/** Build the BitSequenceRG with a sampling factor <code>factor</code>
		 * The <code>factor</code> value has to be either 2,3,4 or 20, being the first one the fastest/bigger.
		 * */
		BitSequenceRG(cds_uint *bitarray, cds_word _n, cds_word _factor);
            
		/** Build the BitSequenceRG with a sampling factor <code>factor</code>
		 * The <code>factor</code> value has to be either 2,3,4 or 20, being the first one the fastest/bigger.
		 * */
           
		~BitSequenceRG();    //destructor
          
		virtual bool access(const cds_word i) const;
				
		//Nivel 1 bin, nivel 2 sec-pop y nivel 3 sec-bit
		virtual cds_word rank1(const cds_word i) const;

		// gives the position of the x:th 1.
		virtual cds_word select0(cds_word x) const;
 
		// gives the position of the x:th 1.
		virtual cds_word select1(cds_word x) const;
        
		virtual cds_word getSize() const;

		/*load-save functions*/
        
		virtual void Save(ofstream & f);
        
		static BitSequenceRG * Load(ifstream & f);
   	
};


#endif
