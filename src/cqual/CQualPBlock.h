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
#ifndef _CQUALPBLOCK_H
#define _CQUALPBLOCK_H

#include <cmath>
#include <algorithm>
#include "CQualBlock.h"


class CQualPBlock : public CQualBlock{

	public:
		/*Create a CQualPblock object
		 * lossyparam: percentange of different accepted (ex: 2)
		 * m: 0- Quality values stored as plain ASCII bytes
		 * 		1- Quality values stored as binary value using range between min and 
		 * 		   max values over the whole file.
		 * 		2- Quality values stored as binary value using range between local min 
		 *			 and max values over each quality line*/
		CQualPBlock(int lossyparam, int m, cds_word sample);
		virtual ~CQualPBlock();
		virtual void ProcessLine(string qual);
		virtual void CreateFile(ofstream &fp);

		static CQualPBlock *Load(ifstream &fp);

	protected:

};


#endif
