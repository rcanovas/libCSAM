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
#ifndef _CQUALLL_H
#define _CQUALLL_H

#include <cmath>
#include <algorithm>
#include "./../basics/extras.h"
#include "CQualBlock.h"


class CQualLL : public CQualBlock{

	public:
		/*Create a CQualLL object*/
		CQualLL(cds_word sample);
		virtual ~CQualLL();
		virtual void ProcessLine(string qual);
		virtual void CreateFile(ofstream &fp);

		virtual void DecompressQual(ifstream &input, ofstream &output);
		virtual vector<string> GetInterval(ifstream &input, cds_word x, cds_word y, cds_word **buff, cds_word *init_buff, size_t MaxBytes = MAXRAMBYTESDECO);

		static CQualLL *Load(ifstream &fp);

	protected:
		string CQual;
};


#endif
