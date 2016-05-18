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

#ifndef _CFIELDS_H
#define _CFIELDS_H

#include <cmath>
#include <algorithm>
#include "./../basics/extras.h"
#include "./../basics/buffer.h"

struct CigOp{
	char op;
	int value;
};

struct FieldBlock{
	vector<string> QNAME;
	vector<string> FLAG;
	vector<string> MAPQ;
	vector<string> CIGAR;
	vector<string> RNEXT;
	vector<string> PNEXT;
	vector<string> TLEN;
	vector<string> OTHER;
};


class CFields{

	public:
		CFields(cds_word sample);
		virtual ~CFields();
		virtual void ProcessLine(string qname, string flag, string mapq, string cigar, string rnext, string pnext, 
														 string tlen, string others, ofstream &fp, cds_word lineNumber);
		virtual void CreateFile(ofstream &fp);
		virtual void ExtractBlock(ifstream &fp, cds_word blockNum, FieldBlock *fields, char **buff, size_t *init_buff, size_t MaxBytes = MAXRAMBYTESDECO);
		virtual void ExtractSelection(ifstream &fp, cds_word blockNum, FieldBlock *fields, char **buff, size_t *init_buff, int selection, size_t MaxBytes = MAXRAMBYTESDECO);
		virtual vector<CigOp> CigarToArray(string cigar);
		static CFields *Load(ifstream &fp);

		virtual void ClearFieldBlock(FieldBlock *fields);

		cds_word rate;



	protected:
		int no_other;
		size_t init_index;
		size_t end_field;
		string Cqname, Cflag, Cmapq, Ccigar, Crnext, Cpnext, Ctlen, Cother;
		vector <size_t> IndexBlock;
		vector <cds_word> IndexFlag;
		vector <cds_word> IndexMapq;
		vector <cds_word> IndexQname;
		vector <cds_word> IndexRnext;
		vector <cds_word> IndexPnext;
		vector <cds_word> IndexTlen;
		vector <cds_word> IndexOthers;
		

	private:
		void SaveBlock(string *qname, string *flag, string *mapq, string *cigar, string *rnext, string *pnext, 
				           string *tlen, string *others, ofstream &fp);


		//stats
	  cds_word memCigar, memFlag, memMapq, memQname, memRnext, memPnext, memTlen, memOther;
};


#endif
