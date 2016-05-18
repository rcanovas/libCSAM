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
#ifndef _CRPSPREMF_H
#define _CRPSPREMF_H

#include "./CRePoSe.h"

class CRPSPreMF: public CRePoSe{

	public:
		CRPSPreMF(cds_word sample, bool create = true);
		virtual ~CRPSPreMF();
		virtual void CreateFile(ifstream &SamFile, ofstream &fp, int SamOrRps, vector<cds_word>* IndexClusterToLines = NULL);
		virtual void DecompressSequence(ifstream &input, ofstream &output);
		virtual cds_word GetInterval(ifstream &input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, size_t MaxBytes = MAXRAMBYTESDECO);
	  //operations for .csam files	
		virtual void DecompressAll(ifstream **input, ofstream &output, CQualBlock *cqual, CFields *cfields);
		virtual cds_word GetIntervalAll(ifstream **input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, CQualBlock *cqual, CFields *cfields, vector<cds_word>* ClusToLine, int selection, size_t MaxBytes = MAXRAMBYTESDECO);
		virtual cds_word CountReads(ifstream **input, string rname, cds_word pos_x, cds_word pos_y, int strand, CFields *cfields, vector<cds_word>* ClusToLine, size_t MaxBytes = MAXRAMBYTESDECO);


		static CRPSPreMF *Load(ifstream &fp);



	protected:
		void ComputeCSeq(ifstream &SamFile, ofstream &fp, int SamOrRps);

	private:
		cds_word lastFieldBlock;
		FieldBlock fieldblock;
		vector<string> vqual;

};

#endif
