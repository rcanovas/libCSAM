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

#ifndef _CQUALBLOCK_H
#define _CQUALBLOCK_H

#include "./../basics/buffer.h"
#include "./../coders/Coder.h"

using namespace cds;
using namespace cds::basic;
using namespace std;


const cds_word LOSSLESS = 1;
const cds_word PBLOCK = 2;
const cds_word RBLOCK = 3;

const cds_word MaxLineSizeByte = kBytesPerWord * 350;  //extreme estimation of max space used to store a line
const cds_word MaxLineSize = kWordSize * 350;

class CQualBlock{

	public:
		virtual ~CQualBlock();
		virtual void ProcessLine(string qual) = 0;
		virtual void CreateFile(ofstream &fp) = 0;
		virtual void DecompressQual(ifstream &input, ofstream &output);
		virtual vector<string> GetInterval(ifstream &input, cds_word x, cds_word y, cds_word **buff, cds_word *init_buff, size_t MaxBytes = MAXRAMBYTESDECO);
		virtual cds_word GetNumberOfLines();

		static CQualBlock *Load(ifstream &fp);



	protected:
		cds_word IndexRate;               //Every IndexRate lines we store a pointer to fast access 
		cds_word *IndexLine;
		cds_word numberOfValues;          //total numer of quality values stored
		cds_word numberOfLines;           //total number of quality lines stored
		int variableSize;                 //0 if all the lines are of the same length. Otherwise 1  
		int mode;                         //mode used to store the sequence
		int lossy;                            //lossy parameter to apply
		cds_word sizeQualLine;            //size of the qual line in case all of them have the same size
		cds_word globalMax, globalMin;
		size_t init_qual; 								//initial point in the file where the qual info start to be stored
		size_t end_qual;

		ofstream fileQual;
		string tmp_name;
		cds_word *buffer;
		cds_word buffer_use;
		string tmp_file_string;

		virtual void storeOnlyQuality(ofstream &fp);
		virtual void storeQualityASCII(ofstream &fp);
		virtual void storeQualityGlobalmm(ofstream &fp);
		virtual void storeQualityLocalmm(ofstream &fp);
		virtual string DecompressLine(cds_word *data, cds_word *PosData);

		cds_word mem_run, mem_rep;
};

#include "CQualPBlock.h"
#include "CQualRBlock.h"
#include "CQualLL.h"

#endif
