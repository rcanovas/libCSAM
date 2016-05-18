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
#ifndef _CSAM_H
#define _CSAM_H

#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "./../crps/CRePoSe.h"
#include "./../cqual/CQualBlock.h"
#include "./../cfields/CFields.h"

const int ALL =  4095;

using namespace cds;
using namespace cds::basic;
using namespace std;


class CSAM{

	public:
		/*CSAM constructor
		 * @param qualmode: Qual compress version that will be use 
		 * @param lossyp: lossy parameter use to compress Qual
		 * @param sample: Rate of alignment lines that will use to set sample points*/
		CSAM(int empty, cds_word qualmode = 1, cds_word encodemode = 1, int lossyp = 0, cds_word sampleRPS = 1000, cds_word sampleQF = 1000);
		
		~CSAM();

		void CompressSAM(string filename);

		/*regenerate the complete SAM file*/
		void DecompressSAM(ofstream &output);

		/*Get interval of only reads, their positions, and reference. Return the number of reads in the interval*/
		virtual cds_word GetInterval(ofstream &output, string rname, cds_word pos_x, cds_word pos_y, bool header = false, int selection = ALL, size_t MaxBytes = MAXRAMBYTESDECO);

		/*Count the number of reads presented in an interval also considering if they are reverse (=1) or not (=0)*/
		virtual cds_word CountReads(string rname, cds_word pos_x, cds_word pos_y, int reverse);
		
		static CSAM *Load(string filename);
																		
		cds_word GetNumberOfLines();

		cds_word max_DN;

	protected:
		ifstream *fileCSAM; //3 file descriptor used to acces the SEQ, QUAL, and other fields structures

		string cheader;
		cds_word sampleRate, posSample;
		CRePoSe *crps;  //Contains the Fields SEQ, RNAME, and POS.
		CQualBlock *cqual;
		CFields *cfields;
		vector <cds_word> IndexClusterToLines; //Index containing the relationship between the CRePoSe cluster and number of lines visited
		cds_word NumberofLines;

	private:

		void CheckCIGAR(string cigar, cds_word len);

};

#endif
