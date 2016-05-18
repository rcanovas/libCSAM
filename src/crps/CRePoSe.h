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
#ifndef _CREPOSE_H
#define _CREPOSE_H

#include "./../basics/buffer.h"
#include "./../coders/Coder.h"
#include "./../coders/Huffman.h"

#include "./../cfields/CFields.h"
#include "./../cqual/CQualBlock.h"


using namespace cds;
using namespace cds::basic;
using namespace std;

const cds_word ROW = 1;
const cds_word WINSIZE = 700; //We assume that the maximum length of a read is 350
const cds_word MAXREADLENGTH = 350;  //We assumed that in the reads are always of length lower than 200, but this is just in case.


class CRePoSe{

	public:
		virtual ~CRePoSe() {};
		/*Read line Seq with reference Ref aligned at position Pos and compute 
		 * data neccesary to generate the pressumed sequence*/
		void ProcessLine(string Ref, cds_word Pos, string Seq, vector<cds_word>* IndexClusterToLines = NULL);
		/*Read Sam file and generate a cseq file containing the compress SEQ, POS, and 
		 * RNAME data. SamOrRps indicates if the input file is a .sam or .rps (rname pos seq) file*/
		virtual void CreateFile(ifstream &SamFile, ofstream &fp, int SamOrRps, vector<cds_word>* IndexClusterToLines = NULL) = 0;
		/*Decompress all the reads*/
		virtual void DecompressSequence(ifstream &input, ofstream &output) = 0;
		/*Decompress the whole SAM file*/
		virtual void DecompressAll(ifstream **input, ofstream &output, CQualBlock *cqual, CFields *cfields) = 0;
		/*Get interval of only reads, their positions, and reference. Return the number of reads in the interval (input can be .cseq or .csam)*/
		virtual cds_word GetInterval(ifstream &input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, size_t MaxBytes = MAXRAMBYTESDECO) = 0;
		/*Get interval of complete alignment reads (input need to be a .csam)*/
		virtual cds_word GetIntervalAll(ifstream **input, ofstream &output, string rname, cds_word pos_x, cds_word pos_y, CQualBlock *cqual, CFields *cfields, vector<cds_word>* ClusToLine, int selection, size_t MaxBytes = MAXRAMBYTESDECO) = 0;
		/*Count the number of reads presented in an interval also considering if they are reverse (=1) or not (=0)*/
		virtual cds_word CountReads(ifstream **input, string rname, cds_word pos_x, cds_word pos_y, int reverse, CFields *cfields, vector<cds_word>* ClusToLine, size_t MaxBytes = MAXRAMBYTESDECO) = 0;

		cds_word GetNumberOfLines();

		static CRePoSe *Load(ifstream &fp);



	protected:
	  cds_word CompressLine(string Ref, cds_word Pos, string Seq, string PSeq, cds_word **buff, cds_word *init_buff);
		string DecompressLine(cds_word *data, cds_word *PosData, string preseq, cds_word *actualPos, bool *first, cds_word *last_lengh, cds_word *rel_pos);
		void initData(cds_word sample, bool create);
		void getValuesPreSeq(int since, cds_word n);
		void restartWinSeq(string seq);
		void CompressPressumeSeq(ofstream &fp, cds_word **buff, cds_word *init_buff);
		string GetPressumeSeq(ifstream &input, cds_word ini, cds_word end);
		cds_word getFirstCluster(string Ref, cds_word *next);

		vector <string> RNames;				// array with reference names 
		vector <cds_word> RNamesCluster;  // first appear line of each reference

		vector <cds_word> ClusterIndex;
		vector <cds_word> ClusterOverlap;
		cds_word IndexRate;           //Every IndexRate lines we store a pointer to fast access

		cds_word numberOfLinesNR;     //number of lines where the reference is "*"
		cds_word numberOfLines;       //total number of lines stored
		int variableSize;             //0 if all the lines are of the same length. Otherwise 1  
		vector <cds_word> sizeLine;            //size of the read line in case all of them have the same size
		size_t init_seq;						  //initial point in the file where sequence info start to be stored
		size_t init_preseq;
		size_t end_seq;


		Huffman *huffA, *huffC, *huffG, *huffT; //huffman codes for the not sets
		cds_word preSeqLength;
		/*auxiliar variable*/
		string last_ref;
		cds_word last_pos, farther_pos;
		cds_word overlapcluster;
		cds_word bitsLength;         // bits used to represent the length
		int **WinSeq;
		int posWin;
		cds_word ClusterMaxPos;       // use to record the index positions
		cds_word ClusterReads;     // use to compute number of reads per cluster
		cds_word MaxClusterReads;  //maximum number of reads find inside of a Cluster
		cds_word MaxDepthCluster;  //maximum number od bases in a column

		/*stats for huffman*/
		cds_uint occ_notA[250], occ_notC[250], occ_notG[250], occ_notT[250]; //number of occurences of each base

		/*buffer used for decompress*/
		cds_word *buffRPS;
		cds_word *buffQual;
		char *buffField;
		cds_word init_buffQual, init_buffRPS, init_buffField;

		cds_word *buffer;
		cds_word buffer_use;
		string tmp_file_string;
		string tmp_name;                 
		ofstream fileSeq;
};

#include "CRPSPreMF.h"

#endif
