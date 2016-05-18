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

#ifndef SRC_BASIC_EXTRAS_H_
#define SRC_BASIC_EXTRAS_H_

#include <cmath>
#include <algorithm>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "./libcds.h"
#include "./io.h"



using namespace cds;
using namespace cds::basic;
using namespace std;

using namespace boost::archive;
using namespace boost::iostreams;

inline string GZcompressString(string s){
	stringstream zippedStream;
	{ //do not erase 
		filtering_ostream filteringStream;
		filteringStream.push(gzip_compressor());
		filteringStream.push(zippedStream);
		filteringStream << s;
	}
	return  zippedStream.str();
}

inline string GZdecompressString(string s){
	stringstream ss, unzippedStream;               
	ss << s;
	filtering_istream filteringStream;
	filteringStream.push(gzip_decompressor());	
	filteringStream.push(ss);
	copy(filteringStream, unzippedStream);
	return unzippedStream.str();
}

inline string NumtoString(cds_word x){
	stringstream ss;
	ss << x;
	return ss.str();
}


#endif




