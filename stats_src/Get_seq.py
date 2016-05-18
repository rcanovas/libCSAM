''' libCSAM
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
'''

import sys
import re

if len(sys.argv)!=2:
	print 'Use: python ./Get_seq.py SAMfile'
	sys.exit('Error')

arg1 = sys.argv[1]
f=open(arg1,'r')
namefile = arg1 + '.seq'
nameref = arg1 + '.ref'
namepos = arg1 + '.pos'
seq = open(namefile, 'w')
ref = open(nameref, 'w')
pos = open(namepos, 'w')

print "parsing"
while(1):
	"""get each line"""
	line = f.readline()
	if not line:
		break
	if line[0]!='@':	
		a=line.split('\t')
		ref.write(a[2] + '\n')
		pos.write(a[3] + '\n')
		seq.write(a[9] + '\n')
		

