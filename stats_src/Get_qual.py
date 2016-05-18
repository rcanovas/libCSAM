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
	print 'Use: python ./Get_qual.py SAMfile'
	sys.exit('Error')

arg1 = sys.argv[1]
f=open(arg1,'r')
namefile = arg1 + '.qual'
qual = open(namefile, 'w')
b=0

print "parsing"
while(1):
	"""get each line"""
	line = f.readline()
	if not line:
		break
	if line[0]!='@':	
		a=line.split('\t')
		qual.write(a[10])
		if len(a)>11:
			qual.write('\n')
			b=12
			l = a[11]
			while 	b< len(a):
	 			l += '\t' + a[b]
				b = b+1	

