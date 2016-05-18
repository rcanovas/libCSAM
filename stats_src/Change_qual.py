'''
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

if len(sys.argv)!=3:
	print 'Use: python ./Change_qual file.sam new_qual.qual'
	sys.exit('Error')

arg1 = sys.argv[1]
f=open(arg1,'r')
arg2 = sys.argv[2]
new_quals=open(arg2,'r')

newSAM = open('newSAM.sam', 'w')


print "parsing"
while(1):
	"""get each line"""
	line = f.readline()
	if not line:
		break
	if line[0]!='@':
		qual = new_quals.readline()	
		a=line.split('\t')
		newLine = a[0]+'\t'+a[1]+'\t'+a[2]+'\t'+a[3]+'\t'+a[4]+'\t'+a[5]+'\t'
		newLine += a[6]+'\t'+a[7]+'\t'+a[8] +'\t'+a[9]+'\t'
		qual=qual.rstrip('\n')
		newLine += qual
		if len(a)>11:
			b=12
			l = a[11]
			while 	b< len(a):
	 			l += '\t' + a[b]
				b = b+1	
			newLine += '\t'+l
		else:
			newLine += '\n'
		newSAM.write(newLine)
	else:	
		newSAM.write(line)

