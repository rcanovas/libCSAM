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

def convertStr(s):
	"""Convert string to either int or float."""
 	try:
 		ret = int(s)
 	except ValueError:
 		"""Try float."""
 		ret = float(s)
 	return ret

if len(sys.argv)!=3:
	print 'Use: python ./get_vcf_stats.py original.vcf lossy.vcf'
	sys.exit('Error')

arg1 = sys.argv[1]
arg2 = sys.argv[2]

original=open(arg1,'r')
lossy=open(arg2,'r')

or_pos=[]
or_ref=[]
or_new=[]
or_qv=[]
or_numSNP = 0

while(1):
	"""get each line from the original file"""
	line = original.readline()
	if not line:
		break
	if line[0]!='#':	
		a=line.split()
		or_pos.append(convertStr(a[1]))
		or_ref.append(a[3])
		or_new.append(a[4])
		or_qv.append(convertStr(a[5]))		
		or_numSNP = or_numSNP + 1 
s='number Original SNP: '+repr(or_numSNP)
print s

lo_pos=[]
lo_ref=[]
lo_new=[]
lo_qv=[]
lo_numSNP = 0

while(1):
	"""get each line from the lossy file"""
	line = lossy.readline()
	if not line:
		break
	if line[0]!='#':	
		a=line.split()
		lo_pos.append(convertStr(a[1]))
		lo_ref.append(a[3])
		lo_new.append(a[4])
		lo_qv.append(convertStr(a[5]))		
		lo_numSNP = lo_numSNP + 1 
s='New number SNP: '+repr(lo_numSNP)
print s

""""Compute true positive, false positive, and false negative"""
Half_TP=0
Equal_TP=0
mse = 0
FP=0
FN=0
j=0
while(j<lo_numSNP):
	if lo_pos[j] in or_pos:
		i=or_pos.index(lo_pos[j])
		mse += pow(or_qv[i] - lo_qv[j], 2)
		if lo_ref[j]==or_ref[i] and lo_new[j]==or_new[i]:
			Equal_TP=Equal_TP+1
		else:
			Half_TP=Half_TP+1			
	else:
		FP=FP+1
	
	j=j+1
s='FP: '+repr(FP)+'  FN: '+repr(or_numSNP-Half_TP-Equal_TP)+'  TP: '+repr(Half_TP+Equal_TP)  
print s
s='Half TP: '+repr(Half_TP)+'  Equal TP: '+repr(Equal_TP)
print s
s='Precision (TP*100.0/(TP+FP)): '+repr((Half_TP+Equal_TP)*100.0/((Half_TP+Equal_TP)+FP))
print s
s='Recall (TP*100.0/(TP+FN)): '+repr((Half_TP+Equal_TP)*100.0/((Half_TP+Equal_TP)+(or_numSNP-Half_TP-Equal_TP)))
print s
mse = mse/(Half_TP+Equal_TP)
s='MSE: '+repr(mse)
print s

