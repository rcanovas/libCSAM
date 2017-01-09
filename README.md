## libCSAM

libCSAM contains several C++ codes for compress,decompress, and 
access each of the fields of any SAM format file. Part of this library 
was taking from Francisco Claude libcds project 
(https://github.com/fclaude/libcds/). Also the boost C++ library must 
be installed in your computer (http://www.boost.org/).

## Methods

-[CompressSAM] (https://github.com/tests):

      Use: ./CompressSAM <arch> <opt>
      opt: 
          -q qm: How the Quality values are stored. q = 0 - lossles, 1 - pblock, 2 - rblock. Default: q = 1
          -m mode: Mode use to store the Representative Array. mode = 0 ASCII, 1 Binary Global, 2 Binary Local. Default mode=1
          -l lossy: lossy parameter use to compress the quality score depending on the mode use. Default: 0
          -s sample: Sample rate used for Fields and Quality structure. Default: s = 1000
          -p position: Sample position rate used for Seq, Rname, and Pos. Default: p = 1000
          
          	output: .cqual File

		Example: ./CompressSAM ./data/file.sam -q 2 -l 60 -s 500 -p 100
		output:  file.sam.csam
        

-[CompressQual] (https://github.com/tests):

	Use: ./CompressQual <arch> <opt>
	arch: 	SAM format file. Note that only the quality field will be use.
	opt: 
		-q mode: 	How the Quality values are stored. 
			  	mode=0 	gzip 
			  	mode=1 	P-Block  
			  	mode=2 	R-Block
 			  	mode=3 	Bins base on the LogBinning Wan et al.
 					2011 paper. Note that UniBinning is 
				       	also implemented but is not included 
				       	in this program. 
				mode=4 	Only one value is stored to represent
				        all the qualities  
			  	default:  0
		-m mode: 	If 'q' is 1 or 2, give the mode use to store
 						the Representative Array. 
				 		mode=0 	ASCII 
				 		mode=1 	Binary Global 
				 		mode=2 	Binary Local
				 		default=1
		-l lossy: 	Lossy parameter use to compress the quality
						score depending on the mode use. 
						P-Block: maximum distance between the values
							and their representative (1,2,3,4)
						R-BLock: Max/Min maximum diference allowed,
							 recieve the extra overhead pbb 
							 (5, 20, 40, 100) 
						Bins:    number of bins to be use (94, 20,
							 10, 5)
		-r :     	Reorder the quality scores by reference and
						position. Also the permutation is stored
		-s sample:  Size of the sample rate that will be use.
						default: no sample. 

		output: .cqual File

		Example: ./CompressQual ./data/file.sam -q 2 -l 60 -s 500 -r
		output:  file.sam.cqual
		Compress the quality scores of file.sam using R-Block with r=1.60 
		storing a sample every 500 lines and reordering the file 
		by reference and position in the reference.


-[CompressSeq] (https://github.com/tests):

      Use: ./CompressSeq <arch> <opt>
      <arch>:  must be a .sam or .rps (rname pos seq) file
      opt: 
          -s sample:  size of the sample rate that will be use. Default: no sample;
          
      output: .cseq file

-[DecompressSAM] (https://github.com/tests):

       Use: ./DecompressSAM <arch>
       arch: .csam File
       output:.sam File containing the SAM information


-[DecompressQual] (https://github.com/tests):
    
       Use: ./DecompressQual <arch>
       arch: .cqual File
       output:.qual File containing the quality scores

-[DescompressSeq] (https://github.com/tests):

       Use: ./DecompressSeq <arch>
       arch: .cseq File
       output:.seq File containing only the sequence field

-[CountReadsSample] (https://github.com/tests): Counts the number of read within each of the intervals in the sample_interval size. Also gives some stats about the interval found.

        Use: ./CountReadsSample <arch>.csam sample_interval_file
        output: On screen


-[GetIntervalSAM] (https://github.com/tests): Extracs from a csam file all the alignment lines wihtin the interval (ref,x,y)
    
        Use: ./GetIntervalSAM <arch> ref_name pos_x pos_y
        arch: .csam File
        ref: reference name
        pos_x, pos_y:  interval positions
        output: file_name + "_inter.sam" File

-[GetIntervalSeq] (https://github.com/tests): Same as before but only extracting only the SEQ fields.

-[GetIntervalSAMSample] (https://github.com/tests): Same as GetIntervalSAM but receive a file containing many intervals to query

           Use: ./GetIntervalSAMSample <arch>.csam sample_interval_file
           Use:	./GetIntervalSAMSample <arch>.csam sample_interval_file BuffSizeInBytes


-[GetIntervalSeqSample] (https://github.com/tests): Same as before but only extracting the SEQ field

-[GetIntervalSSN] (https://github.com/tests): Same as GetIntervalSAMSample but extracting only a selection of the Fields and replacing the rest with empty values. For the moment it extrac a minimal set (QNAME FLAG RNAME POS MAPQ SEQ). Modify line 113 
of the file to change this option (TODO: do it by command line)

           Use: ./GetIntervalSSN <arch>.csam sample_interval_file
           output: file_name + "_inter.sam" File



## Stats Methods

Also this library contains in the  stats_src the following programs:


-[Change_qual] (https://github.com/stats_src):	Changes the quality field of a SAM file with the quality file given	

        Use: python ./Change_qual.py file.sam new_qual.qual
        output: newSAM.sam

-[Change_qual_letter] (https://github.com/stats_src): Changes the quality field of a SAM file to only one quality score value

        Use: python ./Change_qual_letter.py file.sam letter name_output.sam
        output: name_output.sam

-[ComputeEntroHist] (https://github.com/stats_src): Computes the Entropy of order 0 of a file return the histogram of each symbol

      (compile first: g++ -o ComputeMetrics ComputeMetrics.cpp)
      Use: ./ComputeEntroHist <arch> <out_arch>
      output: In screen prints the entropy of the file, and in <out_arch> returns the histogram of the symbols

-[ComputeMetrics] (https://github.com/stats_src): Compares two quality files and compute some distance metrics

        (compile first: g++ -o ComputeMetrics ComputeMetrics.cpp)
        Use: ./ComputeMetrics qualityFile.qual referenceFile.qual
        output : In screen prints the Manhattan, Max:Min, MSE, Chebyshev, Soergel and Lorentzian metrics.

-[Get_qual] (https://github.com/stats_src):	Extracs the quality field of a SAM file

        Use: python ./Get_qual.py file.sam
        output:  file.sam.qual


-[Get_seq] (https://github.com/stats_src): Extracs the reference, positon and sequence field of a SAM file

        Use: python ./Get_seq.py file.sam
        output:  file.sam.rps

-[getVCF] (https://github.com/stats_src): Simple example of how to generate the vcf file of a BAM file using mpileup and bcftools.

        Use: ./getVCF reference_file file.bam
        output:  file.bam.vcg

-[get_vcf_stats] (https://github.com/stats_src): Compares two vcf files computing true positive, false positive, false negative, precision, recall, and MSE.

        Use: python ./get_vcf_stats.py original.vcf second.vcf
        output: Returns stats in screen


Note: These codes assume that the computer have enough RAM memory to read and store the complete input.
