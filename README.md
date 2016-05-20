## libCSAM

libCSAM contains several C++ codes for compress,decompress, and 
access each of the fields of any SAM format file. Part of this library 
was taking from Francisco Claude libcds project 
(https://github.com/fclaude/libcds/). Also the boost C++ library must 
be installed in your computer (http://www.boost.org/).

## Methods
-[CompressQual] (https://github.com/tests)

-[CompressSAM] (https://github.com/tests)

-[CompressSeq] (https://github.com/tests)

-[CountReadsSample] (https://github.com/tests)

-[DecompressQual] (https://github.com/tests):
    
       Use: ./DecompressQual <arch>
       arch: .cqual File
       output:.qual File containing the quality scores*

-[DecompressSAM] (https://github.com/tests)

-[DescompressSeq] (https://github.com/tests)

-[GetIntervalSAM] (https://github.com/tests)

-[GetIntervalSAMSample] (https://github.com/tests)

-[GetIntervalSeq] (https://github.com/tests)

-[GetIntervalSeqSample] (https://github.com/tests)

-[GetIntervalSSN] (https://github.com/tests)


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
        output:  file.sam.seq

-[getVCF] (https://github.com/stats_src): Simple example of how to generate the vcf file of a BAM file using mpileup and bcftools.

        Use: ./getVCF reference_file file.bam
        output:  file.bam.vcg

-[get_vcf_stats] (https://github.com/stats_src): Compares two vcf files computing true positive, false positive, false negative, precision, recall, and MSE.

        Use: python ./get_vcf_stats.py original.vcf second.vcf
        output: Returns stats in screen




Note: These codes assume that the computer have enough RAM memory to read and store the complete input.
