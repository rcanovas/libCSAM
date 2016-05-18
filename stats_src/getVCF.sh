#!/bin/bash

ref=$1
bam=$2
vcf=$bam".vcf"
samtools mpileup -ugf $ref $bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | /usr/share/samtools/vcfutils.pl varFilter -D100 > $vcf


