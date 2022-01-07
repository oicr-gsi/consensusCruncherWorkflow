#!/bin/bash
set -o nounset
set -o pipefail

cd $1

module load samtools/1.9 vcftools/0.1.16 2>/dev/null

ls | sed 's/.*\.//' | sort | uniq -c

for i in *.stats.txt; do md5sum ${i}; done | sort

for i in *.bam; do samtools flagstat ${i} | grep "in total" | md5sum; done | sort

for i in *.vcf.gz; do vcf-query -l ${i}; done | sort

for i in *.vcf.gz; do zcat ${i} | grep -v ^# | cut -f 1 | uniq| md5sum; done | sort
