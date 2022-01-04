#!/bin/bash
set -o nounset
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sed 's/.*\.//' | sort | uniq -c

find -name *.bam -exec samtools flagstat {} \; | sort

find -name *.bam -exec /bin/bash -c "samtools view {} | md5sum" \; | sort

for v in *.vcf;do cat $v | grep -v ^# | md5sum; done | sort -V

