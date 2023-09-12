# consensusCruncher

Workflow to run extract UMIs from fastq and generate consensus Bams as well as run it thru mutect2 task and combinevariants task

## Overview

## Dependencies

* [hg19-bwa-index 0.7.12](http://bio-bwa.sourceforge.net/)
* [samtools 1.9](http://www.htslib.org/)
* [python 3.6](https://www.python.org/downloads/)
* [picard 2.21.2](https://broadinstitute.github.io/picard/)
* [rstats 3.6](https://www.r-project.org/)
* [consensuscruncer-5.0](https://github.com/pughlab/ConsensusCruncher)


## Usage

### Cromwell
```
java -jar cromwell.jar run consensusCruncherWorkflow.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Prefix to use for output file
`intervalFile`|String|interval file to subset variant calls
`inputIntervalsToParalellizeBy`|String|intervals for parallelization
`tumorName`|String|Name of the tumor sample
`reference`|String|reference version


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputGroups`|Array[InputGroup]?|None|Array of fastq files to concatenate if a top-up
`sortedBam`|File?|None|Bam file from bwamem
`sortedBai`|File?|None|Bai file from bwamem


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`concat.threads`|Int|4|Number of threads to request
`concat.jobMemory`|Int|16|Memory allocated for this job
`concat.timeout`|Int|72|Hours before task timeout
`concat.modules`|String|"tabix/0.2.6"|Required environment modules
`align.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`align.bwa`|String|"$BWA_ROOT/bin/bwa"|Path to bwa binary
`align.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`align.threads`|Int|4|Number of threads to request
`align.jobMemory`|Int|16|Memory allocated for this job
`align.timeout`|Int|72|Hours before task timeout
`consensus.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`consensus.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`consensus.ccDir`|String|basePrefix + ".consensuscruncher"|Placeholder
`consensus.cutoff`|Float|0.7|Cutoff to use to call a consenus of reads
`consensus.threads`|Int|8|Number of threads to request
`consensus.jobMemory`|Int|32|Memory allocated for this job
`consensus.timeout`|Int|72|Hours before task timeout
`mutectRunDCSSC.filter_timeout`|Int|12|Hours before task timeout
`mutectRunDCSSC.filter_memory`|Int|16|Memory allocated for job
`mutectRunDCSSC.filter_filterExtraArgs`|String?|None|Extra arguments
`mutectRunDCSSC.mergeStats_timeout`|Int|5|Hours before task timeout
`mutectRunDCSSC.mergeStats_memory`|Int|4|Memory allocated for job
`mutectRunDCSSC.mergeStats_modules`|String|"gatk/4.1.6.0"|Names and versions of modules to load
`mutectRunDCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunDCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunDCSSC.runMutect2_timeout`|Int|24|Hours before task timeout
`mutectRunDCSSC.runMutect2_memory`|Int|32|Memory allocated for job
`mutectRunDCSSC.runMutect2_threads`|Int|4|Number of threads to request
`mutectRunDCSSC.runMutect2_mutect2ExtraArgs`|String?|None|Extra arguments
`mutectRunDCSSC.runMutect2_mutectTag`|String|"mutect2"|Tag
`mutectRunDCSSC.splitStringToArray_modules`|String|""|Names and versions of modules to load
`mutectRunDCSSC.splitStringToArray_timeout`|Int|1|Hours before task timeout
`mutectRunDCSSC.splitStringToArray_memory`|Int|1|Memory allocated for job
`mutectRunDCSSC.splitStringToArray_lineSeparator`|String|","|line separator
`mutectRunDCSSC.normalBam`|File?|None|Input normal file (bam or sam)
`mutectRunDCSSC.normalBai`|File?|None|Index file for normal bam
`mutectRunDCSSC.pon`|File?|None|pon
`mutectRunDCSSC.ponIdx`|File?|None|pon ID
`mutectRunDCSSC.gnomad`|File?|None|gnomad
`mutectRunDCSSC.gnomadIdx`|File?|None|gnomad ID
`mutectRunSSCSSC.filter_timeout`|Int|12|Hours before task timeout
`mutectRunSSCSSC.filter_memory`|Int|16|Memory allocated for job
`mutectRunSSCSSC.filter_filterExtraArgs`|String?|None|Extra arguments
`mutectRunSSCSSC.mergeStats_timeout`|Int|5|Hours before task timeout
`mutectRunSSCSSC.mergeStats_memory`|Int|4|Memory allocated for job
`mutectRunSSCSSC.mergeStats_modules`|String|"gatk/4.1.6.0"|Names and versions of modules to load
`mutectRunSSCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunSSCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunSSCSSC.runMutect2_timeout`|Int|24|Hours before task timeout
`mutectRunSSCSSC.runMutect2_memory`|Int|32|Memory allocated for job
`mutectRunSSCSSC.runMutect2_threads`|Int|4|Number of threads to request
`mutectRunSSCSSC.runMutect2_mutect2ExtraArgs`|String?|None|Extra arguments
`mutectRunSSCSSC.runMutect2_mutectTag`|String|"mutect2"|Tag
`mutectRunSSCSSC.splitStringToArray_modules`|String|""|Names and versions of modules to load
`mutectRunSSCSSC.splitStringToArray_timeout`|Int|1|Hours before task timeout
`mutectRunSSCSSC.splitStringToArray_memory`|Int|1|Memory allocated for job
`mutectRunSSCSSC.splitStringToArray_lineSeparator`|String|","|line separator
`mutectRunSSCSSC.normalBam`|File?|None|Input normal file (bam or sam)
`mutectRunSSCSSC.normalBai`|File?|None|Index file for normal bam
`mutectRunSSCSSC.pon`|File?|None|pon
`mutectRunSSCSSC.ponIdx`|File?|None|pon ID
`mutectRunSSCSSC.gnomad`|File?|None|gnomad
`mutectRunSSCSSC.gnomadIdx`|File?|None|gnomad ID
`mutectRunAllUnique.filter_timeout`|Int|12|Hours before task timeout
`mutectRunAllUnique.filter_memory`|Int|16|Memory allocated for job
`mutectRunAllUnique.filter_filterExtraArgs`|String?|None|Extra arguments
`mutectRunAllUnique.mergeStats_timeout`|Int|5|Hours before task timeout
`mutectRunAllUnique.mergeStats_memory`|Int|4|Memory allocated for job
`mutectRunAllUnique.mergeStats_modules`|String|"gatk/4.1.6.0"|Names and versions of modules to load
`mutectRunAllUnique.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunAllUnique.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunAllUnique.runMutect2_timeout`|Int|24|Hours before task timeout
`mutectRunAllUnique.runMutect2_memory`|Int|32|Memory allocated for job
`mutectRunAllUnique.runMutect2_threads`|Int|4|Number of threads to request
`mutectRunAllUnique.runMutect2_mutect2ExtraArgs`|String?|None|Extra arguments
`mutectRunAllUnique.runMutect2_mutectTag`|String|"mutect2"|Tag
`mutectRunAllUnique.splitStringToArray_modules`|String|""|Names and versions of modules to load
`mutectRunAllUnique.splitStringToArray_timeout`|Int|1|Hours before task timeout
`mutectRunAllUnique.splitStringToArray_memory`|Int|1|Memory allocated for job
`mutectRunAllUnique.splitStringToArray_lineSeparator`|String|","|line separator
`mutectRunAllUnique.normalBam`|File?|None|Input normal file (bam or sam)
`mutectRunAllUnique.normalBai`|File?|None|Index file for normal bam
`mutectRunAllUnique.pon`|File?|None|pon
`mutectRunAllUnique.ponIdx`|File?|None|pon ID
`mutectRunAllUnique.gnomad`|File?|None|gnomad
`mutectRunAllUnique.gnomadIdx`|File?|None|gnomad ID
`hsMetricsRunDCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunDCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunDCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunDCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunDCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunDCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunDCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunSSCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunSSCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunSSCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunSSCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunAllUnique.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunAllUnique.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunAllUnique.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunAllUnique.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunAllUnique.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunAllUnique.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunAllUnique.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`combineVariants.workflows`|Array[String]|["mutect2-dcsSc", "mutect2-sscsSc"]|array of ids of producer workflows
`combineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours
`combineVariants.threads`|Int|8|number of cpu threads to be used
`annotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`annotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`annotation.timeout`|Int|20|timeout in hours
`annotation.threads`|Int|8|number of cpu threads to be used
`variantEffectPredictor.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeVcfs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`variantEffectPredictor.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`variantEffectPredictor.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`variantEffectPredictor.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeMafs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`variantEffectPredictor.vcf2maf_timeout`|Int|48|Hours before task timeout
`variantEffectPredictor.vcf2maf_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vcf2maf_bufferSize`|Int|200|The buffer size
`variantEffectPredictor.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`variantEffectPredictor.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vcf2maf_species`|String|"homo_sapiens"|Species name
`variantEffectPredictor.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.vep_timeout`|Int|16|Hours before task timeout
`variantEffectPredictor.vep_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vep_species`|String|"homo_sapiens"|Species name
`variantEffectPredictor.vep_addParam`|String?|None|Additional vep parameters
`variantEffectPredictor.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.subsetVcf_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`variantEffectPredictor.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`variantEffectPredictor.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.chromosomeArray_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`variantEffectPredictor.getSampleNames_timeout`|Int|1|Hours before task timeout
`variantEffectPredictor.getSampleNames_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.targetBedTask_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.normalName`|String?|None|Name of the normal sample


### Outputs

Output | Type | Description
---|---|---
`rawBam`|File?|raw bame file
`rawBamIndex`|File?|raw bam file index
`dcsScBam`|File|dcSc bam
`dcsScBamIndex`|File|dcSc bam index
`allUniqueBam`|File|all unique bam
`allUniqueBamIndex`|File|all unique bam index
`sscsScBam`|File|sscsSc bam
`sscsScBamIndex`|File|sscsSc bam index
`outputCCStats`|File|CC stats
`outputCCReadFamilies`|File|CC read families
`ccFolder`|File|cc folder
`dcsScVcf`|File|dcsSc vcf
`dcsScVcfIndex`|File|dcsSc vcf index
`allUniqueVcf`|File|all unique vcf
`allUniqueVcfIndex`|File|all unique vcf index
`sscsScVcf`|File|sscsSc vcf
`sscsScVcfIndex`|File|sscsSc vcf index
`vepVcf`|File|vep vcf
`vepVcfIndex`|File|vep vcf index
`mafOutput`|File?|maf output
`dcsScHsMetrics`|File|dcsSc metrics
`sscsScHsMetrics`|File|sscsSc metrics
`allUniqueHsMetrics`|File|all unique metrics


## Commands
 This section lists command(s) run by consensusCruncher workflow
 
 * Running consensusCruncher
 
 === Description here ===.
 
 ```
     set -euo pipefail
 
     zcat READ1S_ARRAY | gzip > OUTPUT_FILE_NAME_R1_001.fastq.gz
 
     zcat READ2S_ARRAY | gzip > OUTPUT_FILE_NAME_R2_001.fastq.gz
 
   ```
 ```
     set -euo pipefail
 
     CONSENSUS_CRUNCHER fastq2bam \
          --fastq1 FASTQR1 \
          --fastq2 FASTQR2 \
          --output . \
          --bwa BWA \
          --ref BWA_REF \
          --samtools SAMTOOLS \
          --skipcheck \
          --blist BLIST
 
     # Necessary for if bam files to be named according to merged library name
     # Additionally if ".sorted" isn't omitted here, file names from align include ".sorted" twice
     mv bamfiles/*.bam bamfiles/"OUTPUT_FILE_NAME.bam"
     mv bamfiles/*.bai bamfiles/"OUTPUT_FILE_NAME.bam.bai"
   ```
 ```
   set -euo pipefail
 
    CONSENSUS_CRUNCHER consensus \
          --input INPUT_BAM \
          --output . \
          --samtools SAMTOOLS \
          --cutoff CUTOFF \
          --genome GENOME \
          --bedfile CYTOBAND \
          --bdelim '|'
 
    tar cf - NAME_PREFiX | gzip --no-name > CC_DIR.tar.gz
   ```
 ```
   python3<<CODE
   import subprocess
   import sys
   inputStrings = []
   v = sep=' ' INPUT_VCFS
   vcfFiles = v.split()
   w = sep=' ' WORKFLOWS
   workflowIds = w.split()
   priority = PRIORITY
   
   if len(vcfFiles) != len(workflowIds):
       print("The arrays with input files and their respective workflow names are not of equal size!")
   else:
       for f in range(0, len(vcfFiles)):
           inputStrings.append("--variant:" + workflowIds[f] + " " + vcfFiles[f])
 
   javaMemory = JOB_MEMORY - 6 
   gatkCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar $GATK_ROOT/GenomeAnalysisTK.jar "
   gatkCommand += "-T CombineVariants "
   gatkCommand += " ".join(inputStrings)
   gatkCommand += " -R REFERENCE_FASTA "
   gatkCommand += "-o OUTPUT_FILE_NAME_combined.vcf.gz "
   gatkCommand += "-genotypeMergeOptions PRIORITIZE "
   gatkCommand += "-priority " + priority
   gatkCommand += " 2>&1"
 
   result_output = subprocess.run(gatkCommand, shell=True)
   sys.exit(result_output.returncode)
   CODE
 ```
 ```
   bcftools annotate -a UNIQUE_VCF \
  -c FMT/AD,FMT/DP MERGED_VCF -Oz \
  -o "OUTPUT_FILE_NAME.merged.vcf.gz"
 
  tabix -p vcf "OUTPUT_FILE_NAME.merged.vcf.gz"
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
