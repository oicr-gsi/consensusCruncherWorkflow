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
`mutectRunDCSSC.gatk`|String|gatk version to be used
`mutectRunDCSSC.outputFileNamePrefix`|String|prefix of output file
`mutectRunSSCSSC.gatk`|String|gatk version to be used
`mutectRunSSCSSC.outputFileNamePrefix`|String|prefix of output file
`mutectRunAllUnique.gatk`|String|gatk version to be used
`mutectRunAllUnique.outputFileNamePrefix`|String|prefix of output file
`hsMetricsRunSSCSSC.reference`|String|the reference genome for input sample
`hsMetricsRunAllUnique.reference`|String|the reference genome for input sample
`combineVariants.workflows`|Array[String]|array of ids of producer workflows


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
`mutectRunDCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunDCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunDCSSC.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mutectRunDCSSC.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`mutectRunDCSSC.runMutect2_threads`|Int|4|Requested CPU threads
`mutectRunDCSSC.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`mutectRunDCSSC.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`mutectRunDCSSC.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`mutectRunDCSSC.getChrCoefficient_memory`|Int|1|Memory allocated for this job
`mutectRunDCSSC.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`mutectRunDCSSC.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`mutectRunDCSSC.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`mutectRunDCSSC.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunDCSSC.normalBai`|File?|None|Index for noramlBam
`mutectRunDCSSC.pon`|File?|None|panel of normal
`mutectRunDCSSC.ponIdx`|File?|None|index of pon
`mutectRunSSCSSC.filter_timeout`|Int|12|Hours before task timeout
`mutectRunSSCSSC.filter_memory`|Int|16|Memory allocated for job
`mutectRunSSCSSC.filter_filterExtraArgs`|String?|None|Extra arguments
`mutectRunSSCSSC.mergeStats_timeout`|Int|5|Hours before task timeout
`mutectRunSSCSSC.mergeStats_memory`|Int|4|Memory allocated for job
`mutectRunSSCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunSSCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunSSCSSC.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mutectRunSSCSSC.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`mutectRunSSCSSC.runMutect2_threads`|Int|4|Requested CPU threads
`mutectRunSSCSSC.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`mutectRunSSCSSC.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`mutectRunSSCSSC.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`mutectRunSSCSSC.getChrCoefficient_memory`|Int|1|Memory allocated for this job
`mutectRunSSCSSC.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`mutectRunSSCSSC.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`mutectRunSSCSSC.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`mutectRunSSCSSC.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunSSCSSC.normalBai`|File?|None|Index for noramlBam
`mutectRunSSCSSC.pon`|File?|None|panel of normal
`mutectRunSSCSSC.ponIdx`|File?|None|index of pon
`mutectRunAllUnique.filter_timeout`|Int|12|Hours before task timeout
`mutectRunAllUnique.filter_memory`|Int|16|Memory allocated for job
`mutectRunAllUnique.filter_filterExtraArgs`|String?|None|Extra arguments
`mutectRunAllUnique.mergeStats_timeout`|Int|5|Hours before task timeout
`mutectRunAllUnique.mergeStats_memory`|Int|4|Memory allocated for job
`mutectRunAllUnique.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunAllUnique.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunAllUnique.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mutectRunAllUnique.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`mutectRunAllUnique.runMutect2_threads`|Int|4|Requested CPU threads
`mutectRunAllUnique.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`mutectRunAllUnique.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`mutectRunAllUnique.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`mutectRunAllUnique.getChrCoefficient_memory`|Int|1|Memory allocated for this job
`mutectRunAllUnique.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`mutectRunAllUnique.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`mutectRunAllUnique.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`mutectRunAllUnique.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunAllUnique.normalBai`|File?|None|Index for noramlBam
`mutectRunAllUnique.pon`|File?|None|panel of normal
`mutectRunAllUnique.ponIdx`|File?|None|index of pon
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
`combineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours
`combineVariants.threads`|Int|8|number of cpu threads to be used
`annotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`annotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`annotation.timeout`|Int|20|timeout in hours
`annotation.threads`|Int|8|number of cpu threads to be used
`variantEffectPredictor.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`variantEffectPredictor.mergeVcfs_jobMemory`|Int|8|Memory allocated to job (in GB).
`variantEffectPredictor.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`variantEffectPredictor.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`variantEffectPredictor.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeMafs_jobMemory`|Int|8|Memory allocated to job (in GB).
`variantEffectPredictor.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`variantEffectPredictor.vcf2maf_timeout`|Int|18|Hours before task timeout
`variantEffectPredictor.vcf2maf_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vcf2maf_jobMemory`|Int|12|Memory allocated for this job (GB)
`variantEffectPredictor.vcf2maf_bufferSize`|Int|200|The buffer size
`variantEffectPredictor.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`variantEffectPredictor.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.tumorOnlyAlign_jobMemory`|Int|12|Memory allocated for this job (GB)
`variantEffectPredictor.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.vep_timeout`|Int|16|Hours before task timeout
`variantEffectPredictor.vep_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vep_jobMemory`|Int|12|Memory allocated for this job (GB)
`variantEffectPredictor.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vep_addParam`|String?|None|Additional vep parameters
`variantEffectPredictor.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.subsetVcf_timeout`|Int|2|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.subsetVcf_jobMemory`|Int|12|Memory allocated to job (in GB).
`variantEffectPredictor.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`variantEffectPredictor.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.getChrCoefficient_timeout`|Int|1|Hours before task timeout
`variantEffectPredictor.getChrCoefficient_memory`|Int|1|Memory allocated for this job
`variantEffectPredictor.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`variantEffectPredictor.getSampleNames_timeout`|Int|1|Hours before task timeout
`variantEffectPredictor.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.targetBedTask_jobMemory`|Int|12|Memory allocated for this job (GB)
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
 
### preprocess Reads
 
```
     set -euo pipefail
 
     zcat READ1S_ARRAY | gzip > OUTPUT_FILE_NAME_R1_001.fastq.gz
 
     zcat READ2S_ARRAY | gzip > OUTPUT_FILE_NAME_R2_001.fastq.gz
 
```

### convert fastq to bam

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

### Call consensus

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

### Run gatk for combining variants

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

### annotate with bcftools

```
   bcftools annotate -a UNIQUE_VCF \
  -c FMT/AD,FMT/DP MERGED_VCF -Oz \
  -o "OUTPUT_FILE_NAME.merged.vcf.gz"
 
  tabix -p vcf "OUTPUT_FILE_NAME.merged.vcf.gz"
```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
