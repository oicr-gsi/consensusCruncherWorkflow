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
`intervalFile`|File|Missing?
`inputRefDict`|String|Missing?
`inputRefFai`|String|Missing?
`inputRefFasta`|String|Missing?
`inputMutectModules`|String|Missing?
`inputIntervalsToParalellizeBy`|String|Missing?
`inputHSMetricsModules`|String|Missing?
`combineVariants.workflows`|Array[String]|array of ids of producer workflows
`combineVariants.modules`|String|modules for running preprocessing
`variantEffectPredictor.vcf2maf_vcfFilter`|String|Filter for the vep module that is used in vcf2maf
`variantEffectPredictor.vcf2maf_vepCacheDir`|String|Directory of vep cache files
`variantEffectPredictor.vcf2maf_vepPath`|String|Path to vep script
`variantEffectPredictor.vcf2maf_ncbiBuild`|String|The assembly version
`variantEffectPredictor.vcf2maf_modules`|String|Required environment modules
`variantEffectPredictor.vep_modules`|String|Required environment modules
`variantEffectPredictor.vep_vepCacheDir`|String|Directory of cache files
`variantEffectPredictor.vep_ncbiBuild`|String|The assembly version


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR1`|File?|None|First Fastq Files
`fastqR2`|File?|None|Second Fastq File
`sortedBam`|File?|None|Bam file from bwamem
`sortedBai`|File?|None|Bai file from bwamem


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`align.modules`|String|"consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"|Names and versions of modules to load
`align.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`align.bwa`|String|"$BWA_ROOT/bin/bwa"|Path to bwa binary
`align.bwaref`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa"|Path to bwa index
`align.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`align.blist`|String|"$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"|Path to blacklist for barcodes
`align.threads`|Int|4|Number of threads to request
`align.jobMemory`|Int|16|Memory allocated for this job
`align.timeout`|Int|72|Hours before task timeout
`consensus.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`consensus.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`consensus.cytoband`|String|"$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt"|Path to cytoband for genome
`consensus.genome`|String|"hg19"|Which genome version to use
`consensus.ccDir`|String|basePrefix + ".consensuscruncher"|Missing?
`consensus.cutoff`|Float|0.7|Cutoff to use to call a consenus of reads
`consensus.threads`|Int|8|Number of threads to request
`consensus.jobMemory`|Int|32|Memory allocated for this job
`consensus.timeout`|Int|72|Hours before task timeout
`consensus.modules`|String|"consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"|Names and versions of modules to load
`mutectRunDCSSC.filter_timeout`|Int|12|Missing?
`mutectRunDCSSC.filter_memory`|Int|16|Missing?
`mutectRunDCSSC.filter_filterExtraArgs`|String?|None|Missing?
`mutectRunDCSSC.mergeStats_timeout`|Int|5|Missing?
`mutectRunDCSSC.mergeStats_memory`|Int|4|Missing?
`mutectRunDCSSC.mergeStats_modules`|String|"gatk/4.1.6.0"|Missing?
`mutectRunDCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunDCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunDCSSC.runMutect2_timeout`|Int|24|Missing?
`mutectRunDCSSC.runMutect2_memory`|Int|32|Missing?
`mutectRunDCSSC.runMutect2_threads`|Int|4|Missing?
`mutectRunDCSSC.runMutect2_mutect2ExtraArgs`|String?|None|Missing?
`mutectRunDCSSC.runMutect2_mutectTag`|String|"mutect2"|Missing?
`mutectRunDCSSC.splitStringToArray_modules`|String|""|Missing?
`mutectRunDCSSC.splitStringToArray_timeout`|Int|1|Missing?
`mutectRunDCSSC.splitStringToArray_memory`|Int|1|Missing?
`mutectRunDCSSC.splitStringToArray_lineSeparator`|String|","|Missing?
`mutectRunDCSSC.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunDCSSC.normalBai`|File?|None|Missing?
`mutectRunDCSSC.pon`|File?|None|Missing?
`mutectRunDCSSC.ponIdx`|File?|None|Missing?
`mutectRunDCSSC.gnomad`|File?|None|Missing?
`mutectRunDCSSC.gnomadIdx`|File?|None|Missing?
`mutectRunSSCSSC.filter_timeout`|Int|12|Missing?
`mutectRunSSCSSC.filter_memory`|Int|16|Missing?
`mutectRunSSCSSC.filter_filterExtraArgs`|String?|None|Missing?
`mutectRunSSCSSC.mergeStats_timeout`|Int|5|Missing?
`mutectRunSSCSSC.mergeStats_memory`|Int|4|Missing?
`mutectRunSSCSSC.mergeStats_modules`|String|"gatk/4.1.6.0"|Missing?
`mutectRunSSCSSC.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunSSCSSC.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunSSCSSC.runMutect2_timeout`|Int|24|Missing?
`mutectRunSSCSSC.runMutect2_memory`|Int|32|Missing?
`mutectRunSSCSSC.runMutect2_threads`|Int|4|Missing?
`mutectRunSSCSSC.runMutect2_mutect2ExtraArgs`|String?|None|Missing?
`mutectRunSSCSSC.runMutect2_mutectTag`|String|"mutect2"|Missing?
`mutectRunSSCSSC.splitStringToArray_modules`|String|""|Missing?
`mutectRunSSCSSC.splitStringToArray_timeout`|Int|1|Missing?
`mutectRunSSCSSC.splitStringToArray_memory`|Int|1|Missing?
`mutectRunSSCSSC.splitStringToArray_lineSeparator`|String|","|Missing?
`mutectRunSSCSSC.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunSSCSSC.normalBai`|File?|None|Missing?
`mutectRunSSCSSC.pon`|File?|None|Missing?
`mutectRunSSCSSC.ponIdx`|File?|None|Missing?
`mutectRunSSCSSC.gnomad`|File?|None|Missing?
`mutectRunSSCSSC.gnomadIdx`|File?|None|Missing?
`mutectRunAllUnique.filter_timeout`|Int|12|Missing?
`mutectRunAllUnique.filter_memory`|Int|16|Missing?
`mutectRunAllUnique.filter_filterExtraArgs`|String?|None|Missing?
`mutectRunAllUnique.mergeStats_timeout`|Int|5|Missing?
`mutectRunAllUnique.mergeStats_memory`|Int|4|Missing?
`mutectRunAllUnique.mergeStats_modules`|String|"gatk/4.1.6.0"|Missing?
`mutectRunAllUnique.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutectRunAllUnique.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutectRunAllUnique.runMutect2_timeout`|Int|24|Missing?
`mutectRunAllUnique.runMutect2_memory`|Int|32|Missing?
`mutectRunAllUnique.runMutect2_threads`|Int|4|Missing?
`mutectRunAllUnique.runMutect2_mutect2ExtraArgs`|String?|None|Missing?
`mutectRunAllUnique.runMutect2_mutectTag`|String|"mutect2"|Missing?
`mutectRunAllUnique.splitStringToArray_modules`|String|""|Missing?
`mutectRunAllUnique.splitStringToArray_timeout`|Int|1|Missing?
`mutectRunAllUnique.splitStringToArray_memory`|Int|1|Missing?
`mutectRunAllUnique.splitStringToArray_lineSeparator`|String|","|Missing?
`mutectRunAllUnique.normalBam`|File?|None|Input normal file (bam or sam).
`mutectRunAllUnique.normalBai`|File?|None|Missing?
`mutectRunAllUnique.pon`|File?|None|Missing?
`mutectRunAllUnique.ponIdx`|File?|None|Missing?
`mutectRunAllUnique.gnomad`|File?|None|Missing?
`mutectRunAllUnique.gnomadIdx`|File?|None|Missing?
`hsMetricsRunDCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunDCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunDCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunDCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunDCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunDCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunSSCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunSSCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunSSCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunAllUnique.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
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
`variantEffectPredictor.vcf2maf_maxfilterAC`|Int|10|The maximum AC filter
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
`variantEffectPredictor.getSampleNames_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.getSampleNames_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.getSampleNames_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.getSampleNames_modules`|String|"vcftools/0.1.16"|Required environment modules
`variantEffectPredictor.getSampleNames_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.targetBedTask_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.targetBedTask_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name


### Outputs

Output | Type | Description
---|---|---
