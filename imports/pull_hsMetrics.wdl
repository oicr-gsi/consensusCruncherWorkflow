version 1.0



struct HsResources {
  String refDict
  String refFasta
  String modules
}

workflow hsMetrics {
  
Map[String,HsResources] resources = {
  "hg19": {
    "refDict": "$HG19_ROOT/hg19_random.dict",
    "refFasta": "$HG19_ROOT/hg19_random.fa",
    "modules": "picard/2.19.2 hg19/p13"
  },
  "hg38": {
    "refDict": "$HG38_ROOT/hg38_random.dict",
    "refFasta": "$HG38_ROOT/hg38_random.fa",
    "modules": "picard/2.19.2 hg38/p12"
  }
}

input {
   Int collectHSmetrics_timeout = 5
   Int collectHSmetrics_maxRecordsInRam = 250000
   Int collectHSmetrics_coverageCap = 500
   Int collectHSmetrics_jobMemory = 18
   String collectHSmetrics_filter = "LENIENT"
   String collectHSmetrics_metricTag = "HS"
   Int bedToBaitIntervals_timeout = 1
   Int bedToBaitIntervals_jobMemory = 16
   Int bedToTargetIntervals_timeout = 1
   Int bedToTargetIntervals_jobMemory = 16
   File    inputBam
   String  baitBed
   String  targetBed
   String outputFileNamePrefix = basename(inputBam, '.bam')
   String reference
}

call bedToIntervals as bedToTargetIntervals { 
  input: 
         timeout = bedToTargetIntervals_timeout,
         
         jobMemory = bedToTargetIntervals_jobMemory,
         inputBed = targetBed,
         refDict = resources[reference].refDict,
         modules = resources[reference].modules
    }

call bedToIntervals as bedToBaitIntervals { 
  input: 
         timeout = bedToBaitIntervals_timeout,
         
         jobMemory = bedToBaitIntervals_jobMemory,
         inputBed = baitBed,
         refDict = resources[reference].refDict,
         modules = resources[reference].modules
    }

call collectHSmetrics{ 
  input: 
         timeout = collectHSmetrics_timeout,
         
         maxRecordsInRam = collectHSmetrics_maxRecordsInRam,
         
         coverageCap = collectHSmetrics_coverageCap,
         
         jobMemory = collectHSmetrics_jobMemory,
         
         filter = collectHSmetrics_filter,
         
         metricTag = collectHSmetrics_metricTag,
         inputBam = inputBam, 
         baitIntervals = bedToBaitIntervals.outputIntervals, 
         targetIntervals = bedToTargetIntervals.outputIntervals, 
         outputPrefix = outputFileNamePrefix,
         refFasta = resources[reference].refFasta,
         modules = resources[reference].modules 
         }

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "HSMetrics 2.0"
  dependencies: [{
    name: "picard/2.21.2",
    url: "https://broadinstitute.github.io/picard/"
  }]
  output_meta: {
   outputHSMetrics: "File with HS metrics"
  } 
}

parameter_meta {
    collectHSmetrics_timeout: "Maximum amount of time (in hours) the task can run for."
    collectHSmetrics_maxRecordsInRam: "Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed."
    collectHSmetrics_coverageCap: "Parameter to set a max coverage limit for Theoretical Sensitivity calculations"
    collectHSmetrics_jobMemory: "Memory allocated to job"
    collectHSmetrics_filter: "Settings for picard filter"
    collectHSmetrics_metricTag: "Extension for metrics file"
    bedToBaitIntervals_timeout: "Maximum amount of time (in hours) the task can run for."
    bedToBaitIntervals_jobMemory: "Memory allocated to job"
    bedToTargetIntervals_timeout: "Maximum amount of time (in hours) the task can run for."
    bedToTargetIntervals_jobMemory: "Memory allocated to job"
 inputBam: "Input bam file"
 baitBed: "Path to input bait bed file"
 targetBed: "Path to input target bed"
 outputFileNamePrefix: "Prefix for output"
 reference: "the reference genome for input sample"
}

output {
  File outputHSMetrics  = collectHSmetrics.outputHSMetrics
}

}

task bedToIntervals {
input {
   String inputBed
   String refDict = "$HG19_ROOT/hg19_random.dict"
   Int    jobMemory = 16
   String modules   = "picard/2.21.2 hg19/p13"
   Int timeout = 1
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BedToIntervalList \
                              INPUT=~{inputBed} \
                              OUTPUT="~{basename(inputBed, '.bed')}.interval_list" \
                              SD="~{refDict}"
>>>

parameter_meta {
 inputBed: "Path to input bed file"
 refDict: "Path to index of fasta reference file"
 jobMemory: "Memory allocated to job"
 modules: "Names and versions of modules needed"
 timeout: "Maximum amount of time (in hours) the task can run for."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputIntervals = "~{basename(inputBed, '.bed')}.interval_list"
}
}

task collectHSmetrics {
input { 
   File   inputBam
   String baitIntervals
   String targetIntervals
   String refFasta   = "$HG19_ROOT/hg19_random.fa"
   String metricTag  = "HS"
   String filter     = "LENIENT"
   String outputPrefix = "OUTPUT"
   Int   jobMemory   = 18
   Int   coverageCap = 500
   Int   maxRecordsInRam = 250000
   String modules    = "picard/2.21.2 hg19/p13"
   Int timeout = 5
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar CollectHsMetrics \
                              TMP_DIR=picardTmp \
                              BAIT_INTERVALS=~{baitIntervals} \
                              TARGET_INTERVALS=~{targetIntervals} \
                              R=~{refFasta} \
                              COVERAGE_CAP=~{coverageCap} \
                              MAX_RECORDS_IN_RAM=~{maxRecordsInRam} \
                              INPUT=~{inputBam} \
                              OUTPUT="~{outputPrefix}.~{metricTag}.txt" \
                              VALIDATION_STRINGENCY=~{filter}
>>>

parameter_meta {
 inputBam: "Input bam file"
 baitIntervals: "path to bed file with bait intervals"
 targetIntervals: "path to bed file with target intervals"
 refFasta: "Path to fasta reference file"
 metricTag: "Extension for metrics file"
 filter: "Settings for picard filter"
 outputPrefix: "prefix to build a name for output file"
 coverageCap: "Parameter to set a max coverage limit for Theoretical Sensitivity calculations"
 maxRecordsInRam: "Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed."
 jobMemory: "Memory allocated to job"
 modules: "Names and versions of modules needed"
 timeout: "Maximum amount of time (in hours) the task can run for."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputHSMetrics = "~{outputPrefix}.~{metricTag}.txt"
}

}

