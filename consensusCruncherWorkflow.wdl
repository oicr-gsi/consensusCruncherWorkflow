version 1.0

import "imports/pull_mutect2.wdl" as mutect2
import "imports/pull_variantEffectPredictor.wdl" as vep
import "imports/pull_hsMetrics.wdl" as hsMetrics

struct InputGroup {
  File fastqR1
  File fastqR2
}

struct RefResources {
  String align_modules
  String align_bwaref
  String align_blist
  String consensus_modules
  String consensus_genome
  String consensus_cytoband
  String combineVariants_modules
  String inputRefFasta
}

workflow consensusCruncher {
  input {
    Array[InputGroup]? inputGroups
    File? sortedBam
    File? sortedBai
    String outputFileNamePrefix
    String intervalFile
    String inputIntervalsToParalellizeBy
    String tumorName
    String reference
  }

  Map[String,RefResources] resources = {
    "hg19": {
      "align_modules": "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9",
      "align_bwaref": "$HG19_BWA_INDEX_ROOT/hg19_random.fa",
      "align_blist": "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list",
      "consensus_modules": "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9",
      "consensus_genome": "hg19",
      "consensus_cytoband": "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg19/p13",
      "inputRefFasta": "$HG19_ROOT/hg19_random.fa"
    },
    "hg38": {
      "align_modules": "consensus-cruncher/5.0 data-hg38-consensus-cruncher/1.0 hg38-bwa-index-with-alt/0.7.12 samtools/1.9",
      "align_bwaref": "$HG38_BWA_INDEX_WITH_ALT_ROOT/hg38_random.fa",
      "align_blist": "$DATA_HG38_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list",
      "consensus_modules": "consensus-cruncher/5.0 data-hg38-consensus-cruncher/1.0 hg38-bwa-index-with-alt/0.7.12 samtools/1.9",
      "consensus_genome": "hg38",
      "consensus_cytoband": "$DATA_HG38_CONSENSUS_CRUNCHER_ROOT/hg38_cytoBand.txt",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg38/p12",
      "inputRefFasta": "$HG38_ROOT/hg38_random.fa"
    }
  }
  
  parameter_meta {
    inputGroups: "Array of fastq files to concatenate if a top-up"
    sortedBam: "Bam file from bwamem"
    sortedBai: "Bai file from bwamem"
    outputFileNamePrefix: "Prefix to use for output file"
    intervalFile: "interval file to subset variant calls"
    inputIntervalsToParalellizeBy: "intervals for parallelization"
    tumorName: "Name of the tumor sample"
    reference: "reference version"
  }

  if (!(defined(sortedBam)) && defined(inputGroups)) {
    Array[InputGroup] inputs = select_first([inputGroups])
    scatter (ig in inputs) {
      File read1s       = ig.fastqR1
      File read2s       = ig.fastqR2
    }
  }
  
  if (!(defined(sortedBam)) && defined(inputGroups)) {
    call concat {
      input:
        read1s = select_first([read1s]),
        read2s = select_first([read2s]),
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  if (!(defined(sortedBam)) && defined(read1s) && defined(read2s)) {
    call align {
      input:
        fastqR1 = select_first([concat.fastqR1]),
        fastqR2 = select_first([concat.fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix,
        modules = resources[reference].align_modules,
        bwaref = resources[reference].align_bwaref,
        blist = resources[reference].align_blist
    }
  }

  call consensus {
    input:
      inputBam = select_first([sortedBam, align.sortedBam]),
      inputBai = select_first([sortedBai, align.sortedBai]),
      basePrefix = outputFileNamePrefix,
      modules = resources[reference].consensus_modules,
      cytoband = resources[reference].consensus_cytoband,
      genome = resources[reference].consensus_genome
  }

  call mutect2.mutect2 as mutectRunDCSSC {
    input:
      tumorBam = consensus.dcsScBam,
      tumorBai = consensus.dcsScBamIndex,
      reference = reference,
      intervalFile = intervalFile,
      intervalsToParallelizeBy = inputIntervalsToParalellizeBy
  }

  call mutect2.mutect2 as mutectRunSSCSSC {
    input:
      tumorBam = consensus.sscsScBam,
      tumorBai = consensus.sscsScBamIndex,
      reference = reference,
      intervalFile = intervalFile,
      intervalsToParallelizeBy = inputIntervalsToParalellizeBy
  }

  call mutect2.mutect2 as mutectRunAllUnique {
    input:
      tumorBam = consensus.allUniqueBam,
      tumorBai = consensus.allUniqueBamIndex,
      reference = reference,
      intervalFile = intervalFile,
      intervalsToParallelizeBy = inputIntervalsToParalellizeBy
  }

  call hsMetrics.hsMetrics as hsMetricsRunDCSSC {
    input: 
      inputBam = consensus.dcsScBam,
      outputFileNamePrefix = "dcsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile,
      reference = reference
  }

  call hsMetrics.hsMetrics as hsMetricsRunSSCSSC {
    input: 
      inputBam = consensus.sscsScBam,
      outputFileNamePrefix = "sscsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile
  }

  call hsMetrics.hsMetrics as hsMetricsRunAllUnique {
    input: 
      inputBam = consensus.allUniqueBam,
      outputFileNamePrefix = "allUnique-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile
  }

  call combineVariants {
    input: 
      inputVcfs = [mutectRunDCSSC.filteredVcfFile,mutectRunSSCSSC.filteredVcfFile],
      inputIndexes = [mutectRunDCSSC.filteredVcfIndex,mutectRunSSCSSC.filteredVcfIndex],
      priority = "mutect2-dcsSc,mutect2-sscsSc",
      outputPrefix = outputFileNamePrefix,
      referenceFasta = resources[reference].inputRefFasta,
      modules = resources[reference].combineVariants_modules
  }

  call annotation {
    input: 
      uniqueVcf = mutectRunAllUnique.filteredVcfFile,
      uniqueVcfIndex = mutectRunAllUnique.filteredVcfIndex,
      mergedVcf = combineVariants.combinedVcf,
      mergedVcfIndex = combineVariants.combinedIndex,
      outputPrefix = outputFileNamePrefix
  }
  
  call vep.variantEffectPredictor {
    input: 
      vcfFile = annotation.annotatedCombinedVcf,
      vcfIndex = annotation.annotatedCombinedIndex,
      toMAF = true,
      onlyTumor = true,
      tumorOnlyAlign_updateTagValue = true,
      vcf2maf_retainInfoProvided = true,
      reference = reference,
      targetBed = intervalFile,
      tumorName = tumorName
  }

  meta {
    author: "Alexander Fortuna and Rishi Shah"
    email: "alexander.fortuna@oicr.on.ca and rshah@oicr.on.ca"
    description: "ConsensusCruncher is a tool developed by Pugh lab. It suppresses errors in next-generation sequencing data by using unique molecular identifers (UMIs) to amalgamate reads derived from the same DNA template into a consensus sequence. The workflow also uses the generated consensus Bams to run it thru mutect2 and combinevariants tasks."
    dependencies: [
     {
      name: "hg19-bwa-index/0.7.12",
      url: "http://bio-bwa.sourceforge.net/"
     },
     {
      name: "samtools/1.9",
      url: "http://www.htslib.org/"
     },
     {
      name: "python/3.6",
      url: "https://www.python.org/downloads/"
     },
     {
      name: "picard/2.21.2",
      url: "https://broadinstitute.github.io/picard/"
     },
     {
      name: "rstats/3.6",
      url: "https://www.r-project.org/"
     },
     {
      name: "consensuscruncer-5.0",
      url: "https://github.com/pughlab/ConsensusCruncher"
     }
    ]
    output_meta: {
    rawBam: {
        description: "raw bam file",
        vidarr_label: "rawBam"
    },
    rawBamIndex: {
        description: "raw bam file index",
        vidarr_label: "rawBamIndex"
    },
    dcsScBam: {
        description: "Duplex Consensus Sequence (DCS) generated from Single-Strand Consensus Sequences (SSCS) + Single Correction (SC)",
        vidarr_label: "dcsScBam"
    },
    dcsScBamIndex: {
        description: "Bam index for DCS generated from SSCS + SC",
        vidarr_label: "dcsScBamIndex"
    },
    allUniqueBam: {
        description: "DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons",
        vidarr_label: "allUniqueBam"
    },
    allUniqueBamIndex: {
        description: "Bam index for DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons",
        vidarr_label: "allUniqueBamIndex"
    },
    sscsScBam: {
        description: "SSCS combined with corrected singletons",
        vidarr_label: "sscsScBam"
    },
    sscsScBamIndex: {
        description: "Bam index for SSCS combined with corrected singletons",
        vidarr_label: "sscsScBamIndex"
    },
    outputCCStats: {
        description: "CC stats",
        vidarr_label: "outputCCStats"
    },
    outputCCReadFamilies: {
        description: "CC read families",
        vidarr_label: "outputCCReadFamilies"
    },
    ccFolder: {
        description: "cc folder",
        vidarr_label: "ccFolder"
    },
    dcsScVcf: {
        description: "DCS vcf",
        vidarr_label: "dcsScVcf"
    },
    dcsScVcfIndex: {
        description: "DCS vcf index",
        vidarr_label: "dcsScVcfIndex"
    },
    allUniqueVcf: {
        description: "vcf of DCS + singletons",
        vidarr_label: "allUniqueVcf"
    },
    allUniqueVcfIndex: {
        description: "vcf index for DCS + singletons",
        vidarr_label: "allUniqueVcfIndex"
    },
    sscsScVcf: {
        description: "SSCS vcf",
        vidarr_label: "sscsScVcf"
    },
    sscsScVcfIndex: {
        description: "SSCS vcf index",
        vidarr_label: "sscsScVcfIndex"
    },
    vepVcf: {
        description: "vep vcf",
        vidarr_label: "vepVcf"
    },
    vepVcfIndex: {
        description: "vep vcf index",
        vidarr_label: "vepVcfIndex"
    },
    mafOutput: {
        description: "maf output",
        vidarr_label: "mafOutput"
    },
    dcsScHsMetrics: {
        description: "DCS metrics",
        vidarr_label: "dcsScHsMetrics"
    },
    sscsScHsMetrics: {
        description: "SSCS metrics",
        vidarr_label: "sscsScHsMetrics"
    },
    allUniqueHsMetrics: {
        description: "metrics for DCS + singletons",
        vidarr_label: "allUniqueHsMetrics"
    }
}
  }
  
  output {
    File? rawBam = align.sortedBam
    File? rawBamIndex = align.sortedBai
    File dcsScBam = consensus.dcsScBam
    File dcsScBamIndex = consensus.dcsScBamIndex
    File allUniqueBam = consensus.allUniqueBam
    File allUniqueBamIndex = consensus.allUniqueBamIndex
    File sscsScBam = consensus.sscsScBam
    File sscsScBamIndex = consensus.sscsScBamIndex
    File outputCCStats = consensus.statsCCFile
    File outputCCReadFamilies = consensus.readFamiliesCCFile
    File ccFolder = consensus.ccFolder
    File dcsScVcf = mutectRunDCSSC.filteredVcfFile
    File dcsScVcfIndex = mutectRunDCSSC.filteredVcfIndex
    File allUniqueVcf = mutectRunAllUnique.filteredVcfFile
    File allUniqueVcfIndex = mutectRunAllUnique.filteredVcfIndex
    File sscsScVcf = mutectRunSSCSSC.filteredVcfFile
    File sscsScVcfIndex = mutectRunSSCSSC.filteredVcfIndex
    File vepVcf = variantEffectPredictor.outputVcf
    File vepVcfIndex = variantEffectPredictor.outputTbi
    File? mafOutput = variantEffectPredictor.outputMaf
    File dcsScHsMetrics = hsMetricsRunDCSSC.outputHSMetrics
    File sscsScHsMetrics = hsMetricsRunSSCSSC.outputHSMetrics
    File allUniqueHsMetrics = hsMetricsRunAllUnique.outputHSMetrics
  }
}

task concat {
  input {
    Array[File]+ read1s
    Array[File]+ read2s
    String outputFileNamePrefix
    Int threads = 4
    Int jobMemory = 16
    Int timeout = 72
    String modules = "tabix/0.2.6"
  }

  parameter_meta {
    read1s: "array of read1s"
    read2s: "array of read2s"
    outputFileNamePrefix: "File name prefix"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    modules: "Required environment modules"
  }

  command <<<
    set -euo pipefail

    zcat ~{sep=" " read1s} | gzip > ~{outputFileNamePrefix}_R1_001.fastq.gz

    zcat ~{sep=" " read2s} | gzip > ~{outputFileNamePrefix}_R2_001.fastq.gz

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
    modules: "~{modules}"

  }

  output {
    File fastqR1 = "~{outputFileNamePrefix}_R1_001.fastq.gz"
    File fastqR2 = "~{outputFileNamePrefix}_R2_001.fastq.gz"
  }
}

task align {
  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
    String modules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String bwa = "$BWA_ROOT/bin/bwa"
    String bwaref = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String blist = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"
    Int threads = 4
    Int jobMemory = 16
    Int timeout = 72
  }

  parameter_meta {
    fastqR1: "Path to left fastq file"
    fastqR2: "Path to right fastq file"
    outputFileNamePrefix: "File name prefix"
    modules: "Names and versions of modules to load"
    consensusCruncherPy: "Path to consensusCruncher binary"
    bwa: "Path to bwa binary"
    bwaref: "Path to bwa index"
    samtools: "Path to samtools binary"
    blist: "Path to blacklist for barcodes"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    ~{consensusCruncherPy} fastq2bam \
         --fastq1 ~{fastqR1} \
         --fastq2 ~{fastqR2}\
         --output . \
         --bwa ~{bwa} \
         --ref ~{bwaref} \
         --samtools ~{samtools} \
         --skipcheck \
         --blist ~{blist}

    # Necessary for if bam files to be named according to merged library name
    # Additionally if ".sorted" isn't omitted here, file names from align include ".sorted" twice
    mv bamfiles/*.bam bamfiles/"~{outputFileNamePrefix}.bam"
    mv bamfiles/*.bai bamfiles/"~{outputFileNamePrefix}.bam.bai"
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"

  }

  output {
    File? sortedBam = "bamfiles/~{outputFileNamePrefix}.bam"
    File? sortedBai = "bamfiles/~{outputFileNamePrefix}.bam.bai"
  }
}

task consensus {
  input {
    File? inputBam
    File? inputBai
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String basePrefix
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String cytoband = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt"
    String genome   = "hg19"
    String ccDir = basePrefix + ".consensuscruncher"
    Float cutoff  = 0.7
    Int threads = 8
    Int jobMemory = 32
    Int timeout = 72
    String modules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"
  }

  parameter_meta {
    inputBam: "Bam file either from bwamem or ConsensusCruncher align."
    inputBai: "Bai file either from bwamem or ConsensusCruncher align."
    consensusCruncherPy: "Path to consensusCruncher binary"
    basePrefix: "Prefix identifier for naming of output files"
    samtools: "Path to samtools binary"
    cytoband: "Path to cytoband for genome"
    genome: "Which genome version to use"
    ccDir: "Placeholder"
    cutoff: "Cutoff to use to call a consenus of reads"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    modules: "Names and versions of modules to load"
  }


  command <<<
  set -euo pipefail

   ~{consensusCruncherPy} consensus \
         --input ~{inputBam} \
         --output . \
         --samtools ~{samtools} \
         --cutoff ~{cutoff} \
         --genome ~{genome} \
         --bedfile ~{cytoband} \
         --bdelim '|'

   tar cf - ~{basePrefix} | gzip --no-name > ~{ccDir}.tar.gz
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File dcsScBam = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam"
    File dcsScBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam.bai"
    File allUniqueBam = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam"
    File allUniqueBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam.bai"
    File sscsScBam = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam"
    File sscsScBamIndex = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam.bai"
    File statsCCFile = "~{basePrefix}/~{basePrefix}.stats.txt"
    File readFamiliesCCFile = "~{basePrefix}/~{basePrefix}.read_families.txt"
    File ccFolder = "~{ccDir}.tar.gz"
    
  }

  meta {
    output_meta: {
      dcsScBam: "DCS generated from SSCS + SC",
      dcsScBamIndex: "Index for DCS SC Bam",
      allUniqueBam: "DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons",
      allUniqueBamIndex: "Index for All Unique Bam",
      sscsScBam: "SSCS combined with corrected singletons (from both rescue strategies)",
      sscsScBamIndex: "Index for SSCS SC Bam",
      statsCCFile: "Stats file",
      readFamiliesCCFile: "read families file",
      ccFolder: "output folder containing files not needed for downstream analysis; info on family size, QC metrics"
    }
  }
}

task combineVariants {
input {
 Array[File] inputVcfs
 Array[File] inputIndexes
 Array[String] workflows
 String referenceFasta
 String outputPrefix 
 String modules
 String priority
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 inputIndexes: "array of tabix indexes for vcf files"
 workflows: "array of ids of producer workflows"
 referenceFasta: "path to the reference FASTA file"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 priority: "Comma-separated list defining priority of workflows when combining variants"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
  python3<<CODE
  import subprocess
  import sys
  inputStrings = []
  v = "~{sep=' ' inputVcfs}"
  vcfFiles = v.split()
  w = "~{sep=' ' workflows}"
  workflowIds = w.split()
  priority = "~{priority}"
  
  if len(vcfFiles) != len(workflowIds):
      print("The arrays with input files and their respective workflow names are not of equal size!")
  else:
      for f in range(0, len(vcfFiles)):
          inputStrings.append("--variant:" + workflowIds[f] + " " + vcfFiles[f])

  javaMemory = ~{jobMemory} - 6 
  gatkCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar $GATK_ROOT/GenomeAnalysisTK.jar "
  gatkCommand += "-T CombineVariants "
  gatkCommand += " ".join(inputStrings)
  gatkCommand += " -R ~{referenceFasta} "
  gatkCommand += "-o ~{outputPrefix}_combined.vcf.gz "
  gatkCommand += "-genotypeMergeOptions PRIORITIZE "
  gatkCommand += "-priority " + priority
  gatkCommand += " 2>&1"

  result_output = subprocess.run(gatkCommand, shell=True)
  sys.exit(result_output.returncode)
  CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File combinedVcf = "~{outputPrefix}_combined.vcf.gz"
  File combinedIndex = "~{outputPrefix}_combined.vcf.gz.tbi"
}
}

task annotation {
input {
 File uniqueVcf 
 File uniqueVcfIndex
 File mergedVcf
 File mergedVcfIndex
 String outputPrefix 
 String modules = "samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 uniqueVcf: "input unique vcf files"
 uniqueVcfIndex: "input unique tabix indexes for vcf files"
 mergedVcf: "input merged vcf"
 mergedVcfIndex: "input merged vcf index"
 outputPrefix: "prefix for output file"
 modules: "module for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
  bcftools annotate -a ~{uniqueVcf} \
 -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
 -o "~{outputPrefix}.merged.vcf.gz"

 tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File annotatedCombinedVcf = "~{outputPrefix}.merged.vcf.gz"
  File annotatedCombinedIndex = "~{outputPrefix}.merged.vcf.gz.tbi"
}
}
