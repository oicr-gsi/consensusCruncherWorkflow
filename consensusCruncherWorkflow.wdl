version 1.0

import "imports/pull_mutect2.wdl" as mutect2
import "imports/pull_variantEffectPredictor.wdl" as vep
import "imports/pull_hsMetrics.wdl" as hsMetrics

workflow consensusCruncher {
  input {
    File? fastqR1
    File? fastqR2
    File? sortedBam
    File? sortedBai
    String outputFileNamePrefix
    String intervalFile
  }

  parameter_meta {
    fastqR1: "First Fastq Files"
    fastqR2: "Second Fastq File"
    sortedBam: "Bam file from bwamem"
    sortedBai: "Bai file from bwamem"
    outputFileNamePrefix: "Prefix to use for output file"
  }

  if (!(defined(sortedBam)) && defined(fastqR1) && defined(fastqR2)) {
    call align {
      input:
        fastqR1 = select_first([fastqR1]),
        fastqR2 = select_first([fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  call consensus {
    input:
      inputBam = select_first([sortedBam, align.sortedBam]),
      inputBai = select_first([sortedBai, align.sortedBai]),
      basePrefix = outputFileNamePrefix
  }

  call mutect2.mutect2 as mutectRun1 {
    input:
      tumorBam = consensus.dcsScBam,
      tumorBai = consensus.dcsScBamIndex
  }

  call mutect2.mutect2 as mutectRun2 {
    input:
      tumorBam = consensus.sscsScBam,
      tumorBai = consensus.sscsScBamIndex
  }

  call mutect2.mutect2 as mutectRun3 {
    input:
      tumorBam = consensus.allUniqueBam,
      tumorBai = consensus.allUniqueBamIndex
  }

  call hsMetrics.hsMetrics as hsMetricsRun1 {
    input: 
      inputBam = consensus.dcsScBam,
      outputFileNamePrefix = "dcsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile
  }

  call hsMetrics.hsMetrics as hsMetricsRun2 {
    input: 
      inputBam = consensus.dcsScBam,
      outputFileNamePrefix = "sscsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile
  }

  call hsMetrics.hsMetrics as hsMetricsRun3 {
    input: 
      inputBam = consensus.dcsScBam,
      outputFileNamePrefix = "allUnique-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile
  }

  call combineVariants {
    input: 
      inputVcfs = [mutectRun1.filteredVcfFile,mutectRun2.filteredVcfFile],
      inputIndexes = [mutectRun1.filteredVcfIndex,mutectRun2.filteredVcfIndex],
      priority = "mutect2-dcsSc,mutect2-sscsSc",
      outputPrefix = outputFileNamePrefix
  }

  call annotation {
    input: 
      uniqueVcf = mutectRun3.filteredVcfFile,
      uniqueVcfIndex = mutectRun3.filteredVcfIndex,
      mergedVcf = combineVariants.combinedVcf,
      mergedVcfIndex = combineVariants.combinedIndex,
      outputPrefix = outputFileNamePrefix
  }
  
  call vep.variantEffectPredictor {
    input: 
      vcfFile = annotation.annotatedCombinedVcf,
      vcfIndex = annotation.annotatedCombinedIndex,
      toMAF = true,
      onlyTumor = true
  }

  meta {
    author: "Alexander Fortuna and Rishi Shah"
    email: "alexander.fortuna@oicr.on.ca and rshah@oicr.on.ca"
    description: "Workflow to run extract UMIs from fastq and generate consensus Bams as well as run it thru mutect2 and combinevariants"
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
  }

  output {
    File dcsScBam = consensus.dcsScBam
    File dcsScBamIndex = consensus.dcsScBamIndex
    File allUniqueBam = consensus.allUniqueBam
    File allUniqueBamIndex = consensus.allUniqueBamIndex
    File sscsScBam = consensus.sscsScBam
    File sscsScBamIndex = consensus.sscsScBamIndex
    File outputCCStats = consensus.statsCCFile
    File outputCCReadFamilies = consensus.readFamiliesCCFile
    File ccFolder = consensus.ccFolder
    File? mafOutput = variantEffectPredictor.outputMaf
    File dcsScHsMetrics = hsMetricsRun1.outputHSMetrics
    File sscsScHsMetrics = hsMetricsRun2.outputHSMetrics
    File allUniqueHsMetrics = hsMetricsRun3.outputHSMetrics
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
    consensusCruncherPy: "Path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
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
    modules: "Names and versions of modules to load"
    samtools: "Path to samtools binary"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    genome: "Which genome version to use"
    cytoband: "Path to cytoband for genome"
    cutoff: "Cutoff to use to call a consenus of reads"
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

