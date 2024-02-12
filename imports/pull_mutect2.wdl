version 1.0

struct M2Resources {
    String refDict
    String refFai
    String refFasta
    String modules
    String gnomad
    String gnomadIdx
}

workflow mutect2 {
  input {
    Int filter_timeout = 12
    Int filter_memory = 16
    String? filter_filterExtraArgs
    Int mergeStats_timeout = 5
    Int mergeStats_memory = 4
    Int mergeVCFs_timeout = 12
    Int mergeVCFs_memory = 4
    Int runMutect2_timeout = 24
    Int runMutect2_memory = 32
    Int runMutect2_minMemory = 6
    Int runMutect2_threads = 4
    String? runMutect2_mutect2ExtraArgs
    String runMutect2_mutectTag = "mutect2"
    Int getChrCoefficient_timeout = 1
    Int getChrCoefficient_memory = 1
    Int splitStringToArray_timeout = 1
    Int splitStringToArray_memory = 1
    String splitStringToArray_lineSeparator = ","
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    String? intervalFile
    String? intervalsToParallelizeBy
    File? pon
    File? ponIdx
    String reference
    String gatk
  }

  parameter_meta {
      filter_timeout: "Hours before task timeout"
      filter_memory: "Memory allocated for job"
      filter_filterExtraArgs: "placehoulder for extra arguments"
      mergeStats_timeout: "Hours before task timeout"
      mergeStats_memory: "Memory allocated for job"
      mergeVCFs_timeout: "Hours before task timeout"
      mergeVCFs_memory: "Memory allocated for job"
      runMutect2_timeout: "Maximum amount of time (in hours) the task can run for."
      runMutect2_memory: "Memory allocated to job (in GB)."
      runMutect2_minMemory: "Minimum RAM allocated to runMutect2 task (in GB)."
      runMutect2_threads: "Requested CPU threads"
      runMutect2_mutect2ExtraArgs: "placehoulder for extra arguments"
      runMutect2_mutectTag: "version tag for mutect"
      getChrCoefficient_timeout: "Hours before task timeout"
      getChrCoefficient_memory: "Memory allocated for this job"
      splitStringToArray_timeout: "Maximum amount of time (in hours) the task can run for."
      splitStringToArray_memory: "Memory allocated to job (in GB)"
      splitStringToArray_lineSeparator: "Interval group separator - these are the intervals to split by."
    tumorBam: "Input tumor file (bam or sam)."
    tumorBai: "Index for tumorBam"
    normalBam: "Input normal file (bam or sam)."
    normalBai: "Index for noramlBam"
    intervalFile: "One or more genomic intervals over which to operate"
    intervalsToParallelizeBy: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)"
    pon: "panel of normal"
    ponIdx: "index of pon"
    gatk: "gatk version to be used"
    reference: "the reference genome for input sample"
  }

  meta {
    author: "Angie Mosquera, Alexander Fortuna"
    email: "amosquera@oicr.on.ca, afortuna@oicr.on.ca"
    description: "Somatic short variant analysis."
    dependencies: [
    {
      name: "samtools/1.9",
      url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
    }]
    output_meta: {
      filteredVcfFile: "the filtered vcf file",
      filteredVcfIndex: "Index of filtered vcf file",
      mergedUnfilteredStats: "Stats for merged unfiltered files",
      filteringStats: "Stats for filtering process"
    }
  }

Map[String, M2Resources] resources = {
  "hg19": {
        "refDict" : "$HG19_ROOT/hg19_random.dict",
    		"refFai" : "$HG19_ROOT/hg19_random.fa.fai",
    		"refFasta" : "$HG19_ROOT/hg19_random.fa",
    		"modules" : "hg19/p13 samtools/1.9",
        "gnomad": "",
        "gnomadIdx": ""
  },
  "hg38": {
        "refDict" : "$HG38_ROOT/hg38_random.dict",
    		"refFai" : "$HG38_ROOT/hg38_random.fa.fai",
    		"refFasta" : "$HG38_ROOT/hg38_random.fa",
        "gnomad": "$HG38_GATK_GNOMAD_ROOT/af-only-gnomad.hg38.vcf.gz",
        "gnomadIdx": "$HG38_GATK_GNOMAD_ROOT/af-only-gnomad.hg38.vcf.gz.tbi",
    		"modules" : "hg38/p12 samtools/1.9 hg38-gatk-gnomad/2.0"
  },
  "mm10": {
        "refDict" : "$MM10_ROOT/mm10.dict",
        "refFai" : "$MM10_ROOT/mm10.fa.fai",
        "refFasta" : "$MM10_ROOT/mm10.fa",
        "modules" : "mm10/p6 samtools/1.9",
        "gnomad": "",
        "gnomadIdx": ""
  }

}

  call splitStringToArray {
    input:
      timeout = splitStringToArray_timeout,
      memory = splitStringToArray_memory,
      lineSeparator = splitStringToArray_lineSeparator,
      intervalsToParallelizeBy = intervalsToParallelizeBy
  }

  String outputBasename = basename(tumorBam, '.bam')
  Boolean intervalsProvided = if (defined(intervalsToParallelizeBy)) then true else false

  scatter(subinterval in flatten(splitStringToArray.out)) {
    call getChrCoefficient {
      input:
        timeout = getChrCoefficient_timeout,
        memory = getChrCoefficient_memory,
        refDict = resources[reference].refDict,
        region = subinterval,
        modules = resources [ reference ].modules
    }

    call runMutect2 {
      input:
        timeout = runMutect2_timeout,
        memory = runMutect2_memory,
        minMemory = runMutect2_minMemory,
        threads = runMutect2_threads,
        mutect2ExtraArgs = runMutect2_mutect2ExtraArgs,
        mutectTag = runMutect2_mutectTag,
        intervals = [subinterval],
        intervalsProvided = intervalsProvided,
        intervalFile = intervalFile,
        tumorBam = tumorBam,
        tumorBai = tumorBai,
        normalBam = normalBam,
        normalBai = normalBai,
        pon = pon,
        ponIdx = ponIdx,
        gnomad = resources [ reference ].gnomad,
        gnomadIdx = resources [ reference ].gnomadIdx,
        outputBasename = outputBasename,
        modules = resources [ reference ].modules + ' ' + gatk,
        refFai = resources[reference].refFai,
        refFasta = resources[reference].refFasta,
        refDict = resources[reference].refDict,
        scaleCoefficient = getChrCoefficient.coeff

    }
  }

  Array[File] unfilteredVcfs = runMutect2.unfilteredVcf
  Array[File] unfilteredVcfIndices = runMutect2.unfilteredVcfIdx
  Array[File] unfilteredStats = runMutect2.stats

  call mergeVCFs {
    input:
      timeout = mergeVCFs_timeout,
      memory = mergeVCFs_memory,
      vcfs = unfilteredVcfs,
      vcfIndices = unfilteredVcfIndices,
      modules = resources [ reference ].modules + ' ' + gatk,
      refFasta = resources[reference].refFasta
  }

  call mergeStats {
    input:
      timeout = mergeStats_timeout,
      memory = mergeStats_memory,
      stats = unfilteredStats,
      modules = resources [ reference ].modules + ' ' + gatk,
  }

  call filter {
    input:
      timeout = filter_timeout,
      memory = filter_memory,
      filterExtraArgs = filter_filterExtraArgs,
      intervalFile = intervalFile,
      unfilteredVcf = mergeVCFs.mergedVcf,
      unfilteredVcfIdx = mergeVCFs.mergedVcfIdx,
      mutectStats = mergeStats.mergedStats,
      modules = resources [ reference ].modules + ' ' + gatk,
      refFasta = resources[reference].refFasta,
      refDict = resources[reference].refDict,
      refFai = resources[reference].refFai
  }


  output {
    File filteredVcfFile = filter.filteredVcfGz
    File filteredVcfIndex = filter.filteredVcfTbi
    File mergedUnfilteredStats = mergeStats.mergedStats
    File filteringStats = filter.filteringStats
  }
}

task splitStringToArray {
  input {
    String? intervalsToParallelizeBy
    String lineSeparator = ","
    Int memory = 1
    Int timeout = 1
  }
  parameter_meta {
    lineSeparator: "Interval group separator - these are the intervals to split by."
    memory: "Memory allocated to job (in GB)"
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
    echo "~{intervalsToParallelizeBy}" | tr '~{lineSeparator}' '\n'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }
}

# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 1
    Int timeout = 1
    String modules
    String region
    String refDict
  }

  parameter_meta {
    refDict: ".dict file for the reference genome, we use it to extract chromosome ids"
    timeout: "Hours before task timeout"
    region: "Region to extract a chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Environment module names and version to load (space separated) before command execution"
  }

  command <<<
    CHROM=$(echo ~{region} | sed 's/:.*//')
    if [[ $CHROM ]]; then
      LARGEST=$(grep SN:chr ~{refDict} | cut -f 3 | sed 's/LN://' | sort -n | tail -n 1)
      grep -w SN:$CHROM ~{refDict} | cut -f 3 | sed 's/.*://' | awk -v largest_chr=$LARGEST '{print int(($1/largest_chr + 0.1) * 10)/10}'
    else
      echo "1.0"
    fi
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
  }
}

task runMutect2 {
  input {
    String modules
    String refFasta
    String refFai
    String refDict 
    String mutectTag = "mutect2"
    String? intervalFile
    Array[String]? intervals
    Boolean intervalsProvided
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    File? pon
    File? ponIdx
    String gnomad
    String gnomadIdx
    String? mutect2ExtraArgs
    String outputBasename
    Int threads = 4
    Int memory = 32
    Int minMemory = 6
    Int timeout = 24
    Float scaleCoefficient = 1.0
  }

  parameter_meta {
    modules: "Modules to be used by the task"    
    refFasta: "Path to the reference fasta"
    refFai: "Reference fasta fai index"
    refDict: "Reference fasta dict file"
    intervalFile: "Interval file, optional" 
    intervals: "Normally, this is a chromosome. This task is parallelized"
    intervalsProvided: "Flag that shows that we have intervals provided"
    tumorBam: "Tumor bam file path" 
    tumorBai: "Tumor bai index"
    normalBam: "Normal bam" 
    normalBai: "Normal bai"
    pon: "Panel of normals"
    ponIdx: "Index for PON"
    gnomad: "Gnomad resource, reference germline variants"
    gnomadIdx: "Index for GNOMAD"
    outputBasename: "basename for output"
    mutectTag: "version tag for mutect"
    mutect2ExtraArgs: "placehoulder for extra arguments"
    threads: "Requested CPU threads"
    memory: "Memory allocated to job (in GB)."
    minMemory: "A minimum amount of memory allocated to the task, overrides the scaled RAM setting"
    timeout: "Maximum amount of time (in hours) the task can run for."
    scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
  }

  String outputVcf = if (defined(normalBam)) then outputBasename + "." + mutectTag + ".vcf" else outputBasename + "." + mutectTag + ".tumor_only.vcf"
  String outputVcfIdx = outputVcf + ".idx"
  String outputStats = outputVcf + ".stats"
  Int allocatedMemory = if minMemory > round(memory * scaleCoefficient) then minMemory else round(memory * scaleCoefficient)

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{tumorBam} -O tumor_name.txt -encode
    tumor_command_line="-I ~{tumorBam} -tumor `cat tumor_name.txt`"

    cp ~{refFai} .
    cp ~{refDict} .

    if [ -f "~{normalBam}" ]; then
      gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{normalBam} -O normal_name.txt -encode
      normal_command_line="-I ~{normalBam} -normal `cat normal_name.txt`"
    else
      normal_command_line=""
    fi

    if [ -f "~{intervalFile}" ]; then
      if ~{intervalsProvided} ; then
        intervals_command_line="-L ~{sep=" -L " intervals} -L ~{intervalFile} -isr INTERSECTION"
      else
        intervals_command_line="-L ~{intervalFile}"
      fi
    else
      if ~{intervalsProvided} ; then
        intervals_command_line="-L ~{sep=" -L " intervals} "
      fi
    fi

    if [[ ! -z "~{gnomad}" ]] ; then
      germline_resource_line="--germline-resource ~{gnomad}"
    else 
      germline_resource_line=""
    fi

    gatk --java-options "-Xmx~{memory-8}g" Mutect2 \
    -R ~{refFasta} \
    $tumor_command_line \
    $normal_command_line \
    $germline_resource_line \
    ~{"-pon " + pon} \
    $intervals_command_line \
    -O "~{outputVcf}" \
    ~{mutect2ExtraArgs}
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{allocatedMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}"
    File unfilteredVcfIdx = "~{outputVcfIdx}"
    File stats = "~{outputStats}"
  }
}

task mergeVCFs {
  input {
    String modules
    String refFasta
    Array[File] vcfs
    Array[File] vcfIndices
    Int memory = 4
    Int timeout = 12
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    vcfs: "Vcf's from scatter to merge together"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      mergedVcf: "Merged vcf, unfiltered.",
      mergedVcfIdx: "Merged vcf index, unfiltered."
    }
  }

  String outputName = basename(vcfs[0])

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{memory-3}g" MergeVcfs \
    -I ~{sep=" -I " vcfs} \
    -O ~{outputName}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "~{outputName}"
    File mergedVcfIdx = "~{outputName}.idx"
  }
}

task mergeStats {
  input {
    String modules
    Array[File]+ stats
    Int memory = 4
    Int timeout = 5
  }

  parameter_meta {
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  String outputStats = basename(stats[0])

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{memory-3}g" MergeMutectStats \
    -stats ~{sep=" -stats " stats} \
    -O ~{outputStats}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedStats = "~{outputStats}"
  }
}

task filter {
  input {
    String modules
    String refFasta
    String refFai
    String refDict
    String? intervalFile
    File unfilteredVcf
    File unfilteredVcfIdx
    File mutectStats
    String? filterExtraArgs
    Int memory = 16
    Int timeout = 12
  }

  parameter_meta {
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
    filterExtraArgs: "placehoulder for extra arguments"
  }

  String unfilteredVcfName = basename(unfilteredVcf)
  String filteredVcfName = basename(unfilteredVcf, ".vcf") + ".filtered.vcf"

  command <<<
    set -euo pipefail

    cp ~{refFai} .
    cp ~{refDict} .

    gatk --java-options "-Xmx~{memory-4}g" FilterMutectCalls \
    -V ~{unfilteredVcf} \
    -R ~{refFasta} \
    -O ~{filteredVcfName} \
    ~{"-stats " + mutectStats} \
    --filtering-stats ~{filteredVcfName}.stats \
    ~{filterExtraArgs}

    bgzip -c ~{filteredVcfName} > ~{filteredVcfName}.gz
    bgzip -c ~{unfilteredVcf} > ~{unfilteredVcfName}.gz

    gatk --java-options "-Xmx~{memory-5}g" IndexFeatureFile -I ~{filteredVcfName}.gz
    gatk --java-options "-Xmx~{memory-5}g" IndexFeatureFile -I ~{unfilteredVcfName}.gz
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcfGz = "~{unfilteredVcfName}.gz"
    File unfilteredVcfTbi = "~{unfilteredVcfName}.gz.tbi"
    File filteredVcfGz = "~{filteredVcfName}.gz"
    File filteredVcfTbi = "~{filteredVcfName}.gz.tbi"
    File filteringStats = "~{filteredVcfName}.stats"
  }
}
