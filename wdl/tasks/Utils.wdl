version 1.0

import "Structs.wdl"

task ShardLongReads {
    input {
        Array[String] unmapped_files
        Int? num_reads_per_split

        RuntimeAttr? runtime_attr_override
    }

    Int nr = select_first([num_reads_per_split, 200000])
    Int disk_size = 4*ceil(size(unmapped_files, "GB"))
    Int num_files = length(unmapped_files)

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{sep=' -I ' unmapped_files} -nr ~{nr} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        Array[File] unmapped_shards = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ShardLongReadsWithCopy {
    input {
        Array[File] unmapped_files
        Int? num_reads_per_split

        RuntimeAttr? runtime_attr_override
    }

    Int nr = select_first([num_reads_per_split, 200000])
    Int disk_size = 4*ceil(size(unmapped_files, "GB"))
    Int num_files = length(unmapped_files)

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{sep=' -I ' unmapped_files} -nr ~{nr} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        Array[File] unmapped_shards = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task PrepareManifest {
    input {
        Array[String] files

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        echo ~{sep=' ' files} | sed 's/ /\n/g' > manifest.txt
    >>>

    output {
        File manifest = "manifest.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task EchoManifest {
    input {
        File manifest

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cat ~{manifest}
    >>>

    output {
        File out = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ChunkManifest {
    input {
        File manifest
        Int manifest_lines_per_chunk

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        split -a 5 -d --additional-suffix=".txt" -l ~{manifest_lines_per_chunk} ~{manifest} chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("chunk_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# TODO: investigate if samtools is better for this task
# Sort BAM file by coordinate order
# copied from "dsde_pipelines_tasks/BamProcessing.wdl", with
# customization on the runtime block, and "preemptible_tries" taken out
task SortSam {
    input {
        File input_bam
        String output_bam_basename
        Int compression_level

        RuntimeAttr? runtime_attr_override
    }

    # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
    # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

    command {
        java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
            SortSam \
            INPUT=~{input_bam} \
            OUTPUT=~{output_bam_basename}.bam \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=true \
            MAX_RECORDS_IN_RAM=300000 \
            VALIDATION_STRINGENCY=SILENT
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_index = "~{output_bam_basename}.bai"
        File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             5,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FastaToSam {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Float fasta_sam_disk_multiplier = 3.25
    Int disk_size = ceil(fasta_sam_disk_multiplier * size(fasta, "GiB")) + 20

    command <<<
        python /usr/local/bin/prepare_run.py ~{fasta}
    >>>

    output {
        File output_bam = "unmapped.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CountFastqRecords {
    input {
        File fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(fastq, "GiB"))

    command <<<
        set -euxo pipefail

        FILE="~{fastq}"
        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat ~{fastq} | awk '{s++}END{print s/4}'
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat ~{fastq} | awk '{s++}END{print s/4}'
        fi
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CountFastaRecords {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(fasta, "GiB"))

    command <<<
        grep -c '>' ~{fasta}

        exit 0
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CountBamRecords {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))

    command <<<
        samtools view ~{bam} | wc -l
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GrepCountBamRecords {
    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-vc" else "-c"

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GrepCountUniqueBamRecords {
    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-v" else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} | > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Sum {
    input {
        Array[Int] ints
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python -c "print(~{sep="+" ints})" > ~{prefix}.txt
    >>>

    output {
        Int sum = read_int("~{prefix}.txt")
        File sum_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BamToTable {
    input {
        String bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    command <<<
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{bam} | perl -n -e '($nm) = $_ =~ /NM:i:(\d+)/; ($as) = $_ =~ /AS:i:(\d+)/; ($za) = $_ =~ /ZA:Z:(\w+|\.)/; ($zu) = $_ =~ /ZU:Z:(\w+|\.)/; ($cr) = $_ =~ /CR:Z:(\w+|\.)/; ($cb) = $_ =~ /CB:Z:(\w+|\.)/; @a = split(/\s+/); print join("\t", $a[0], $a[1], $a[2], $a[3], $a[4], length($a[9]), $nm, $as, $za, $zu, $cr, $cb, $a[1], ($a[1] & 0x1 ? "paired" : "unpaired"), ($a[1] & 0x4 ? "unmapped" : "mapped"), ($a[1] & 0x10 ? "rev" : "fwd"), ($a[1] & 0x100 ? "secondary" : "primary"), ($a[1] & 0x800 ? "supplementary" : "non_supplementary")) . "\n"' | gzip > ~{prefix}.txt.gz
    >>>

    output {
        File table = "~{prefix}.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConvertReads {
    input {
        File reads
        String output_format
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        filename=~{reads}
        input_filetype=${filename##*.}
        output_filetype=~{output_format}

        if [[ ($input_filetype == "fastq" || $input_filetype == "fq") && $output_filetype == "fasta" ]]; then
            echo "Converting $input_filetype to $output_filetype"
            seqkit fq2fa $filename -o tmp.out
        elif [ $input_filetype == $output_filetype ]; then
            echo "Input filetype is the output filetype"
            mv $filename tmp.out
        else
            echo "ConvertReads does not know how to convert $input_filetype to $output_filetype"
            exit 1
        fi

        mv tmp.out converted_reads.$output_filetype
    >>>

    output {
        File converted_reads = "converted_reads.~{output_format}"
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task BamToBed {
    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail
        bedtools genomecov -ibam ~{bam} -bg > ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# A utility to merge several input BAMs into a single BAM.
task MergeBams {
    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams:   "input array of BAMs to be merged"
        prefix: "[default-valued] prefix for output BAM"
    }

    Int disk_size = 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge -p -c -@2 --no-PG ~{prefix}.bam ~{sep=" " bams}
        samtools index ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeFiles {
    # ------------------------------------------------
    # Input args:
    input {
        Array[File] files_to_merge

        String merged_file_name = "merged_file"

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String timing_output_file = "merge_fastas.timingInformation.txt"
    Int split_size = 5400

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(files_to_merge, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 2048
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command <<<

        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        cat ~{sep=" " files_to_merge} > ~{merged_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Get and compute timing information:
        set +e
        elapsedTime=""
        which bc &> /dev/null ; bcr=$?
        which python3 &> /dev/null ; python3r=$?
        which python &> /dev/null ; pythonr=$?
        if [[ $bcr -eq 0 ]] ; then elapsedTime=`echo "scale=6;$endTime - $startTime" | bc`;
        elif [[ $python3r -eq 0 ]] ; then elapsedTime=`python3 -c "print( $endTime - $startTime )"`;
        elif [[ $pythonr -eq 0 ]] ; then elapsedTime=`python -c "print( $endTime - $startTime )"`;
        fi
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
    >>>

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: select_first([preemptible_attempts, 0])
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
    output {
      File merged_file   = "${merged_file_name}"
      File timing_info    = "${timing_output_file}"
    }
}

# Get the current timestamp as a string.
# Levergaes the unix `date` command.
# You can enter your own date format string.
# The default date string is:
#     %Y%m%d_%H%M%S_%N
# which corresponds to a date of the following form:
# For August 10th, 2020 at 16:06:32.7091 EDT (20:06:32.7091 UTC):
#     20200810_200632_709100000
#
task GetCurrentTimestampString {

    meta {
        # The volatile keyword forces call caching to be turned off, which is
        # exactly what we want for this task.
        # For more info see: https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks/
        volatile: true
    }

    input {
        String date_format = "%Y%m%d_%H%M%S_%N"
    }

    String date_file = "the_date_file.txt"

    command {
        date +~{date_format} > ~{date_file}
        cat ~{date_file}
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: "512 MB"
         disks: "local-disk 10 HDD"
         bootDiskSizeGb: "15"
         preemptible: 0
         cpu: 1
     }

    output {
        String timestamp_string   = read_string(date_file)
    }
}

task ConvertFastaToOneLineSequences {
    # ------------------------------------------------
    # Input args:
    input {
        File reads_fasta

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String out_file_name = sub(reads_fasta, "\\.fasta$", "one_line.fasta")
    String timing_output_file = "convert_fasta_to_one_line_sequences.timingInformation.txt"
    Int split_size = 5400

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_fasta, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 2048
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command <<<

        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ~{reads_fasta} > ~{out_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Get and compute timing information:
        set +e
        elapsedTime=""
        which bc &> /dev/null ; bcr=$?
        which python3 &> /dev/null ; python3r=$?
        which python &> /dev/null ; pythonr=$?
        if [[ $bcr -eq 0 ]] ; then elapsedTime=`echo "scale=6;$endTime - $startTime" | bc`;
        elif [[ $python3r -eq 0 ]] ; then elapsedTime=`python3 -c "print( $endTime - $startTime )"`;
        elif [[ $pythonr -eq 0 ]] ; then elapsedTime=`python -c "print( $endTime - $startTime )"`;
        fi
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
    >>>

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: select_first([preemptible_attempts, 0])
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
    output {
      File one_line_fasta   = "${out_file_name}"
      File timing_info      = "${timing_output_file}"
    }
}
