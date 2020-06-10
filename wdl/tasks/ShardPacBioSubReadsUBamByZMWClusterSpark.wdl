version 1.0

import "Structs.wdl"

workflow ShardPacBioSubReadsUBamByZMWClusterSpark {
	input {
		File input_ubam
		String output_prefix

        Int? shard_size
	}

    call CreateSparkIndex {
        input:
        input_ubam = input_ubam
    }

	call SparkShard {
		input:
		input_ubam = input_ubam,
        input_ubam_splittingindex = CreateSparkIndex.spark_sbi_index,
		split_prefix = output_prefix,
        shard_size = shard_size
	}

    output {
        File spark_sbi_index = CreateSparkIndex.spark_sbi_index

        Array[File] split_bam = SparkShard.split_bam
        File original_header_hd_line = SparkShard.original_header_hd_line
    }
}

task CreateSparkIndex {
    input {
        File input_ubam

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_ubam, ".bam")

    Int disk_size = 10 + ceil(size(input_ubam, "GB")) # assumes the sbi index will be large

    command <<<
        set -euo pipefail

        gatk \
            --java-options "-Xms35G -Xmx45G" \
            CreateHadoopBamSplittingIndex \
            -I ~{input_ubam} \
            -O ~{prefix}.sbi \
            --splitting-index-granularity 1
    >>>

    output {
        File spark_sbi_index = "~{prefix}.sbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "shuangbroad/gatk:custom"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SparkShard {
	input {
        File input_ubam
        File input_ubam_splittingindex

        Int? shard_size

        String split_prefix

        RuntimeAttr? runtime_attr_override
    }

    String sparse_zmws_per_shards = if (defined(shard_size)) then "--shard-size ~{shard_size}" else " "

    Int cpu_guess = 96
    Int disk_size = if ( 600 < ceil(size(input_ubam, "GB")) ) then 3000 else 2250 # 8 vs 6 local SSDs

    Int mem_guess = if ( 600 < ceil(size(input_ubam, "GB")) ) then 600 else 400

    Int jvm_start_mem = ceil (0.85 * mem_guess)
    Int jvm_max_mem   = ceil (0.90 * mem_guess)
    String java_options = "\"-Xms" + jvm_start_mem + "G -Xmx" + jvm_max_mem + "G\""

    Int spark_driver_mem = ceil (0.87 * mem_guess)

    command <<<
        set -euo pipefail

        mkdir -p split_dir
        gatk \
            --java-options ~{java_options} \
            ShardPacBioSubReadsUBamByZMWClusterSpark \
            -I ~{input_ubam} \
            --read-index ~{input_ubam_splittingindex} \
            -O split_dir/~{split_prefix} \
            ~{sparse_zmws_per_shards} \
            --use-jdk-deflater \
            --use-jdk-inflater \
            -- \
            --conf spark.master="local[*]" \
            --conf spark.driver.memory=~{spark_driver_mem}g \
            --conf spark.memory.fraction=0.85 \
            --conf spark.memory.storageFraction=0.25

        # some header hack, otherwise downstream CCS will fail
        # because the Spark tool rips off the extra (illegal?) "pb:3.0.5" from the @HD line
        samtools view -H ~{input_ubam} | head -n 1 | sed -E 's/SO:unknown/SO:queryname/'> save_hd_line.txt

        tot_bam_cnt=$(ls split_dir/*.bam | wc -l)
        echo -e "\nTotal number of bams: ${tot_bam_cnt}\n\nDisk space usage:"
        df -h .
        echo
    >>>

    output {
        Array[File] split_bam = glob("split_dir/~{split_prefix}*.bam")
        File original_header_hd_line = "save_hd_line.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu_guess,
        mem_gb:             mem_guess,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "shuangbroad/gatk:custom"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
