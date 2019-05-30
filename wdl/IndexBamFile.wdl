# Create a bam index for a given bam file.
#
# Description of inputs:
#
#   Required:
#     String gatk_docker                  - GATK Docker image in which to run
#     File input_bam_file                 - BAM file to index
#
#   Optional:
#     gatk4_jar_override                  - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow IndexBamFile {
    String gatk_docker
    File input_bam_file

    File? gatk4_jar_override

    call IndexBamFileTask {
        input:
            gatk_docker    = gatk_docker,
            input_bam_file = input_bam_file,

            gatk_override  = gatk4_jar_override
    }

    output {
        File bam_index     = IndexBamFileTask.bam_index
        File timing_info   = IndexBamFileTask.timing_info
    }
}

#----------------------------------

task IndexBamFileTask {

    # ------------------------------------------------
    # Required Input Args:
    File input_bam_file

    # ------------------------------------------------
    # Required Runtime Args:
    String gatk_docker

    # ------------------------------------------------
    # Optional Runtime Args:
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Process Input args:
    String timing_output_file = basename(input_bam_file) + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 1024 * 3
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 1024

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Run our command:
    command <<<

        set -e

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        samtools index -b ${input_bam_file}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File bam_index   = "${input_bam_file}.bai"
        File timing_info = timing_output_file
    }
}