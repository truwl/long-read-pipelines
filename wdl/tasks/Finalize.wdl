version 1.0

import "Structs.wdl"

task FinalizeToFile {

    meta {
        description : "Copy a file to an output google cloud location."
        author : "Kiran Garimella"
        email : "kiran@broadinstitute.org"
    }

    input {
        String file
        String outfile

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        file : "Google cloud path to the source file."
        outfile : "Google cloud path to the destination file"

        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."

        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    # This idiom ensures that we don't accidentally have double-slashes in our GCS paths
    String gcs_output_file = sub(sub(outfile, "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        gsutil -m cp ~{file} ~{gcs_output_file}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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

task FinalizeToDir {

    meta {
        description : "Copy a set of files to an output google cloud location."
        author : "Kiran Garimella"
        email : "kiran@broadinstitute.org"
    }

    input {
        Array[String] files
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        files : "Array of Google cloud paths to the source files."
        outdir : "Google cloud path to the destination folder."

        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."

        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    # This idiom ensures that we don't accidentally have double-slashes in our GCS paths
    String gcs_output_dir = sub(sub(outdir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        gsutil -m cp ~{sep=' ' files} ~{gcs_output_dir}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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

task WriteCompletionFile {

    meta {
        description : "Write a file to the given directory indicating the run has completed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String outdir
        File? keyfile
    }

    parameter_meta {
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    command <<<
        set -euxo pipefail

        completion_file="COMPLETED_AT_$(date +%Y%m%dT%H%M%S).txt"
        touch $completion_file

        gsutil cp $completion_file ~{outdir}
    >>>

    #########################

    runtime {
        cpu:                    1
        memory:                 2 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}
