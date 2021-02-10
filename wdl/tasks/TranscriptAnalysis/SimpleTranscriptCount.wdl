version 1.0

workflow QuantifyTranscriptReads {

    meta {
        description : "Quantify transcript isoforms simple counting of aligned reads to transcripts.  Assumes each read is a single, complete transcript."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_bam
        File transcript_isoforms_fasta
        String prefix = ""

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        reads_bam : "Bam file containing transcript reads to be aligned to the given transcript isoforms file."
        transcript_isoforms_fasta : "FASTA file containing isoforms sequences to quantify."
        prefix : "[optional] Prefix to prepend to the output file name."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    call QuantifyTranscriptReadsTask {
        input:
            reads_bam                 = reads_bam,
            transcript_isoforms_fasta = transcript_isoforms_fasta,
            prefix                    = prefix
    }

    # ------------------------------------------------
    # Outputs:
    output {
      # Default output file name:
      File aligned_reads  = QuantifyTranscriptReadsTask.aligned_reads
      File unmapped_reads = QuantifyTranscriptReadsTask.unmapped_reads
      File count_table    = QuantifyTranscriptReadsTask.count_table

      File timing_info   = QuantifyTranscriptReadsTask.timing_info
      File memory_log    = QuantifyTranscriptReadsTask.memory_log
    }
}


task QuantifyTranscriptReadsTask {

    meta {
        description : "Quantify transcript isoforms simple counting of aligned reads to transcripts.  Assumes each read is a single, complete transcript."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_bam
        File transcript_isoforms_fasta
        String prefix = ""

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        reads_bam : "Bam file containing transcript reads to be aligned to the given transcript isoforms file."
        transcript_isoforms_fasta : "FASTA file containing isoforms sequences to quantify."
        prefix : "[optional] Prefix to prepend to the output file name."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # ------------------------------------------------
    # Process input args:

    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    String out_base_name="~{prefix}quantified_transcript_reads"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # Triple our input file size for our disk size - should be plenty.
    Float input_files_size_gb = 10*(size(reads_bam, "GiB") + size(transcript_isoforms_fasta, "GiB"))

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 4096
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command <<<
        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -e

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        ##########################################################################################
        # Do the real work here:

        # Create a python file to run containing our info:
        cat >count_sam_transcripts.py <<'END_SCRIPT'
#!/usr/bin/env python
import sys

tx_qual_list_dict = dict()
for line in sys.stdin.readlines():
    sp = line.split("\t")

    tx_name = sp[2]
    qual = int(sp[4])

    try:
        tx_qual_list_dict[tx_name].append(qual)
    except KeyError:
        tx_qual_list_dict[tx_name] = [qual]

print("Transcript_Name\tCount\tMean_Quality")
for tx_name, quals in dict(sorted(tx_qual_list_dict.items(), key=lambda item: len(item[1]))).items():
    print(f"{tx_name}\t{len(quals)}\t{(sum(quals)/len(quals)):2.3f}")

END_SCRIPT
        chmod +x count_sam_transcripts.py

        echo "Indexing input transcript isoforms fasta..."
        samtools faidx ~{transcript_isoforms_fasta}

        echo "Aligning reads and filtering for only primary alignments..."
        samtools fastq ~{reads_bam} \
            | ~/minimap2-2.17_x64-linux/minimap2 -ayYL --MD --eqx -x asm20 -t4 ~{transcript_isoforms_fasta} - \
            | samtools view -hb -F 2304 - \
        > ~{out_base_name}.bam

        echo "Getting unmapped reads..."
        samtools view -hb -f4 ~{out_base_name}.bam > ~{out_base_name}.unmapped.bam

        echo "Quantifying alignments..."
        samtools view -F 2308 ~{out_base_name}.bam | ./count_sam_transcripts.py > ~{out_base_name}transcript_count_table.tsv


        ##########################################################################################

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi

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
         docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: select_first([preemptible_attempts, 0])
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
    output {
      # Default output file name:
      File aligned_reads  = "~{out_base_name}.bam"
      File unmapped_reads = "~{out_base_name}.unmapped.bam"
      File count_table    = "~{out_base_name}transcript_count_table.tsv"

      File timing_info    = "${timing_output_file}"
      File memory_log     = "${memory_log_file}"
    }
}