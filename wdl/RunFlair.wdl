workflow flair_alignment {
    File fastq
    File genome_fa
    String prefix

    call flair_align {
        input:
            fastq = fastq,
            genome_fa = genome_fa,
            prefix = prefix
    }
}

task flair_align {
    File fastq
    File genome_fa
    String prefix

    command {
        python3 /flair-1.4/flair.py align -r ${fastq} -g ${genome_fa} -t 8 -o ${prefix}
        tar -cvzf ${prefix}.tar.gz ${prefix}
    }

    output {
        File aligned_result = "${prefix}.tar.gz"
        Array[File] all = glob("*")
    }

    runtime {
        docker: "kgarimella/lr-flair:0.01.00"
        cpu: 8
        memory: "60G"
        disks: "local-disk 100 SSD"
        preemptible: 1
    }
}

