version 1.0

task SubsetCountsMatrixByGenes {

    meta {
        description : "Subsets a count matrix TSV file to contain only the transcripts from the given list of genes.  Assumes the count matrix was produced by comparison with Gencode (due to data formatting) and that the table is a TSV with samples as rows and transcripts as columns."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        Array[String] gene_names
    }

    parameter_meta {
        count_matrix_tsv : "TSV file containing the counts of each transcript expressed in a sample.  (Transcripts in columns.  Samples in rows.  One header row.)"
        gene_names : "Array of gene names for which to keep data from the given count matrix TSV."
    }

    # We're subsetting the file, so we should be able to get away with very little space here.
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size = 3 * ceil(size(count_matrix_tsv, "GB"))

    command {
        /python_scripts/subset_count_matrix_by_gene.py ~{count_matrix_tsv} ~{sep=' ' gene_names}
    }

    output {
        File subset_count_matrix_tsv = "count_matrix_subset_by_gene.tsv"
        File subset_count_matrix_h5ad = "count_matrix_subset_by_gene.h5ad"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.1"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}
