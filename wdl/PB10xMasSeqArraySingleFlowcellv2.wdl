version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF
import "tasks/AlignReads.wdl" as AR
import "tasks/Cartographer.wdl" as CART
import "tasks/TranscriptAnalysis/Flair_Tasks.wdl" as ISO
import "tasks/ReadsMetrics.wdl" as RM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/JupyterNotebooks.wdl" as JUPYTER
import "tasks/Annmas.wdl" as ANNMAS

workflow PB10xMasSeqSingleFlowcellv2 {

    meta {
        description : "This workflow is designed to process data from the MASSeq v2 protocol and produce aligned reads that are ready for downstream analysis (e.g. transcript isoform identification).  It takes in a raw PacBio run folder location on GCS and produces a folder containing the aligned reads and other processed data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/PB10xMasSeqSingleFlowcellv2"

        File segments_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/cDNA_array_15x.unique_seqs_for_cartographer.fasta"
        File boundaries_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/bounds_file_for_extraction.txt"

        File head_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/10x_adapter.fasta"
        File tail_adapter_fasta = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/tso_adapter.fasta"

        File ref_fasta =  "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.fasta.fai"
        File ref_fasta_dict = "gs://broad-dsde-methods-long-reads/resources/references/grch38/Homo_sapiens_assembly38.dict"

        File jupyter_template_static = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/MAS-seq_QC_report_template-static.ipynb"
        File workflow_dot_file = "gs://broad-dsde-methods-long-reads/resources/MASseq_0.0.2/PB10xMasSeqArraySingleFlowcellv2.dot"

        Int min_ccs_passes = 2

        String? sample_name
    }

    parameter_meta {
        gcs_input_dir : "Input folder on GCS in which to search for BAM files to process."
        gcs_out_root_dir : "Root output GCS folder in which to place results of this workflow."

        segments_fasta : "FASTA file containing unique segments for which to search in the given BAM files.   These segments are used as delimiters in the reads.  Read splitting uses these delimiters and the boundaries file."
        boundaries_file : "Text file containing two comma-separated segment names from the segments_fasta on each line.  These entries define delimited sections to be extracted from the reads and treated as individual array elements."

        head_adapter_fasta : "FASTA file containing the sequence that each transcript should start with.  Typically this will be the 10x adapter sequence from the 10x library prep."
        tail_adapter_fasta : "FASTA file containing the sequence that each transcript should end with.  Typically this will be the Template Switch Oligo (TSO) sequence from the 10x library prep."

        ref_fasta : "FASTA file containing the reference sequence to which the input data should be aligned."
        ref_fasta_index : "FASTA index file for the given ref_fasta file."
        ref_fasta_dict : "Sequence dictionary file for the given ref_fasta file."

        jupyter_template_static : "Jupyter notebook / ipynb file containing a template for the QC report which will contain static plots.  This should contain the same information as the jupyter_template_interactive file, but with static images."
        workflow_dot_file : "DOT file containing the representation of this WDL to be included in the QC reports.  This can be generated with womtool."

        min_ccs_passes : "[optional] Minimum number of passes required for CCS to take place on a ZMW (Default: 2)."

        sample_name : "[optional] The name of the sample to associate with the data in this workflow."
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {

        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String RG_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String RG_array_elements = "@RG\\tID:~{ID}.array_elements\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # Get statistics on polymerase reads:
        call PB.CollectPolymeraseReadLengths {
            input:
                gcs_input_dir = gcs_input_dir,
                subreads = subread_bam,
                prefix = SM + "_polymerase_read_lengths"
        }

        call Utils.ShardLongReadsWithCopy { input: unmapped_files = [ subread_bam ], num_reads_per_split = 100000 }

        scatter (subreads in ShardLongReadsWithCopy.unmapped_shards) {

            ## No more preemption on this sharding - takes too long otherwise.
            RuntimeAttr disable_preemption_runtime_attrs = object {
                preemptible_tries: 0
            }

            # Call CCS on the subreads from the sequencer:
            # No preepting because these take long enough that it doesn't seem to save $
            # 16 gigs of memory because I had a crash at 8
            RuntimeAttr ccs_runtime_attrs = object {
                mem_gb: 16,
                preemptible_tries: 0
            }
            call PB.CCS {
                input:
                    subreads = subreads,
                    min_passes = min_ccs_passes,
                    disk_space_scale_factor = 4,
                    runtime_attr_override = ccs_runtime_attrs
            }

            # Get our uncorrected / CCS Rejected reads:
            call PB.ExtractUncorrectedReads {
                input:
                    subreads = subreads,
                    consensus = CCS.consensus,
                    prefix = SM + "_ccs_rejected",
                    runtime_attr_override = disable_preemption_runtime_attrs
            }

            call Utils.ShardLongReadsWithCopy as SubshardRawSubreads {
                input:
                    unmapped_files = [ subreads ],
                    num_reads_per_split = 50000,
                    runtime_attr_override = disable_preemption_runtime_attrs
            }
            scatter (subsharded_subreads in SubshardRawSubreads.unmapped_shards) {

                # Get ZMW Subread stats here to shard them out wider and make it faster:
                call PB.CollectZmwSubreadStats as CollectZmwSubreadStats_subsharded {
                    input:
                        subreads = subsharded_subreads,
                        prefix = SM + "_zmw_subread_stats"
                }

                # Get approximate subread array lengths here:
                call CART.GetApproxRawSubreadArrayLengths as GetApproxRawSubreadArrayLengths_subsharded {
                    input:
                        reads_file = subsharded_subreads,
                        delimiters_fasta = segments_fasta,
                        min_qual = 7.0,
                        ignore_seqs = ["Poly_A", "Poly_T", "3_prime_TSO", "5_prime_TSO"],
                        prefix = SM + "_approx_raw_subread_array_lengths"
                }
            }

            # Merge our micro-shards of subread stats:
            call Utils.MergeTsvFiles as MergeMicroShardedZmwSubreadStats {
                input:
                    tsv_files = CollectZmwSubreadStats_subsharded.zmw_subread_stats
            }

            # Merge the micro-sharded raw subread array element counts:
            call Utils.MergeTsvFiles as MergeMicroShardedRawSubreadArrayElementCounts {
                input:
                    tsv_files = GetApproxRawSubreadArrayLengths_subsharded.approx_subread_array_lengths
            }

            # Shard these reads even wider so we can make sure we don't run out of memory:
            call Utils.ShardLongReadsWithCopy as ShardCorrectedReads { input: unmapped_files = [ CCS.consensus ], num_reads_per_split = 20000 }

            scatter (corrected_shard in ShardCorrectedReads.unmapped_shards) {

                call ANNMAS.Annotate as AnnotateReads {
                    input:
                        reads = corrected_shard
                }

                call ANNMAS.Segment as SegmentAnnotatedReads {
                    input:
                        annotated_reads = AnnotateReads.annotated_bam
                }
            }

            # Merge all outputs of Annmas Annotate / Segment:
            call Utils.MergeBams as MergeAnnotatedReads_1 {
                input:
                    bams = AnnotateReads.annotated_bam,
                    prefix = SM + "_AnnotatedReads_intermediate_1"
            }
            call Utils.MergeBams as MergeArrayElements_1 {
                input:
                    bams = SegmentAnnotatedReads.segmented_bam,
                    prefix = SM + "_ArrayElements_intermediate_1"
            }

            ## For some reason we need a LOT of memory for this.
            ## Need to debug it or remove the alignment of CCS (non-split) reads:
            RuntimeAttr align_ccs_reads_runtime_attrs = object {
                mem_gb: 60,
                preemptible_tries: 0
            }
            call AR.Minimap2 as AlignCCSReads {
                input:
                    reads      = [ CCS.consensus ],
                    ref_fasta  = ref_fasta,
                    RG         = RG_consensus,
                    map_preset = "splice",
                    runtime_attr_override = align_ccs_reads_runtime_attrs
            }

            call TENX.AnnotateBarcodesAndUMIs as TenxAnnotateArrayElements {
                input:
                    bam_file = MergeAnnotatedReads_1.merged_bam,
                    bam_index = MergeAnnotatedReads_1.merged_bai,
                    head_adapter_fasta = head_adapter_fasta,
                    tail_adapter_fasta = tail_adapter_fasta,
                    read_end_length = 200,
                    poly_t_length = 31,
                    barcode_length = 16,
                    umi_length = 10,
                    mem_gb = 32
            }

            call AR.Minimap2 as AlignArrayElements {
                input:
                    reads      = [ TenxAnnotateArrayElements.output_bam ],
                    ref_fasta  = ref_fasta,
                    RG         = RG_array_elements,
                    map_preset = "splice"
            }
        }

        # Merge all sharded merged outputs from annotating / splitting:
        call Utils.MergeBams as MergeAnnotatedReads_2 {
            input:
                bams = MergeAnnotatedReads_1.merged_bam,
                prefix = SM + "_AnnotatedReads_intermediate_2"
        }
        call Utils.MergeBams as MergeArrayElements_2 {
            input:
                bams = MergeArrayElements_1.merged_bam,
                prefix = SM + "_ArrayElements_intermediate_2"
        }

        # Merge the sharded zmw subread stats:
        call Utils.MergeTsvFiles as MergeShardedZmwSubreadStats {
            input:
                tsv_files = MergeMicroShardedZmwSubreadStats.merged_tsv,
                prefix = SM + "_zmw_subread_stats"
        }

        # Merge the sharded raw subread array element counts:
        call Utils.MergeTsvFiles as MergeShardedRawSubreadArrayElementCounts {
            input:
                tsv_files = MergeMicroShardedRawSubreadArrayElementCounts.merged_tsv,
                prefix = SM + "_approx_subread_array_lengths"
        }

        # Merge the 10x stats:
        call Utils.MergeCountTsvFiles as Merge10XStats_1 { input: count_tsv_files = TenxAnnotateArrayElements.stats }

        # Merge all Aligned array elements together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedAnnoatatedArrayElementChunk { input: bams = TenxAnnotateArrayElements.output_bam }

        # Merge all Aligned array elements together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedArrayElementChunk { input: bams = AlignArrayElements.aligned_bam }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus }

        # Merge all CCS Rejected bams together:
        call Utils.MergeBams as MergeCCSRejectedChunks { input: bams = ExtractUncorrectedReads.uncorrected }

        # Merge all CCS bams together for this Subread BAM:
        call Utils.MergeBams as MergeAlignedChunks { input: bams = AlignCCSReads.aligned_bam }

        # Merge all CCS reports together for this Subread BAM:
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }

        # Collect metrics on the subreads bam:
        RuntimeAttr subreads_sam_stats_runtime_attrs = object {
            cpu_cores:          2,
            mem_gb:             8,
            disk_gb:            ceil(3 * size(subread_bam, "GiB")),
            boot_disk_gb:       10,
            preemptible_tries:  0,
            max_retries:        1,
            docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.8"
        }
        call AM.SamtoolsStats as CalcSamStatsOnInputBam {
            input:
                bam = subread_bam,
                runtime_attr_override = subreads_sam_stats_runtime_attrs
        }
        call FF.FinalizeToDir as FinalizeSamStatsOnInputBam {
            input:
                # an unfortunate hard-coded path here:
                outdir = outdir + "/" + DIR + "/" + WdlExecutionStartTimestamp.timestamp_string + "/metrics/input_bam_stats",
                files = [
                    CalcSamStatsOnInputBam.raw_stats,
                    CalcSamStatsOnInputBam.summary_stats,
                    CalcSamStatsOnInputBam.first_frag_qual,
                    CalcSamStatsOnInputBam.last_frag_qual,
                    CalcSamStatsOnInputBam.first_frag_gc_content,
                    CalcSamStatsOnInputBam.last_frag_gc_content,
                    CalcSamStatsOnInputBam.acgt_content_per_cycle,
                    CalcSamStatsOnInputBam.insert_size,
                    CalcSamStatsOnInputBam.read_length_dist,
                    CalcSamStatsOnInputBam.indel_distribution,
                    CalcSamStatsOnInputBam.indels_per_cycle,
                    CalcSamStatsOnInputBam.coverage_distribution,
                    CalcSamStatsOnInputBam.gc_depth
                ]
        }
    }

    # Merge all sharded merged sharded outputs for Annotation / Segmentation:
    # TODO: Fix SM[0] to be the right sample name if multiple samples are found.
    call Utils.MergeBams as MergeAnnotatedReads_3 {
        input:
            bams = MergeAnnotatedReads_2.merged_bam,
            prefix = SM[0] + "_AnnotatedReads"
    }
    call Utils.MergeBams as MergeArrayElements_3 {
        input:
            bams = MergeArrayElements_2.merged_bam,
            prefix = SM[0] + "_ArrayElements"
    }
    call PB.PBIndex as PbIndexAnnotatedReads {
        input:
            bam = MergeAnnotatedReads_3.merged_bam
    }
    call PB.PBIndex as PbIndexArrayElements {
        input:
            bam = MergeArrayElements_3.merged_bam
    }

    # Merge the 10x merged stats:
    call Utils.MergeCountTsvFiles as Merge10XStats_2 { input: count_tsv_files = Merge10XStats_1.merged_tsv }

    # Merge all array element bams together for this flowcell:
    call Utils.MergeBams as MergeAnnotatedAlignedArrayElements { input: bams = MergeAlignedAnnoatatedArrayElementChunk.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all array element bams together for this flowcell:
    call Utils.MergeBams as MergeAlignedArrayElements { input: bams = MergeAlignedArrayElementChunk.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all CCS Rejected chunks:
    call Utils.MergeBams as MergeAllCCSRejectedBams { input: bams = MergeCCSRejectedChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}.ccs_rejected" }

    # Merge all aligned CCS bams together for this flowcell:
    call Utils.MergeBams as MergeAllAlignedCCSBams { input: bams = MergeAlignedChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all CCS bams together for this flowcell:
    call Utils.MergeBams as MergeAllCCSBams { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    # Merge all CCS reports together for this flowcell:
    call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

    ##########
    # Metrics and plots
    ##########

    String base_out_dir = outdir + "/" + DIR[0] + "/" + WdlExecutionStartTimestamp.timestamp_string
    String metrics_out_dir = base_out_dir + "/metrics"

    # Aligned CCS Metrics:
    call RM.CalculateAndFinalizeReadMetrics as AlignedCCSMetrics {
        input:
            bam_file = MergeAllAlignedCCSBams.merged_bam,
            bam_index = MergeAllAlignedCCSBams.merged_bai,
            ref_dict = ref_fasta_dict,

            base_metrics_out_dir = metrics_out_dir + "/aligned_ccs_metrics"
    }

    # Aligned Array Element Metrics:
    call RM.CalculateAndFinalizeAlternateReadMetrics as AlignedArrayElementMetrics {
        input:
            bam_file = MergeAlignedArrayElements.merged_bam,
            bam_index = MergeAlignedArrayElements.merged_bai,
            ref_dict = ref_fasta_dict,

            base_metrics_out_dir = metrics_out_dir + "/aligned_array_element_metrics"
    }

    ##########
    # Create Report:
    ##########

    ## NOTE: This assumes ONE file for both the raw input and the 10x array element stats!
    ##       This should be fixed in version 2.
    RuntimeAttr create_report_runtime_attrs = object {
            preemptible_tries:  0
    }
    ## TODO: THIS WILL CURRENTLY FAIL.  THIS IS EXPECTED BECAUSE I HAVE NOT FIXED IT YET.
    call JUPYTER.PB10xMasSeqSingleFlowcellReport as GenerateStaticReport {
        input:
            notebook_template                = jupyter_template_static,

            subreads_stats                   = CalcSamStatsOnInputBam.raw_stats[0],
            ccs_reads_stats                  = AlignedCCSMetrics.sam_stats_raw_stats,
            array_elements_stats             = AlignedArrayElementMetrics.sam_stats_raw_stats,
            ccs_report_file                  = MergeAllCCSReports.report,

            ccs_bam_file                     = MergeAllCCSBams.merged_bam,
            array_element_bam_file           = MergeAlignedArrayElements.merged_bam,
            ccs_rejected_bam_file            = MergeAllCCSRejectedBams.merged_bam,

            annotated_bam_file               = MergeAnnotatedReads_3.merged_bam,

            zmw_subread_stats_file           = MergeShardedZmwSubreadStats.merged_tsv[0],
            polymerase_read_lengths_file     = CollectPolymeraseReadLengths.polymerase_read_lengths_tsv[0],
            approx_raw_subread_array_lengths = MergeShardedRawSubreadArrayElementCounts.merged_tsv[0],

            ten_x_metrics_file               = Merge10XStats_2.merged_tsv,

            workflow_dot_file                = workflow_dot_file,
            prefix                           = SM[0] + "_masseq_",

            runtime_attr_override            = create_report_runtime_attrs,
    }

    ##########
    # Finalize
    ##########

    # NOTE: We key all finalization steps on the static report.
    #       This will prevent incomplete runs from being placed in the output folders.

    # TODO: Should make this iterate through all found bam files, not just the first one.  This seems to be the case for all workflows right now though...

    # Finalize the notebook:
    String static_report_dir = metrics_out_dir + "/report_static"
    call FF.FinalizeToDir as FinalizeStaticReport {
        input:
            files = [
                GenerateStaticReport.populated_notebook,
                GenerateStaticReport.html_report,
                GenerateStaticReport.pdf_report
            ],
            outdir = static_report_dir,
            keyfile = GenerateStaticReport.html_report
    }

    # Finalize all the Extracted Bounded Regions data:
    String annotatedReadsDir = base_out_dir + "/annmas"
    call FF.FinalizeToDir as FinalizeAnnotatedReads {
        input:
            files = [
                MergeAnnotatedReads_3.merged_bam,
                MergeAnnotatedReads_3.merged_bai,
                PbIndexAnnotatedReads.pbindex,

                MergeArrayElements_3.merged_bam,
                MergeArrayElements_3.merged_bai,
                PbIndexArrayElements.pbindex
            ],
            outdir = annotatedReadsDir,
            keyfile = GenerateStaticReport.html_report
    }

    # Finalize all the 10x metrics here:
    String tenXToolMetricsDir = metrics_out_dir + "/ten_x_tool_metrics"
    scatter ( i in range(length(TenxAnnotateArrayElements.output_bam[0]))) {
        call FF.FinalizeToDir as FinalizeTenXRgStats {
            input:
                files = [
                    TenxAnnotateArrayElements.barcode_stats[0][i],
                    TenxAnnotateArrayElements.starcode[0][i],
                    TenxAnnotateArrayElements.stats[0][i],
                    TenxAnnotateArrayElements.timing_info[0][i]
                ],
                outdir = tenXToolMetricsDir + "/" + i,
                keyfile = GenerateStaticReport.html_report
        }
    }

    call FF.FinalizeToDir as FinalizeAnnotatedArrayElements {
        input:
            files = [ MergeAnnotatedAlignedArrayElements.merged_bam, MergeAnnotatedAlignedArrayElements.merged_bai ],
            outdir = base_out_dir + "/annotated_aligned_array_elements",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeArrayElements {
        input:
            files = [ MergeAlignedArrayElements.merged_bam, MergeAlignedArrayElements.merged_bai ],
            outdir = base_out_dir + "/aligned_array_elements",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeAlignedCCSBams {
        input:
            files = [ MergeAllAlignedCCSBams.merged_bam, MergeAllAlignedCCSBams.merged_bai ],
            outdir = base_out_dir + "/merged_bams/aligned",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeAlignedCCSRejectedBams {
        input:
            files = [ MergeAllCCSRejectedBams.merged_bam, MergeAllCCSRejectedBams.merged_bai ],
            outdir = base_out_dir + "/merged_bams/ccs_rejected",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeCCSBams {
        input:
            files = [ MergeAllCCSBams.merged_bam, MergeAllCCSBams.merged_bai ],
            outdir = base_out_dir + "/merged_bams/unaligned",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ MergeAllCCSReports.report ],
            outdir = metrics_out_dir + "/ccs_metrics",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeZmwSubreadStats {
        input:
            files = MergeShardedZmwSubreadStats.merged_tsv,
            outdir = metrics_out_dir + "/ccs_metrics",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizeRawSubreadArrayElementCounts {
        input:
            files = MergeShardedRawSubreadArrayElementCounts.merged_tsv,
            outdir = metrics_out_dir + "/array_stats",
            keyfile = GenerateStaticReport.html_report
    }

    call FF.FinalizeToDir as FinalizePolymeraseReadLengths {
        input:
            files = CollectPolymeraseReadLengths.polymerase_read_lengths_tsv,
            outdir = metrics_out_dir + "/array_stats",
            keyfile = GenerateStaticReport.html_report
    }

    # Write out completion file so in the future we can be 100% sure that this run was good:
    call FF.WriteCompletionFile {
        input:
            outdir = base_out_dir + "/",
            keyfile = GenerateStaticReport.html_report
    }
}
