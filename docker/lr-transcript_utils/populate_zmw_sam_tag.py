import pysam
import time

ZMW_TAG = 'zm'

if __name__ == '__main__':

    # Right now we accept only 1 argument - the directory in which the quant files live:
    if len(sys.argv) != 2:
        print(f"{sys.argv[0]} PAC_BIO_BAM", file=sys.stderr)
        print(f"Error: you must specify at least one gene name.", file=sys.stderr)
        sys.exit(1)

    reads_in_intermediate_file = 0
    with pysam.AlignmentFile(analysis_name + '.intermediate.bam', 'rb', check_sq=False) as bam_file:
        with pysam.AlignmentFile(analysis_name + '.bam', 'wb', check_sq=False, header=bam_file.header) as output_file:
            for read in bam_file.fetch(until_eof=True):
                if reads_in_intermediate_file % 1000 == 0:
                    current_time = time.time()
                    print('Processed reads: {:,}. Time elapsed: {:.2f}s'.format(reads_in_intermediate_file, current_time - last_timing))
                    last_timing = current_time

                read.set_tag(BARCODE_TAG, read.get_tag(RAW_BARCODE_TAG), value_type='Z')

                output_file.write(read)

                reads_in_intermediate_file += 1