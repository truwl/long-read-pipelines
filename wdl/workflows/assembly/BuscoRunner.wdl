version 1.0

import "../../tasks/metrics/Busco.wdl" as Busco
import "../../tasks/utils/Finalize.wdl" as FF

workflow BuscoRunner {
    input {
        File assembly
        String lineage

        String out_dir
    }

    call Busco.Busco {
        input:
            assembly = assembly,
            lineage = lineage
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [Busco.output_tar],
            outdir = out_dir
    }
}
