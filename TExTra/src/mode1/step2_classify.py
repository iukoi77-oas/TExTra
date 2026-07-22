"""Mode1 TE-exon classification step for prepared metaexon inputs."""

import os
from util.common.write_logs import log_message
from util.qual.run_hitindex import HITindex_pipeline
from util.common.define_layout import CLASSIFY_DIR


def classify_func(buffer_exon_bed, bamfiles, samples, args):
    """Run HITindex-based exon classification from buffered metaexons."""

    log_message("[INFO]", "Step 3-3: Compute HITindex of exons and classify", color="step")
    log_message("[INFO]", "Running the HITindex pipeline.", color="info")
    classify_dir = os.path.join(args.out_dir, CLASSIFY_DIR)
    output_hitindex_dir = os.path.join(classify_dir, "HITindex")
    os.makedirs(output_hitindex_dir, exist_ok=True)
    exon_outnames = HITindex_pipeline(
        buffer_exon_bed,
        bamfiles,
        samples,
        output_hitindex_dir,
        args,
        combined_out_dir=classify_dir,
    )
    log_message("[SUCCESS]", f"The results of the HITindex-based calculations and classifications are stored at {output_hitindex_dir}.", color="success")

    return exon_outnames
