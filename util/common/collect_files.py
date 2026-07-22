"""Collect common file paths, sample names, and replicate BAM inputs."""

import os


def parse_sample_csv(sample_csv):
    """Parse comma-separated sample names and drop empty tokens."""
    return [s.strip() for s in str(sample_csv).split(",") if s.strip()]


def _replicate_folder_belongs_to_sample(replicate_folder, sample):
    """Return True when a prep replicate folder exactly belongs to sample."""
    prefix = f"{sample}_rep"
    if not replicate_folder.startswith(prefix):
        return False
    return replicate_folder[len(prefix):].isdigit()


def collect_bam_and_replicates(alignment_dir, sample_list):
    """
    Collect accepted-hit BAM paths and replicate folders for each sample.

    The returned dictionaries are keyed by sample name:
      - bamfiles_dict[sample] -> list[str]
      - replicates_dict[sample] -> list[str]
    """
    bamfiles_dict = {sample: [] for sample in sample_list}
    replicates_dict = {sample: [] for sample in sample_list}
    if not os.path.isdir(alignment_dir):
        return bamfiles_dict, replicates_dict

    for replicate_folder in sorted(os.listdir(alignment_dir)):
        folder_path = os.path.join(alignment_dir, replicate_folder)
        if not os.path.isdir(folder_path):
            continue
        for sample in sample_list:
            if not _replicate_folder_belongs_to_sample(replicate_folder, sample):
                continue
            bam_path = os.path.join(folder_path, f"{replicate_folder}_accepted_hits.bam")
            if os.path.isfile(bam_path):
                bamfiles_dict[sample].append(bam_path)
                replicates_dict[sample].append(replicate_folder)
    return bamfiles_dict, replicates_dict
