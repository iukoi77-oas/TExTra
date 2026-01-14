# TExTra

TExTra is a unified computational framework designed to identify and quantify TE-derived exonization events with locus-level precision from short-read RNA sequencing. The workflow integrates four coordinated modules: (i) consensus reconstruction of TE-chimeric transcripts, (ii) high-confidence identification and positional classification of TE-derived exons, (iii) locus-specific quantification through reallocation of multi-mapping reads, and (iv) downstream functional interpretation. This design enables recovery of complete exon boundaries and accurate quantification of exon usage across samples, overcoming limitations of existing approaches that detect only TE-gene junction fragments.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/iukoi77-oas/TExTra.git
cd TExTra
```

### 2. Create and activate the conda environment

```bash
conda env create -n TExTra -f environment.yml
conda activate TExTra
```

Then install TExTra in editable mode:

```bash
pip install -e .
```

### 3. Common installation issues and solutions

#### (1) Conda environment creation fails due to strict channel priority

If conda fails to resolve dependencies because `channel_priority` is set to `strict`, run:

```bash
conda config --set channel_priority flexible
```

Then retry environment creation.

#### (2) PLEK2 dependency issues

PLEK2 requires several Python packages that may not be correctly resolved by conda.
If you encounter import or runtime errors related to PLEK2, install the following packages manually:

```bash
pip install keras==2.4.3
pip install tensorflow==2.4.1
pip install regex
pip install bio
```

#### (3) External tool paths (PLEK2 and TACO)

TExTra expects external tools to be placed under the `util/` directory:

```text
TExTra/
├── util/
│   ├── taco*/          # TACO binary directory
│   └── PLEK*/            # PLEK2 directory
```

> **Note**: The TACO and PLEK2 tools are already provided as compressed packages under the util/ directory.
They can be used directly by extracting the corresponding archives in place, and no additional download or reinstallation is required.

* Ensure that **TACO** and **PLEK2** are installed and executable.
* The `taco_run` executable should be located under `util/taco*/`.
* PLEK2.py should be located under `util/PLEK*/`.

If installed elsewhere, either move them into the corresponding directories or update the paths in the TExTra configuration accordingly.

### 4. Verify installation

After installation, you should be able to run:

```bash
TExTra --help
```

and see the available subcommands:

```text
prep     Read mapping and transcriptome assembly
qual     TE-derived exon identification and classification
quant    Quantification of TE-derived exons
diff     Downstream analysis
```
