# Test data: paired-end, reverse-stranded (rf), simulated mouse heart development dataset

# Step 1: Read Mapping and Transcriptome Assembly
TExTra prep -i input.tsv -t 8 \
    --index test_result/alignment/index \
    -g GRCm38.primary_assembly.genome.fa \
    -G gencode.vM21.primary_assembly.annotation.gtf \
    -o test_result \
    -r GRCm38_GENCODE_rmsk_TE.gtf \
    --strand rf --readtype paired \
    --de-novo-disable --keep-temp

# Step 2: TE-derived Exon Identification and Classification
TExTra qual -s heart_postnatal_14d,heart_adult_2m,heart_adult_18-20m \
    -t 12 --prep test_result \
    -o test_result \
    --strand rf --readtype paired

# Step 3: Quantification of TE-derived Exons
TExTra quant -o test_result \
    -r GRCm38_GENCODE_rmsk_TE.gtf \
    --strand rf --readtype paired \
    --prep test_result --qual test_result \
    -s heart_postnatal_14d,heart_adult_2m,heart_adult_18-20m \
    -t 12 -e 2

# Step 4: Downstream Analysis of TE-derived Exons (e.g. heart_postnatal_14d vs heart_adult_2m)
TExTra diff -o test_result \
    -g GRCm38.primary_assembly.genome.fa \
    --prep test_result --quant test_result \
    -s heart_postnatal_14d,heart_adult_2m \
    --ncpred --table-only \
    -t 12