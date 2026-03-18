# 1. Download Uniref90 database

```
source ~/miniforge3/bin/activate mmseqs2

# 1. Set the destination in your home directory
DB_HOME_DIR="/home/ruchita.solanki/uniref90_new"
mkdir -p "$DB_HOME_DIR"
cd "$DB_HOME_DIR"

# 2. Define a local tmp folder 
# This is where MMseqs2 will store intermediate calculations
mkdir -p "$DB_HOME_DIR/tmp"

# This will download about 60-70GB of data
mmseqs databases UniRef90 uniref90_db "$DB_HOME_DIR/tmp"

# This converts the 735GB RAM requirement into a one-time disk-based index
# This is the step that fixes the "No k-mer could be extracted" error
mmseqs createindex uniref90_db "$DB_HOME_DIR/tmp" --threads 32

echo "Step 3: Verification"
ls -lh

# Optional: Remove the tmp folder to save space after indexing is done
# rm -rf "$DB_HOME_DIR/tmp"
```
# 2. Metaeuk easy predict
```
source ~/miniforge3/bin/activate metaeuk_env

# 1. Paths to your NEW indexed database
DB_PATH="/home/ruchita.solanki/uniref90_new/uniref90_db"
# Use the taxonomy folder created during your successful download
TAX_PATH="/home/ruchita.solanki/uniref90_new/uniref90_db_taxonomy"

# 2. Your project files
WORKDIR="/work/ebg_lab/eb/diatom_consortia/metaeuk"
QUERY_FASTA="$WORKDIR/sr_pypolca_corrected.fasta"
RESULT_DIR="$WORKDIR/easy_predict_results"

mkdir -p "$RESULT_DIR"

# 3. Run Easy-Predict
# This performs: Gene prediction -> Protein Search -> Taxonomic Assignment
echo "Starting MetaEuk easy-predict..."

metaeuk easy-predict \
    "$QUERY_FASTA" \
    "$DB_PATH" \
    "$RESULT_DIR/diatom_results" \
    "$WORKDIR/tmp_folder" \
    --threads 32 \
    --tax-lineage 1 \
    --lca-mode 3 \
    --min-ungapped-score 35 \
    --overlap 1
```
Outputs needed: Polished_MetaEUK.fas, Polished_MetaEUK.headersMap.tsv, Polished_MetaEUK.codon.fas, Polished_MetaEUK.gff. 
# 3. Metaeuk tax to contig
```
source ~/miniforge3/bin/activate metaeuk_env

# 1. Paths - DOUBLE CHECK THESE
DB_PATH="/home/ruchita.solanki/uniref90_new/uniref90_db"
WORKDIR="/work/ebg_lab/eb/diatom_consortia/metaeuk"

# 2. Setup Scratch
SCRATCH_DIR="/scratch/${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"

echo "Copying query files to scratch..."
cp "$WORKDIR/sr_pypolca_corrected.fasta" .
cp "$WORKDIR/Polished_MetaEUK.fas" .
cp "$WORKDIR/Polished_MetaEUK.headersMap.tsv" .

# Step A: Create Query DB from your assembly
echo "Creating query database..."
metaeuk createdb sr_pypolca_corrected.fasta assembly_db

# Step B: Run taxtocontig
# Since uniref90_db.idx now exists in your home folder
# this will skip k-mer counting and start searching immediately!
echo "Running MetaEuk taxtocontig..."
metaeuk taxtocontig \
    assembly_db \
    Polished_MetaEUK.fas \
    Polished_MetaEUK.headersMap.tsv \
    "$DB_PATH" \
    UniRef90_Contig_Taxonomy \
    "$SCRATCH_DIR" \
    --threads 32 \
    --tax-lineage 1 \
    --lca-mode 3

# Step C: Create TSV
echo "Converting results to TSV..."
metaeuk createtsv assembly_db UniRef90_Contig_Taxonomy UniRef90_Contig_Taxonomy.tsv

# 3. Transfer results back
echo "Copying results back to work directory..."
cp UniRef90_Contig_Taxonomy.tsv "$WORKDIR/"

echo "SUCCESS. Results are in $WORKDIR/UniRef90_Contig_Taxonomy.tsv"
```
