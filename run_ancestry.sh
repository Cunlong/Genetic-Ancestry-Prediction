#!/usr/bin/env bash
#
# run_ancestry.sh -- genetic ancestry inference with PLINK 1.9

set -euo pipefail

# settings 
TARGET=${TARGET:-ATP_data}        # your bfile prefix, build 37, no extension
REFDIR=${REFDIR:-./1000G}
OUTDIR=${OUTDIR:-./ancestry_out}
PLINK=${PLINK:-plink}             # PLINK 1.9 binary
NPC=${NPC:-10}                    # PCs for ancestry assignment
NPC_GWAS=${NPC_GWAS:-20}          # PCs computed within your cohort, for GWAS
MEMORY=${MEMORY:-8000}            # MB
THREADS=${THREADS:-4}
MIN_SHARED=${MIN_SHARED:-20000}
HWE_P=${HWE_P:-1e-6}        # tested WITHIN each reference population, never pooled

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP="$OUTDIR/tmp"
mkdir -p "$TMP"
PK="$PLINK --memory $MEMORY --threads $THREADS --allow-no-sex"

link_or_copy () { ln -f "$1" "$2" 2>/dev/null || cp -f "$1" "$2"; }

echo "=============================================================="
echo " cohort    : $TARGET"
echo " reference : $REFDIR/1000G_ref"
echo " output    : $OUTDIR"
echo "=============================================================="

command -v "$PLINK" >/dev/null || {
  echo "ERROR: '$PLINK' not on PATH. See README.md section 1." >&2; exit 1; }
for e in bed bim fam; do
  [[ -f "$TARGET.$e" ]] || { echo "ERROR: $TARGET.$e not found." >&2; exit 1; }
done
[[ -f "$REFDIR/1000G_ref.bed" ]] || {
  echo "ERROR: reference not found. Run 00_setup_reference.sh first." >&2; exit 1; }
[[ -f "$OUTDIR/ref_labels.tsv" ]] || {
  echo "ERROR: $OUTDIR/ref_labels.tsv not found. Run 00_setup_reference.sh first." >&2; exit 1; }

# 1. CLEAN 

echo
echo "[1/7] cleaning"
$PK --bfile "$TARGET" --autosome --snps-only just-acgt --biallelic-only strict \
    --make-bed --out "$TMP/target_clean" > "$TMP/log.1.target.txt"
if [[ -f "$REFDIR/1000G_ref.ready" ]]; then
  # 00_setup_reference.sh already applied exactly these filters; re-running them
  # would only write out a second multi-GB copy of the same data.
  echo "      reference already filtered by 00_setup_reference.sh, linking"
  for e in bed bim fam; do link_or_copy "$REFDIR/1000G_ref.$e" "$TMP/ref_clean.$e"; done
else
  $PK --bfile "$REFDIR/1000G_ref" --autosome --snps-only just-acgt --biallelic-only strict \
      --make-bed --out "$TMP/ref_clean" > "$TMP/log.1.ref.txt"
fi
echo "      cohort:    $(wc -l < "$TMP/target_clean.fam") samples, $(wc -l < "$TMP/target_clean.bim") variants"
echo "      reference: $(wc -l < "$TMP/ref_clean.fam") samples, $(wc -l < "$TMP/ref_clean.bim") variants"

# 2. SAMPLE IDS 

echo
echo "[2/7] checking sample IDs"
awk '{print $2}' "$TMP/target_clean.fam" | sort > "$TMP/ids_target.txt"
awk '{print $2}' "$TMP/ref_clean.fam"    | sort > "$TMP/ids_ref.txt"
comm -12 "$TMP/ids_target.txt" "$TMP/ids_ref.txt" > "$TMP/ids_clash.txt" || true

if [[ -s "$TMP/ids_clash.txt" ]]; then
  N=$(wc -l < "$TMP/ids_clash.txt")
  echo "      $N ID(s) clash with the reference -- prefixing the cohort with 'COHORT_'"
  awk '{print $1, $2, $1, "COHORT_"$2}' "$TMP/target_clean.fam" > "$TMP/rename.txt"
  $PK --bfile "$TMP/target_clean" --update-ids "$TMP/rename.txt" \
      --make-bed --out "$TMP/target_ids" > "$TMP/log.2.rename.txt"
else
  echo "      no clashes"
  for e in bed bim fam; do link_or_copy "$TMP/target_clean.$e" "$TMP/target_ids.$e"; done
fi

# 3. HARMONISE 

echo
echo "[3/7] harmonising alleles"
python3 "$HERE/harmonise.py" \
    --ref-bim        "$TMP/ref_clean.bim" \
    --target-bim     "$TMP/target_ids.bim" \
    --out-ref-bim    "$TMP/ref_h.bim" \
    --out-target-bim "$TMP/target_h.bim" \
    --out-shared     "$TMP/shared.snps" \
    --min-shared     "$MIN_SHARED"
# only the .bim is rewritten by harmonise.py, so the .bed/.fam are linked
link_or_copy "$TMP/ref_clean.bed" "$TMP/ref_h.bed"
link_or_copy "$TMP/ref_clean.fam" "$TMP/ref_h.fam"
link_or_copy "$TMP/target_ids.bed" "$TMP/target_h.bed"
link_or_copy "$TMP/target_ids.fam" "$TMP/target_h.fam"

# 4. LD PRUNE 

echo
echo "[4/7] LD pruning (reference only)"
cat > "$TMP/lrld_b37.txt" <<'EOF'
1 48000000 52000000 lrld1
2 86000000 100500000 lrld2
2 134500000 138000000 lrld3
2 183000000 190000000 lrld4
3 47500000 50000000 lrld5
3 83500000 87000000 lrld6
3 89000000 97500000 lrld7
5 44500000 50500000 lrld8
5 98000000 100500000 lrld9
5 129000000 132000000 lrld10
5 135500000 138500000 lrld11
6 25000000 35000000 MHC
6 57000000 64000000 lrld12
6 140000000 142500000 lrld13
7 55000000 66000000 lrld14
8 7000000 13000000 inv8p23
8 43000000 50000000 lrld15
8 112000000 115000000 lrld16
10 37000000 43000000 lrld17
11 46000000 57000000 lrld18
11 87500000 90500000 lrld19
12 33000000 40000000 lrld20
12 109500000 112000000 lrld21
17 40900000 44900000 inv17q21
20 32000000 34500000 lrld22
EOF

$PK --bfile "$TMP/ref_h" --extract "$TMP/shared.snps" \
    --make-bed --out "$TMP/ref_shared" > "$TMP/log.4.subset.txt"

POPS=$(awk 'NR>1{print $2}' "$OUTDIR/ref_labels.tsv" | sort -u)
N_POPS=0
: > "$TMP/hwe_all.txt"
for POP in $POPS; do
  awk -v p="$POP" 'NR==FNR { if (FNR>1 && $2==p) k[$1]; next } ($2 in k){print $1, $2}' \
      "$OUTDIR/ref_labels.tsv" "$TMP/ref_shared.fam" > "$TMP/keep.$POP"
  N_IN=$(wc -l < "$TMP/keep.$POP")
  if (( N_IN < 20 )); then
    echo "      skipping $POP for HWE (only $N_IN samples)"
    continue
  fi
  $PK --bfile "$TMP/ref_shared" --keep "$TMP/keep.$POP" --hwe "$HWE_P" \
      --write-snplist --out "$TMP/hwe.$POP" > /dev/null
  cat "$TMP/hwe.$POP.snplist" >> "$TMP/hwe_all.txt"
  N_POPS=$((N_POPS + 1))
done
# a SNP is kept only if it passes in every population tested
sort "$TMP/hwe_all.txt" | uniq -c | awk -v n="$N_POPS" '$1==n{print $2}' > "$TMP/hwe_pass.snps"
echo "      HWE within $N_POPS populations: $(wc -l < "$TMP/hwe_pass.snps") of $(wc -l < "$TMP/shared.snps") SNPs pass everywhere"

$PK --bfile "$TMP/ref_h" --extract "$TMP/hwe_pass.snps" \
    --exclude range "$TMP/lrld_b37.txt" \
    --maf 0.05 --geno 0.02 \
    --indep-pairwise 1000 50 0.1 \
    --out "$TMP/prune" > "$TMP/log.4.prune.txt"

N_PRUNE=$(wc -l < "$TMP/prune.prune.in")
echo "      SNPs kept: $N_PRUNE"
(( N_PRUNE < 10000 )) && echo "      WARNING: under 10,000 SNPs -- PCs will be noisy." >&2

$PK --bfile "$TMP/ref_h"    --extract "$TMP/prune.prune.in" \
    --make-bed --out "$TMP/ref_sub"    > "$TMP/log.4.sub_ref.txt"
$PK --bfile "$TMP/target_h" --extract "$TMP/prune.prune.in" \
    --make-bed --out "$TMP/target_sub" > "$TMP/log.4.sub_target.txt"

# 5. SAMPLE QC 
echo
echo "[5/7] sample QC on the cohort"

$PK --bfile "$TMP/target_h" --missing --out "$OUTDIR/qc" > "$TMP/log.5.missing.txt"

$PK --bfile "$TMP/target_sub" --het --out "$OUTDIR/qc" > "$TMP/log.5.het.txt"
echo "      wrote $OUTDIR/qc.imiss (all SNPs) and $OUTDIR/qc.het (pruned)"

#  6. MERGE AND PCA echo
echo "[6/7] merging and running PCA"
if ! $PK --bfile "$TMP/ref_sub" --bmerge "$TMP/target_sub" \
         --make-bed --out "$TMP/merged" > "$TMP/log.6.merge.txt" 2>&1; then
  if [[ -f "$TMP/merged-merge.missnp" ]]; then
    echo "      $(wc -l < "$TMP/merged-merge.missnp") conflicting SNPs dropped, retrying"
    $PK --bfile "$TMP/ref_sub"    --exclude "$TMP/merged-merge.missnp" \
        --make-bed --out "$TMP/ref_sub2"    > /dev/null
    $PK --bfile "$TMP/target_sub" --exclude "$TMP/merged-merge.missnp" \
        --make-bed --out "$TMP/target_sub2" > /dev/null
    $PK --bfile "$TMP/ref_sub2" --bmerge "$TMP/target_sub2" \
        --make-bed --out "$TMP/merged" > "$TMP/log.6.merge2.txt"
  else
    echo "ERROR: merge failed -- see $TMP/log.6.merge.txt" >&2
    exit 1
  fi
fi
echo "      merged: $(wc -l < "$TMP/merged.fam") samples, $(wc -l < "$TMP/merged.bim") SNPs"

$PK --bfile "$TMP/merged" --maf 0.01 --pca "$NPC" header tabs \
    --out "$OUTDIR/merged_pca" > "$TMP/log.6.pca.txt"

$PK --bfile "$TMP/target_h" \
    --extract "$TMP/shared.snps" --exclude range "$TMP/lrld_b37.txt" \
    --maf 0.05 --geno 0.02 --indep-pairwise 1000 50 0.1 \
    --out "$TMP/prune_cohort" > "$TMP/log.6.prune_cohort.txt"
N_COH=$(wc -l < "$TMP/prune_cohort.prune.in")
echo "      cohort-specific pruning kept $N_COH SNPs"

$PK --bfile "$TMP/target_h" --extract "$TMP/prune_cohort.prune.in" \
    --pca "$NPC_GWAS" header tabs \
    --out "$OUTDIR/cohort_pca" > "$TMP/log.6.pca_cohort.txt"
echo "      wrote merged_pca.eigenvec (ancestry) and cohort_pca.eigenvec (GWAS covariates)"

# 7. CLASSIFY 
echo
echo "[7/7] classifying"
python3 "$HERE/classify_ancestry.py" \
    --outdir   "$OUTDIR" \
    --eigenvec "$OUTDIR/merged_pca.eigenvec" \
    --qc       "$OUTDIR/qc"

cat <<EOF

==============================================================
 Done.

   $OUTDIR/target_ancestry.tsv      per-sample call + probabilities + QC
   $OUTDIR/ancestry_pca.png         cohort plotted over the reference
   $OUTDIR/cohort_pca.eigenvec      <- GWAS covariates come from HERE
   $OUTDIR/merged_pca.eigenvec      ancestry assignment only, not covariates
   $TMP/                            PLINK logs and intermediates

 Check the CV accuracy printed above (want >= 0.98) and look at the PNG
 before trusting the table.
==============================================================
EOF
