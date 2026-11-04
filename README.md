# Genetic-Ancestry-Prediction
genetic ancestry prediction based on genetic data

## 1. Files

**Run these three by hand:**

| | |
|---|---|
| `00_setup_reference.sh` | once: download 1000 Genomes, convert, build labels |
| `run_ancestry.sh` | the analysis, steps 1–7 |
| `mds_stratification.sh` | optional: MDS output in the Marees-tutorial format |

**Called automatically — never run these directly:**

| | called by |
|---|---|
| `make_labels.py` | `00_setup_reference.sh` |
| `harmonise.py` | `run_ancestry.sh` |
| `classify_ancestry.py` | `run_ancestry.sh`|

All the shell scripts find their Python helpers through their own directory, so
**keep everything in one folder**.

---

## 2. Downloads

### Software

Single static binaries — unzip, put on `PATH`, no install.

| | Version | Link |
|---|---|---|
| **PLINK 1.9** — runs everything | 1.9.0, 27 Sep 2026 | <https://www.cog-genomics.org/plink/1.9/> |
| — Linux 64-bit | | <https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip> |
| — macOS | | <https://s3.amazonaws.com/plink1-assets/plink_mac_20260927.zip> |
| — Windows 64-bit | | <https://s3.amazonaws.com/plink1-assets/plink_win64_20260927.zip> |
| **PLINK 2.0** — setup step only | beta 1.3, 31 Oct 2026 | <https://www.cog-genomics.org/plink/2.0/> |
| — Linux 64-bit | | <https://s3.amazonaws.com/plink2-assets/beta1/plink2_linux_x86_64_20261031.zip> |
| — Linux AVX2 (faster) | | <https://s3.amazonaws.com/plink2-assets/beta1/plink2_linux_avx2_20261031.zip> |
| — macOS ARM (M1+) | | <https://s3.amazonaws.com/plink2-assets/beta1/plink2_mac_arm64_20261031.zip> |
| — Windows 64-bit | | <https://s3.amazonaws.com/plink2-assets/beta1/plink2_win64_20261031.zip> |

PLINK 2 is needed only because PLINK 1.9 cannot read `.pgen`. After setup you
never need it again.

```bash
pip install numpy pandas scikit-learn matplotlib scipy
```

### Reference panel — 1000 Genomes phase 3, b37, 2,504 samples

`00_setup_reference.sh` fetches these. Canonical source:
<https://www.cog-genomics.org/plink/2.0/resources>

| File | Size | Link |
|---|---|---|
| `all_phase3.pgen.zst` | 2.25 GiB | <https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1> |
| `all_phase3.pvar.zst` | 1.26 GiB | <https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1> |
| `phase3_corrected.psam` | small | <https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1> |

Build 37 is right here: `human_g1k_v37.fasta` and `HRC.r1-1.GRCh37` say the
cohort is b37. A build mismatch is caught in step 3 with a clear message.

If `phase3_corrected.psam` lacks `SuperPop`/`Population` columns,
`make_labels.py` says so; drop
`integrated_call_samples_v3.20130502.ALL.panel` from
<https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/> into `1000G/`
instead.

---

## 3. Run

```bash
bash 00_setup_reference.sh            # once, 15-30 min, mostly download
TARGET=ATP_data bash run_ancestry.sh  # a few minutes
```

Settings are environment variables:

| | default | |
|---|---|---|
| `TARGET` | `ATP_data` | bfile prefix, no extension |
| `REFDIR` / `OUTDIR` | `./1000G` / `./ancestry_out` | |
| `PLINK` / `PLINK2` | `plink` / `plink2` | |
| `NPC` / `NPC_GWAS` | 10 / 20 | PCs for ancestry / for GWAS covariates |
| `MEMORY` / `THREADS` | 8000 / 4 | MB; match your job's `--mem` |
| `HWE_P` | `1e-6` | applied **within** each reference population |
| `MIN_SHARED` | 20000 | abort below this many shared SNPs |

### On a cluster

```bash
cd path
export REFDIR="$PWD/1000G" OUTDIR="$PWD/out_file"
export TARGET="$PWD/data_file"
export PLINK="$(command -v plink)" PLINK2="$(command -v plink2)"
export MEMORY=12000 THREADS=4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LOKY_MAX_CPU_COUNT=4

DOWNLOAD_ONLY=1 bash 00_setup_reference.sh      # login node: network only

sbatch --job-name=ATP_ancestry --ntasks=1 --cpus-per-task=4 --mem=16G \
       --time=04:00:00 --chdir="$PWD" --export=ALL \
       --output="$OUTDIR/slurm-%j.log" \
       --wrap='SKIP_DOWNLOAD=1 PLINK2_OPTS="--memory 12000 --threads 4" \
               bash 00_setup_reference.sh && bash run_ancestry.sh'

## 4. Output

| | |
|---|---|
| `target_ancestry.tsv` | **the result** — one row per cohort sample |
| `ancestry_pca.png` | cohort over the reference clusters |
| `cohort_pca.eigenvec` | **GWAS covariates come from here** |
| `merged_pca.eigenvec` | ancestry assignment only |
| `qc.imiss` / `qc.het` | call rate / heterozygosity |
| `ref_labels.tsv` | 1000 Genomes sample → population → superpopulation |
| `tmp/` | PLINK logs and intermediates |

### `target_ancestry.tsv` columns

| column | meaning |
|---|---|
| `best_superpopulation` | the random forest's pick |
| `probability` | its probability for that group |
| `nearest_by_distance` | the closest group by Mahalanobis distance |
| `call` | **the usable answer**: a group, `Admixed/Uncertain`, or `Outlier` |
| `group_fit` | `inside` / `edge` / `outside` |
| `note` | set when the forest and the distances disagree, or when the sample is inside more than one group |
| `d_AFR` … `d_SAS` | Mahalanobis distance to **every** group |
| `rel_AFR` … `rel_SAS` | the same, divided by that group's cut-off. **< 1 means within it** |
| `p_AFR` … `p_SAS` | full probability vector |
| `missing_rate`, `low_callrate` | genotyping call rate |
| `het_F`, `het_z`, `het_outlier` | heterozygosity, standardised within ancestry group |
| `PC1` … `PCn` | joint-PCA coordinates |

Groups: **AFR** African, **AMR** Admixed American, **EAS** East Asian,
**EUR** European, **SAS** South Asian.

Parse by column **name**, not position — the layout changed between versions:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
            {n[$h["call"]]++} END{for(k in n) printf "  %-20s %4d\n", k, n[k]}' \
    ancestry_out/target_ancestry.tsv | sort -k2 -rn
```


## Run on MobaXterm

```
module load apps/binapps/conda/miniforge3/26.9.4-0
conda activate atp_ancestry
(
set -euo pipefail

cd ...

export TARGET="$PWD/..."

for ext in bed bim fam; do
    test -s "$TARGET.$ext" || {
        echo "miss：$TARGET.$ext" >&2
        exit 1
    }
done

for file in 00_setup_reference.sh run_ancestry.sh \
            make_labels.py harmonise.py classify_ancestry.py; do
    test -s "$file" || {
        echo "miss：$file" >&2
        exit 1
    }
done

command -v plink
command -v plink2
python3 -c "import numpy, pandas, sklearn, matplotlib"

export PLINK="$(command -v plink)"
export PLINK2="$(command -v plink2)"

export REFDIR="$PWD/1000G"
export OUTDIR="$(mktemp -d "$PWD/ancestry_trial_XXXXXX")"

export THREADS=4
export MEMORY=12000
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LOKY_MAX_CPU_COUNT=4

sed 's/--pfile all_phase3 vzs/--pfile all_phase3 vzs --memory 12000 --threads 4/' \
00_setup_reference.sh > 00_setup_reference_csf.sh

JOB_ID=$(sbatch --parsable \
    --job-name=ATP_ancestry \
    --partition=multicore_small \
    --ntasks=1 \
    --cpus-per-task=4 \
    --mem=16G \
    --time=04:00:00 \
    --chdir="$PWD" \
    --export=ALL \
    --output="$OUTDIR/slurm-%j.log" \
    --wrap='bash 00_setup_reference_csf.sh && bash run_ancestry.sh')

echo "NO：$JOB_ID"
echo "Outdic：$OUTDIR"
echo "LOGDIC：$OUTDIR/slurm-$JOB_ID.log"
)
```
