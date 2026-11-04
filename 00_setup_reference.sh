#!/usr/bin/env bash
#
# 00_setup_reference.sh -- one-off preparation of the 1000 Genomes reference panel.
# Needs plink2 for this step only: PLINK 1.9 cannot read .pgen. Everything after this runs on PLINK 1.9.

set -euo pipefail

REFDIR=${REFDIR:-./1000G}
OUTDIR=${OUTDIR:-./ancestry_out}
PLINK2=${PLINK2:-plink2}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$REFDIR" "$OUTDIR"
cd "$REFDIR"

command -v "$PLINK2" >/dev/null || {
  echo "ERROR: '$PLINK2' not on PATH." >&2
  echo "       Download the single binary from https://www.cog-genomics.org/plink/2.0/" >&2
  echo "       It is needed for this script only." >&2
  exit 1
}

# download 
# Links from https://www.cog-genomics.org/plink/2.0/resources
get () {  # get <url> <outfile>
  if [[ -f "$2" ]]; then
    echo "  $2 already here, skipping"
  else
    echo "  downloading $2"
    wget --continue --show-progress -O "$2" "$1" || curl -L -o "$2" "$1"
  fi
}

echo "[1/3] downloading the reference panel (~3.5 GB)"
get "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1" all_phase3.pgen.zst
get "https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1" all_phase3.pvar.zst
get "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1" all_phase3.psam

# convert 
echo
echo "[2/3] converting to PLINK 1 binary"
if [[ ! -f all_phase3.pgen ]]; then
  "$PLINK2" --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
fi

"$PLINK2" --pfile all_phase3 vzs \
          --autosome --snps-only just-acgt --max-alleles 2 --min-alleles 2 \
          --maf 0.01 \
          --make-bed --out 1000G_ref

echo "      1000G_ref: $(wc -l < 1000G_ref.fam) samples, $(wc -l < 1000G_ref.bim) variants"

# labels 
echo
echo "[3/3] population labels"
cd "$HERE"
python3 "$HERE/make_labels.py" --refdir "$REFDIR" --out "$OUTDIR/ref_labels.tsv"

cat <<EOF

==============================================================
 Reference ready.
   $REFDIR/1000G_ref.bed/.bim/.fam
   $OUTDIR/ref_labels.tsv

 Expected label counts: AFR 661, AMR 347, EAS 504, EUR 503, SAS 489.

 Next:  TARGET=ATP_data bash run_ancestry.sh
==============================================================
EOF
