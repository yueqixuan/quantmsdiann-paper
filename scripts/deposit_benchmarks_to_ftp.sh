#!/bin/bash
# Deposit the pending quantmsdiann reanalysis reports to the public PRIDE FTP,
# in the layout the figure generators expect:
#   quantmsdiann-benchmarks/<class>/<dataset>/v<diann_version>/quant_tables/...
#
# Run on a host that has BOTH the cluster scratch (read) and the FTP area
# (write) mounted. DRY-RUN by default; set APPLY=1 to actually copy.
#
#   APPLY=1 WHAT=quant_tables bash scripts/deposit_benchmarks_to_ftp.sh
#
# WHAT=quant_tables (default) copies only quant_tables/ + run_metadata.json
#   (the reproducibility-essential set — keeps ProCan from shipping 79 GB).
# WHAT=full copies the entire results dir (mirrors the single-cell deposition).
set -euo pipefail

SRC=/hps/nobackup/juan/pride/reanalysis/absolute-expression
DST=/nfs/ftp/public/databases/pride/resources/proteomes/quantmsdiann-benchmarks
APPLY="${APPLY:-0}"
WHAT="${WHAT:-quant_tables}"

# "<src-relative-to-$SRC> | <dst-relative-to-$DST>"   (only trees that exist)
ENTRIES=(
  "single-cell-proteomics/MSV000093870/results-plexDIA|single-cell/MSV000093870/v2_5_1"
  "cell-lines-proteomes/PXD030304/results|cell-lines/PXD030304/v2_5_1"
  "spatial-proteomics/PXD064049/results-CHP212-MYCN-DVP-diaPASEF-carbamidomethyl|spatial/PXD064049/v2_5_0"
  "phospho-proteomics/PXD049692/results-NK-fibrin-IL15-phospho-diaPASEF|phospho/PXD049692/v2_5_1"
  "_v2_5_1_enterprise/phospho-proteomics/PXD049692/results-NK-fibrin-IL15-phospho-diaPASEF|phospho/PXD049692/v2_5_1_enterprise"
  "phospho-proteomics/PXD034128/results-phospho-biological-study|phospho/PXD034128-biological-study/v2_5_1"
  "_v2_5_1_enterprise/phospho-proteomics/PXD034128/results-phospho-biological-study|phospho/PXD034128-biological-study/v2_5_1_enterprise"
  "phospho-proteomics/PXD034128/results-phospho-highspeed-DIA|phospho/PXD034128-highspeed-DIA/v2_5_1"
  "_v2_5_1_enterprise/phospho-proteomics/PXD034128/results-phospho-highspeed-DIA|phospho/PXD034128-highspeed-DIA/v2_5_1_enterprise"
  "phospho-proteomics/PXD034623/results-M2-Galectin1-phospho-DIA|phospho/PXD034623/v2_5_1"
  "_v2_5_1_enterprise/phospho-proteomics/PXD034623/results-M2-Galectin1-phospho-DIA|phospho/PXD034623/v2_5_1_enterprise"
  # plexDIA enterprise tree is only ~3 MB (incomplete run) -- intentionally omitted.
)

DRY="--dry-run"; [ "$APPLY" = "1" ] && DRY=""
echo "MODE: $([ "$APPLY" = 1 ] && echo APPLY || echo DRY-RUN)   WHAT: $WHAT"

for e in "${ENTRIES[@]}"; do
  s="$SRC/${e%%|*}"; d="$DST/${e##*|}"
  if [ ! -d "$s" ]; then echo "SKIP (missing): $s"; continue; fi
  echo "=== $s  ->  $d ==="
  if [ "$WHAT" = "full" ]; then
    [ "$APPLY" = 1 ] && mkdir -p "$d"
    rsync -a $DRY --info=progress2 --exclude='work' "$s/" "$d/"
  else
    [ "$APPLY" = 1 ] && mkdir -p "$d/quant_tables"
    rsync -a $DRY --info=progress2 "$s/quant_tables/" "$d/quant_tables/"
    [ -f "$s/run_metadata.json" ] && rsync -a $DRY "$s/run_metadata.json" "$d/"
  fi
done
echo "Done ($([ "$APPLY" = 1 ] && echo copied || echo dry-run))."
