#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd -P)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd -P)
cd "$PROJECT_ROOT"

mkdir -p analysis/meta

GTDB=gtdbtk.bac120.summary.tsv
TARGET_OUT=analysis/step1_select_target_mags
PRODIGAL_OUT=analysis/step2_prodigal
CARVEME_OUT=analysis/step3_carveme
SMETANA_OUT=analysis/step5_smetana_no_split

### Step1: select target MAGs
mkdir -p ${TARGET_OUT}/{genomes,logs}

awk -F '\t' 'BEGIN{OFS="\t"} NR==1{next} $2 ~ /s__Lactobacillus delbrueckii$/ || $2 ~ /s__Streptococcus thermophilus$/ {sp=$2; sub(/^.*s__/,"",sp); print $1,sp,$2,$3,$6,$7}' "$GTDB" | sort -k2,2 -k1,1 > analysis/meta/target_mags.tmp.tsv

: > "${TARGET_OUT}/logs/missing_fasta.txt"
{
  printf 'mag_id\tspecies\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\tsource_fasta\tselected_fasta\n'
  while IFS=$'\t' read -r mid sp class ref ani af; do
    src=''
    for suf in fa fna fasta fas; do
      if [ -s "dereplicated_genomes/${mid}.${suf}" ]; then
        src=$(readlink -f "dereplicated_genomes/${mid}.${suf}")
        break
      fi
    done
    if [ -z "$src" ]; then
      echo "$mid" >> "${TARGET_OUT}/logs/missing_fasta.txt"
      continue
    fi
    dst="${TARGET_OUT}/genomes/${mid}.fa"
    [ -s "$dst" ] || cp "$src" "$dst"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$mid" "$sp" "$class" "$ref" "$ani" "$af" "$src" "$(readlink -f "$dst")"
  done < analysis/meta/target_mags.tmp.tsv
} > analysis/meta/target_mags.tsv

cut -f8 analysis/meta/target_mags.tsv | tail -n +2 > analysis/meta/target_mag_fasta.list
cut -f1,2 analysis/meta/target_mags.tsv | tail -n +2 > analysis/meta/target_mag_species.tsv
{
  echo -e 'species\tn_mags'
  awk -F '\t' 'NR>1{c[$2]++} END{for (s in c) print s"\t"c[s]}' analysis/meta/target_mags.tsv | sort
} > analysis/meta/target_mag_counts.tsv

cat analysis/meta/target_mag_counts.tsv
echo "Selected $(tail -n +2 analysis/meta/target_mags.tsv | wc -l) MAG fasta files"

### Step2: run Prodigal
source activate mags

mkdir -p ${PRODIGAL_OUT}/{faa,ffn,gbk,gff,logs}
find ${TARGET_OUT}/genomes -maxdepth 1 -type f -name '*.fa' | sort > analysis/meta/target_mag_fasta.list

parallel -j ${PRODIGAL_JOBS:-8} '
  fa={}
  mid=$(basename "$fa" .fa)
  if [ -s '"${PRODIGAL_OUT}"'/faa/${mid}.faa ]; then
    echo "[SKIP] ${mid}"
  else
    echo "[RUN] ${mid}"
    prodigal -i "$fa" -a '"${PRODIGAL_OUT}"'/faa/${mid}.faa -d '"${PRODIGAL_OUT}"'/ffn/${mid}.ffn -o '"${PRODIGAL_OUT}"'/gbk/${mid}.gbk -f gbk -p meta > '"${PRODIGAL_OUT}"'/logs/${mid}.prodigal.log 2>&1
    prodigal -i "$fa" -o '"${PRODIGAL_OUT}"'/gff/${mid}.gff -f gff -p meta >> '"${PRODIGAL_OUT}"'/logs/${mid}.prodigal.log 2>&1
  fi
' :::: analysis/meta/target_mag_fasta.list > ${PRODIGAL_OUT}/logs/prodigal.parallel.log 2>&1

find ${PRODIGAL_OUT}/faa -maxdepth 1 -type f -name '*.faa' | sort > analysis/meta/target_faa.list
wc -l analysis/meta/target_faa.list > analysis/meta/target_faa.count.txt
awk '/^>/ {n++} END {print FILENAME"\t"n}' ${PRODIGAL_OUT}/faa/*.faa > analysis/meta/prodigal_gene_counts.tsv
ln -sfn step2_prodigal/faa analysis/faa_local
cat analysis/meta/target_faa.count.txt

### Step3: build CarveMe models
source activate gems

mkdir -p ${CARVEME_OUT}/models ${CARVEME_OUT}/logs
find ${PRODIGAL_OUT}/faa -maxdepth 1 -type f -name '*.faa' | sort > analysis/meta/faa_local.list

parallel -j ${CARVEME_JOBS:-8} '
  faa={}
  mid=$(basename "$faa" .fa)
  mid=${mid%.faa}
  if [ -s '"${CARVEME_OUT}"'/models/${mid}.xml ]; then
    echo "[SKIP] ${mid}"
  else
    echo "[RUN] ${mid}"
    carve "$faa" -g M9 -i M9 -o '"${CARVEME_OUT}"'/models/${mid}.xml > '"${CARVEME_OUT}"'/logs/${mid}.carve.log 2>&1
  fi
' :::: analysis/meta/faa_local.list > ${CARVEME_OUT}/logs/carveme.parallel.log 2>&1

find ${CARVEME_OUT}/models -maxdepth 1 -type f -name '*.xml' | sed 's#.*/##; s#.xml$##' | sort > analysis/meta/models_built.list
wc -l analysis/meta/models_built.list > analysis/meta/models_built.count.txt
ln -sfn step3_carveme/models analysis/models
cat analysis/meta/models_built.count.txt

### Step4: run SMETANA
mkdir -p ${SMETANA_OUT}/{logs,results}

find ${CARVEME_OUT}/models -maxdepth 1 -type f -name '*.xml' | sort > analysis/meta/model_xml.list
n=$(wc -l < analysis/meta/model_xml.list)
[ "$n" -ge 2 ] || { echo "Need at least 2 models, found $n" >&2; exit 1; }

mapfile -t models < analysis/meta/model_xml.list
out=${SMETANA_OUT}/results/smetana_no_split
if ls "${out}"* >/dev/null 2>&1; then
  echo "[SKIP] existing SMETANA outputs"
else
  /usr/bin/time -v smetana -d --molweight --flavor fbc2 -o "$out" "${models[@]}" > ${SMETANA_OUT}/logs/smetana.log 2> ${SMETANA_OUT}/logs/smetana.time.log
fi

ln -sfn step5_smetana_no_split/results analysis/smetana
find ${SMETANA_OUT}/results -maxdepth 1 -type f -printf '%f\t%s bytes\n' | sort > analysis/meta/smetana_outputs.tsv
cat analysis/meta/smetana_outputs.tsv
