#!/usr/bin/env bash
# Enchaîne toutes les étapes du pipeline, du téléchargement miRBase au rapport final.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "### 1. Données miRBase ###"
if [ ! -f data/mature.fa ]; then
  curl -sL -o data/mature.fa https://www.mirbase.org/download/mature.fa
else
  echo "data/mature.fa déjà présent, téléchargement ignoré."
fi
if [ ! -f data/hairpin.fa ]; then
  curl -sL -o data/hairpin.fa https://www.mirbase.org/download/hairpin.fa
else
  echo "data/hairpin.fa déjà présent, téléchargement ignoré."
fi

echo
echo "### 2. Extraction par famille ###"
python3 scripts/00_extract_sequences.py

echo
echo "### 3. Classification taxonomique (API NCBI) ###"
python3 scripts/01_fetch_taxonomy.py

echo
echo "### 4. Alignement MAFFT + score de conservation (global + seed) ###"
python3 scripts/02_align_and_score.py

echo
echo "### 5. Recommandation de modèle animal ###"
python3 scripts/03_propose_model.py

echo
echo "### 6. Croisement orthologie Ensembl Compara ###"
python3 scripts/04_compara_crosscheck.py

echo
echo "### 7. Rapport final ###"
docker run --rm --user "$(id -u):$(id -g)" -v "$ROOT:/proj" -w /proj/scripts rocker/tidyverse:4.3.1 \
  Rscript -e "rmarkdown::render('05_report.Rmd')"
mv scripts/05_report.html results/report.html

echo
echo "Terminé -- voir results/report.html"
