#!/usr/bin/env python3
"""
Classifie chaque espèce rencontrée dans un grand clade (Mammalia, Aves,
Actinopterygii, Insecta, Nematoda, autre) via l'API NCBI Taxonomy -- pas de
table de correspondance codée en dur, entièrement automatisé pour n'importe
quelle espèce présente dans miRBase.
"""
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLADES_OF_INTEREST = ["Mammalia", "Aves", "Actinopterygii", "Amphibia", "Insecta", "Nematoda"]


def fetch_lineage(species_name: str, retries: int = 3) -> str:
    search_url = f"{EUTILS}/esearch.fcgi?db=taxonomy&term={urllib.parse.quote(species_name)}&retmode=json"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(search_url, timeout=30) as resp:
                ids = json.loads(resp.read())["esearchresult"]["idlist"]
            if not ids:
                return "NA"
            fetch_url = f"{EUTILS}/efetch.fcgi?db=taxonomy&id={ids[0]}&retmode=xml"
            with urllib.request.urlopen(fetch_url, timeout=30) as resp:
                xml_data = resp.read()
            root = ET.fromstring(xml_data)
            lineage = root.findtext(".//Lineage") or ""
            return lineage
        except Exception:
            time.sleep(2)
    return "ECHEC"


def classify_clade(lineage: str) -> str:
    for clade in CLADES_OF_INTEREST:
        if clade in lineage:
            return clade
    return "Autre"


def main():
    species_file = DATA / "species_encountered.tsv"
    out_file = DATA / "species_taxonomy.tsv"

    species_list = []
    with open(species_file) as fh:
        next(fh)
        for line in fh:
            name, code = line.rstrip("\n").split("\t")
            species_list.append((name, code))

    print(f"Classification taxonomique de {len(species_list)} espèces via NCBI...")
    with open(out_file, "w") as out:
        out.write("species_full\tspecies_code\tclade\n")
        for i, (name, code) in enumerate(species_list):
            lineage = fetch_lineage(name.replace("_", " "))
            clade = classify_clade(lineage) if lineage not in ("NA", "ECHEC") else lineage
            out.write(f"{name}\t{code}\t{clade}\n")
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(species_list)} traitées")
            time.sleep(0.35)  # courtoisie API NCBI (max ~3 req/s sans clé)

    print(f"Terminé -> {out_file}")


if __name__ == "__main__":
    main()
