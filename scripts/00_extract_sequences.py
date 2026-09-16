#!/usr/bin/env python3
"""
Extrait, pour une liste de familles de microARN donnée en entrée, les
séquences matures de toutes les espèces disponibles dans miRBase. Pipeline
paramétré par MIRNA_FAMILIES -- accepte n'importe quel miARN sans
modification du code (objectif du plan initial).
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SEQ_DIR = DATA / "sequences"
SEQ_DIR.mkdir(exist_ok=True)

# Panel de miARN musculaires/cardiaques (myomiRs) + let-7 comme référence de
# conservation profonde -- voir PLAN.md pour la justification biologique.
MIRNA_FAMILIES = [
    "miR-1-3p",
    "miR-133a-3p",
    "miR-208a-3p",
    "miR-208b-3p",
    "miR-499a-5p",
    "let-7a-5p",
]


def parse_mirbase_fasta(path: Path):
    """Parse mature.fa/hairpin.fa -> liste de (species_code, mirna_name, species_full, seq)."""
    entries = []
    header, seq_lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, "".join(seq_lines)))

    parsed = []
    for header, seq in entries:
        # Format : "<code>-<mirna_name> MIMATxxxxxxx <Species full name> <mirna_name>"
        # Le code espèce miRBase fait toujours 3-4 lettres minuscules -- on ne peut
        # pas juste découper au premier "-" car les noms de miARN en contiennent aussi.
        m = re.match(r"^([a-z]{3,4})-(\S+)\s+(MIMAT\d+)\s+(.+?)\s+\S+$", header)
        if not m:
            continue
        species_code, mirna_name, mimat_id, species_full = m.groups()
        parsed.append({
            "species_code": species_code, "mirna_name": mirna_name,
            "mimat_id": mimat_id, "species_full": species_full, "seq": seq, "header": header,
        })
    return parsed


def find_human_paralog_variant(all_entries, family: str):
    """
    Repli pour les miARN dupliqués chez l'humain (ex: MIR1-1/MIR1-2) : miRBase
    nomme alors le produit mature avec un suffixe de locus
    (hsa-miR-1-2-3p plutôt que hsa-miR-1-3p). Sans ce repli, l'espèce de
    référence (humaine) serait absente du fichier de la famille.
    """
    m = re.match(r"^(.+?)-(\d[a-z]?p)$", family)  # base (ex: "miR-1") + bras (ex: "3p")
    if not m:
        return None
    base, arm = m.groups()
    pattern = re.compile(rf"^{re.escape(base)}-\d+-{re.escape(arm)}$")
    candidates = [e for e in all_entries if e["species_code"] == "hsa" and pattern.match(e["mirna_name"])]
    return candidates[0] if candidates else None


def main():
    all_entries = parse_mirbase_fasta(DATA / "mature.fa")
    print(f"{len(all_entries)} séquences matures totales dans miRBase")

    species_seen = {}
    for family in MIRNA_FAMILIES:
        matches = [e for e in all_entries if e["mirna_name"] == family]

        if not any(e["species_code"] == "hsa" for e in matches):
            fallback = find_human_paralog_variant(all_entries, family)
            if fallback:
                matches.append(fallback)
                print(f"! {family} : pas de 'hsa-{family}' exact -- utilisation du variant "
                      f"par locus '{fallback['species_code']}-{fallback['mirna_name']}' comme référence humaine")

        out_path = SEQ_DIR / f"{family}.fa"
        with open(out_path, "w") as out:
            for e in matches:
                out.write(f">{e['species_code']}|{e['species_full'].replace(' ', '_')}\n{e['seq']}\n")
                species_seen[e["species_full"]] = e["species_code"]
        print(f"{family} : {len(matches)} espèces -> {out_path}")

    # Table des espèces uniques rencontrées, pour la classification taxonomique (étape suivante)
    species_table = DATA / "species_encountered.tsv"
    with open(species_table, "w") as out:
        out.write("species_full\tspecies_code\n")
        for name, code in sorted(species_seen.items()):
            out.write(f"{name}\t{code}\n")
    print(f"{len(species_seen)} espèces uniques -> {species_table}")


if __name__ == "__main__":
    main()
