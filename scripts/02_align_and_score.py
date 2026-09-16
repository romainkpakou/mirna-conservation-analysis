#!/usr/bin/env python3
"""
Aligne (MAFFT) les séquences matures de chaque famille de miARN à travers les
espèces disponibles, puis calcule un score de conservation par identité de
séquence par rapport à l'orthologue humain -- pas de dN/dS ici : les miARN
matures sont non-codants, la notion de substitutions synonymes/non-synonymes
ne s'applique pas (écart documenté au plan initial, voir PLAN.md).
"""
import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SEQ_DIR = DATA / "sequences"
ALIGN_DIR = DATA / "alignments"
RESULTS = ROOT / "results"
ALIGN_DIR.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

MAFFT_IMAGE = "quay.io/biocontainers/mafft:7.525--h031d066_1"


def run_mafft(input_fa: Path, output_fa: Path):
    cmd = [
        "docker", "run", "--rm", "--user", f"{__import__('os').getuid()}:{__import__('os').getgid()}",
        "-v", f"{input_fa.parent}:/data",
        MAFFT_IMAGE, "mafft", "--auto", "--quiet", f"/data/{input_fa.name}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    output_fa.write_text(result.stdout)


def read_fasta(path: Path) -> dict:
    seqs, header, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                seqs[header] = "".join(buf)
            header = line[1:]
            buf = []
        else:
            buf.append(line)
    if header:
        seqs[header] = "".join(buf)
    return seqs


def percent_identity(seq_a: str, seq_b: str) -> tuple[float, int]:
    """Identité en % sur les colonnes où aucune des deux séquences n'a de gap."""
    compared, matches = 0, 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a.upper() == b.upper():
            matches += 1
    return (100 * matches / compared if compared else 0.0), compared


# Positions 2-8 de la séquence mature (1-indexé) -- la "seed region", segment le
# plus déterminant pour la reconnaissance des ARNm cibles. Sa conservation
# peut différer de l'identité globale de la séquence mature.
SEED_START, SEED_END = 2, 8


def seed_alignment_columns(human_aligned_seq: str) -> list[int]:
    """Colonnes de l'alignement correspondant aux positions non-gappées 2-8
    de la séquence humaine (référence)."""
    columns, ungapped_pos = [], 0
    for col, base in enumerate(human_aligned_seq):
        if base == "-":
            continue
        ungapped_pos += 1
        if SEED_START <= ungapped_pos <= SEED_END:
            columns.append(col)
        if ungapped_pos > SEED_END:
            break
    return columns


def seed_percent_identity(human_seq: str, other_seq: str, seed_columns: list[int]) -> float:
    """Identité en % sur la seed uniquement (un gap chez l'espèce comparée
    compte comme un mismatch : l'absence de base à cette position n'est pas
    une conservation)."""
    if not seed_columns:
        return float("nan")
    matches = sum(
        1 for col in seed_columns
        if other_seq[col] != "-" and other_seq[col].upper() == human_seq[col].upper()
    )
    return 100 * matches / len(seed_columns)


def load_taxonomy() -> dict:
    tax = {}
    with open(DATA / "species_taxonomy.tsv") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tax[row["species_code"]] = row
    return tax


def main():
    taxonomy = load_taxonomy()
    summary_rows = []

    for fasta_path in sorted(SEQ_DIR.glob("*.fa")):
        family = fasta_path.stem
        n_seqs = fasta_path.read_text().count(">")
        if n_seqs < 3:
            print(f"{family} : trop peu de séquences ({n_seqs}) pour un alignement informatif, ignoré")
            continue

        aligned_path = ALIGN_DIR / f"{family}_aligned.fa"
        print(f"{family} : alignement MAFFT de {n_seqs} séquences...")
        run_mafft(fasta_path, aligned_path)

        seqs = read_fasta(aligned_path)
        human_header = next((h for h in seqs if h.startswith("hsa|")), None)
        if human_header is None:
            print(f"! {family} : pas de séquence humaine trouvée, impossible de calculer l'identité de référence")
            continue
        human_seq = seqs[human_header]
        seed_cols = seed_alignment_columns(human_seq)
        if len(seed_cols) < (SEED_END - SEED_START + 1):
            print(f"! {family} : séquence humaine plus courte que la seed region attendue (positions {SEED_START}-{SEED_END})")

        rows = []
        for header, seq in seqs.items():
            code, species_full = header.split("|", 1)
            pct_id, n_compared = percent_identity(human_seq, seq)
            pct_seed = seed_percent_identity(human_seq, seq, seed_cols)
            tax_row = taxonomy.get(code, {})
            rows.append({
                "species_code": code,
                "species_full": species_full.replace("_", " "),
                "clade": tax_row.get("clade", "NA"),
                "pct_identity_vs_human": round(pct_id, 1),
                "pct_identity_seed_vs_human": round(pct_seed, 1) if pct_seed == pct_seed else "NA",  # NaN check
                "positions_comparees": n_compared,
            })

        rows.sort(key=lambda r: -r["pct_identity_vs_human"])
        out_path = RESULTS / f"conservation_{family}.tsv"
        with open(out_path, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=[
                "species_code", "species_full", "clade", "pct_identity_vs_human",
                "pct_identity_seed_vs_human", "positions_comparees"
            ], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        clades_present = sorted(set(r["clade"] for r in rows if r["clade"] not in ("NA", "hsa")))
        seed_values = [r["pct_identity_seed_vs_human"] for r in rows if r["pct_identity_seed_vs_human"] != "NA"]
        summary_rows.append({
            "famille": family, "n_especes": len(rows),
            "clades_presents": ";".join(clades_present),
            "identite_moyenne": round(sum(r["pct_identity_vs_human"] for r in rows) / len(rows), 1),
            "identite_seed_moyenne": round(sum(seed_values) / len(seed_values), 1) if seed_values else "NA",
        })
        print(f"{family} : {len(rows)} espèces alignées -> {out_path}")

    with open(RESULTS / "conservation_summary.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "famille", "n_especes", "clades_presents", "identite_moyenne", "identite_seed_moyenne"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print("Terminé.")


if __name__ == "__main__":
    main()
