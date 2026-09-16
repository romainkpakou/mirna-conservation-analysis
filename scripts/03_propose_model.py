#!/usr/bin/env python3
"""
Pour chaque famille de miARN, recommande le meilleur modèle animal de
laboratoire parmi un panel standard, sur la base du score de conservation
(identité de séquence vs humain) calculé à l'étape précédente.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# Organismes de laboratoire courants et leur code miRBase
LAB_MODELS = {
    "mmu": "Souris (Mus musculus)",
    "rno": "Rat (Rattus norvegicus)",
    "dre": "Poisson-zèbre (Danio rerio)",
    "gga": "Poulet (Gallus gallus)",
    "xtr": "Xénope (Xenopus tropicalis)",
    "dme": "Drosophile (Drosophila melanogaster)",
    "cel": "Nématode (Caenorhabditis elegans)",
}


def main():
    recommendations = []

    conservation_files = [
        f for f in sorted(RESULTS.glob("conservation_*.tsv")) if f.name != "conservation_summary.tsv"
    ]
    for conservation_file in conservation_files:
        family = conservation_file.stem.replace("conservation_", "")
        with open(conservation_file) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))

        available_models = [r for r in rows if r["species_code"] in LAB_MODELS]
        available_models.sort(key=lambda r: -float(r["pct_identity_vs_human"]))

        if available_models:
            best = available_models[0]
            recommendation = LAB_MODELS[best["species_code"]]
            identity = best["pct_identity_vs_human"]
            alternatives = ", ".join(
                f"{LAB_MODELS[r['species_code']]} ({r['pct_identity_vs_human']}%)"
                for r in available_models[1:4]
            )
        else:
            recommendation, identity, alternatives = "Aucun modèle de laboratoire courant disponible", "NA", ""

        recommendations.append({
            "famille": family,
            "n_modeles_disponibles": len(available_models),
            "modele_recommande": recommendation,
            "identite_pct": identity,
            "alternatives": alternatives,
        })
        print(f"{family} : {recommendation} ({identity}% identité)"
              + (f" -- alternatives : {alternatives}" if alternatives else ""))

    out_path = RESULTS / "model_recommendations.tsv"
    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "famille", "n_modeles_disponibles", "modele_recommande", "identite_pct", "alternatives"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(recommendations)
    print(f"\nTerminé -> {out_path}")


if __name__ == "__main__":
    main()
