# mirna-conservation-analysis

[![Docker](https://img.shields.io/badge/container-Docker-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Conservation inter-espèces de microARN musculaires/cardiaques et choix de modèle animal**

> Auteur : Romain KPAKOU | Master 2 Bioinformatique, Nantes Université
> GitHub : [github.com/romainkpakou](https://github.com/romainkpakou)

---

## Contexte

Étude de la conservation d'un panel de microARN musculaires/cardiaques ("myomiRs") —
*miR-1*, *miR-133a*, *miR-208a*, *miR-208b*, *miR-499a* — et de *let-7a* comme référence
de conservation profonde, avec établissement d'un profil automatisé et recommandation de
modèle animal. Panel cohérent avec le fil cardiaque des projets précédents de cette série
(*MYH7* dans le panel [cnv-clinical-pipeline](https://github.com/romainkpakou/cnv-clinical-pipeline),
RNA-seq cardiaque dans [splicing-cardio-encode](https://github.com/romainkpakou/splicing-cardio-encode)) :
*miR-208a*, *miR-208b* et *miR-499a* sont nichés dans des introns des gènes de myosine
cardiaque *MYH6*, *MYH7* et *MYH7B* respectivement (vérifié par coordonnées Ensembl
GRCh38, voir [PLAN.md](PLAN.md)).

Toutes les données sont publiques (miRBase, NCBI Taxonomy). Pipeline **entièrement
automatisé et paramétré par liste de miARN** — reproductible pour n'importe quel autre
miARN présent dans miRBase sans modification du code.

---

## Pipeline

```mermaid
flowchart TD
    A["miRBase (mature.fa)<br/>74 000+ séquences, toutes espèces"] --> B["1. Extraction par famille<br/>paramétrée (n'importe quel miARN)"]
    B --> C["2. Classification taxonomique<br/>API NCBI Taxonomy (grands clades)"]
    B --> D["3. Alignement MAFFT<br/>par famille, toutes espèces"]
    D --> E["4. Score de conservation<br/>identité globale + seed (positions 2-8)"]
    C --> E
    E --> F["5. Recommandation de modèle animal<br/>parmi organismes de laboratoire courants"]
    E --> G["6. Croisement orthologie<br/>API Ensembl Compara"]
    F --> H["7. Rapport final"]
    G --> H
```

---

## Prérequis

| Outil | Installation |
|---|---|
| Docker | [docs.docker.com](https://docs.docker.com) (MAFFT en conteneur) |
| Python 3 | déjà présent sur la plupart des systèmes |

---

## Utilisation

```bash
./scripts/run_all.sh
```

Télécharge les données miRBase si absentes, puis enchaîne toutes les étapes. Détail :

```bash
# Télécharger les données miRBase (~15 Mo)
curl -sL -o data/mature.fa https://www.mirbase.org/download/mature.fa
curl -sL -o data/hairpin.fa https://www.mirbase.org/download/hairpin.fa

# 1. Extraction par famille (paramétrer MIRNA_FAMILIES en tête de script pour tout autre miARN)
python3 scripts/00_extract_sequences.py

# 2. Classification taxonomique (API NCBI)
python3 scripts/01_fetch_taxonomy.py

# 3. Alignement MAFFT + score de conservation (identité globale + seed region)
python3 scripts/02_align_and_score.py

# 4. Recommandation de modèle animal
python3 scripts/03_propose_model.py

# 5. Croisement orthologie Ensembl Compara (locus uniques uniquement)
python3 scripts/04_compara_crosscheck.py

# 6. Rapport final
docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd):/proj" -w /proj/scripts rocker/tidyverse:4.3.1 \
  Rscript -e "rmarkdown::render('05_report.Rmd')"
mv scripts/05_report.html results/report.html
```

---

## Résultats

| Famille | Espèces (miRBase) | Clades | Identité globale moy. | Identité seed moy. |
|---|---|---|---|---|
| let-7a-5p | 74 | Nématodes → mammifères | 97,7 % | 99,8 % |
| miR-1-3p | 119 | Nématodes → mammifères | 95,0 % | 99,4 % |
| miR-133a-3p | 55 | Poissons → mammifères | 99,8 % | 100,0 % |
| miR-208a-3p | 20 | Mammifères uniquement | 97,5 % | 100,0 % |
| miR-208b-3p | 23 | Mammifères uniquement | 97,0 % | 100,0 % |
| miR-499a-5p | 6 | Mammifères + 1 poisson | 99,2 % | 100,0 % |

**La seed region (positions 2-8, déterminante pour la reconnaissance des cibles) est
systématiquement plus conservée que la séquence mature entière** — cohérent avec une
pression de sélection purificatrice plus forte sur la seed.

**Souris recommandée pour 4 des 6 familles** (100 % d'identité). Cas notable :
*miR-1-3p* atteint 100 % d'identité même chez *C. elegans* (conservation sur >600
millions d'années).

**Croisement Ensembl Compara** (3 miARN à locus unique) : orthologie confirmée
(`ortholog_one2one`) pour *MIR208A* jusqu'au poisson-zèbre et *MIR208B* jusqu'à la souris.
**Aucun orthologue Compara trouvé pour *MIR499A*** chez la souris, le rat, le poulet ou
le poisson-zèbre — confirme indépendamment l'absence déjà constatée dans miRBase.

Voir [`results/report.html`](https://htmlpreview.github.io/?https://github.com/romainkpakou/mirna-conservation-analysis/blob/main/results/report.html)
pour le rapport complet.

## Limites

- Séquences matures très courtes (20-23 nt) : chaque substitution pèse 4-5 points de
  pourcentage d'identité — mesure grossière à cette échelle.
- Absence d'un orthologue dans miRBase ≠ absence biologique (lacune d'annotation possible) —
  corroboré indépendamment par Ensembl Compara pour *miR-499a-5p* (voir ci-dessus).
- Identité de séquence utilisée plutôt que dN/dS : les miARN matures sont non-codants,
  la notion de substitutions synonymes/non-synonymes ne s'applique pas (écart documenté
  au plan initial, voir [PLAN.md](PLAN.md)).
- Le croisement Compara ne couvre que 3 des 6 familles (locus génomique humain unique
  requis) et 4 espèces cibles — pas une confirmation exhaustive.

## Licence

MIT — voir [LICENSE](LICENSE).
