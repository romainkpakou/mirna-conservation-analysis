# Plan — Conservation inter-espèces de microARN et choix de modèle animal

## Contexte

Reprise du projet original sur l'étude de microARN absents et de leur conservation
inter-espèces, avec établissement d'un profil automatisé et sélection d'un modèle animal
optimal. Toutes les données nécessaires sont publiques.

## Objectif scientifique

Étudier la conservation d'un ensemble de microARN d'intérêt à travers plusieurs espèces,
automatiser l'établissement d'un profil de conservation, et proposer un modèle animal
pertinent sur cette base.

## Données

| Source | Rôle | Accès |
|---|---|---|
| miRBase | Référence des miARN matures et précurseurs par espèce | `https://mirbase.org` |
| TargetScan | Conservation des sites cibles | `https://www.targetscan.org` |
| Ensembl Compara | Orthologie inter-espèces | `https://www.ensembl.org/info/genome/compara` |

## Méthodologie

1. Sélection d'un ensemble de microARN d'intérêt (ex : liés à une voie biologique donnée)
2. Extraction des séquences précurseurs/matures par espèce (dump/API miRBase)
3. Alignement multiple inter-espèces (`MAFFT` ou `ClustalO`) et calcul d'un score de
   conservation de séquence
4. Croisement avec l'orthologie Ensembl Compara pour confirmer les relations d'orthologie
5. Score de pression de sélection (conservation de séquence ; dN/dS si les régions codantes
   flanquantes s'y prêtent)
6. Sélection du modèle animal optimal selon le profil de conservation obtenu
7. **Automatisation** : le pipeline doit accepter n'importe quel miARN (ou liste) en entrée
   et produire le profil de conservation sans intervention manuelle

## Stack technique

`Python` · `R` · `MAFFT` / `ClustalO` · miRBase · Ensembl Compara · TargetScan

## Structure de dépôt prévue

```
mirna-conservation-analysis/
├── README.md
├── scripts/
│   ├── 01_fetch_mirbase.py
│   ├── 02_align_conservation.py
│   └── 03_report.py
├── data/
├── results/
└── LICENSE
```

## Livrables

- Pipeline automatisé et réutilisable pour tout microARN donné en entrée
- Rapport de conservation avec recommandation de modèle animal, justifiée par les scores
  obtenus

## Panel de miARN retenu (choisi le 2026-09-17)

miARN musculaires/cardiaques ("myomiRs") : *miR-1*, *miR-133a*, *miR-208a*, *miR-208b*,
*miR-499a*, plus *let-7a* comme référence de conservation profonde. Choisi pour sa
cohérence thématique avec les projets cardiaques précédents de la série
(cnv-clinical-pipeline, splicing-cardio-encode).

**Vérifié par coordonnées Ensembl (GRCh38)** : *MIR208A* (chr14:23 388 596-23 388 666) est
intronique dans *MYH6* (chr14:23 380 206-23 408 945) ; *MIR208B* (chr14:23 417 987-
23 418 063) est intronique dans *MYH7* (chr14:23 412 732-23 436 137) ; *MIR499A*
(chr20:34 990 376-34 990 497) est intronique dans *MYH7B* (chr20:34 955 810-35 002 440).
Confirmé par coordonnées réelles, pas une affirmation de mémoire non vérifiée.

## Écart au plan initial : identité de séquence plutôt que dN/dS

Les miARN matures sont non-codants ; la notion de substitutions synonymes/non-synonymes
(dN/dS) ne s'applique pas à leur séquence. Le score de conservation retenu est
l'**identité de séquence par rapport à l'orthologue humain** sur l'alignement MAFFT
(colonnes sans gap), mesure standard et appropriée pour des séquences non codantes courtes.

## Bug rencontré et corrigé : nomenclature miRBase par locus

*miR-1* existe en deux copies génomiques chez l'humain (*MIR1-1*/*MIR1-2*). miRBase
annote le produit mature du second locus `hsa-miR-1-2-3p` plutôt que `hsa-miR-1-3p` —
sans traitement particulier, la séquence de référence humaine aurait été absente du
fichier de la famille. Un repli automatique (`find_human_paralog_variant` dans
`00_extract_sequences.py`) détecte et gère ce cas, avec message explicite dans les logs.

## Statut

- [x] Sélection du panel de miARN (myomiRs + let-7a)
- [x] Extraction automatisée multi-espèces (miRBase, 160 espèces uniques rencontrées)
- [x] Classification taxonomique (API NCBI Taxonomy, automatisée par grand clade)
- [x] Alignement MAFFT + score de conservation (identité vs humain)
- [x] Recommandation de modèle animal par famille
- [x] Rapport final

## Améliorations post-livraison (2026-09-17)

Trois améliorations apportées après une relecture critique du projet terminé :

1. **Conservation de la seed region** (positions 2-8 de la séquence mature,
   `scripts/02_align_and_score.py`) : calculée séparément de l'identité globale. Résultat
   réel : la seed est systématiquement plus conservée (99,4-100 %) que la séquence
   entière (95,0-99,8 %) sur les 6 familles — cohérent avec la pression de sélection plus
   forte sur le segment déterminant pour la reconnaissance des cibles.
2. **Croisement Ensembl Compara** (`scripts/04_compara_crosscheck.py`) : comble un écart
   entre le plan initial (qui prévoyait Ensembl Compara) et la première implémentation
   (purement basée sur l'identité de séquence miRBase). Limité aux 3 miARN à locus
   génomique humain unique (*MIR208A*, *MIR208B*, *MIR499A*) — *miR-1*, *miR-133a* et
   *let-7a* ont plusieurs loci paralogues chez l'humain, un choix de gène "canonique"
   serait arbitraire pour ces trois-là, non traité plutôt que tranché arbitrairement.
   Résultat notable : confirme indépendamment l'absence d'orthologue pour *MIR499A* chez
   les modèles de laboratoire courants (0/4 espèces cibles), déjà constatée via miRBase.
3. **`scripts/run_all.sh`** : orchestration de bout en bout, cohérent avec les autres
   dépôts de la série.

### Statut

- [x] Seed region implémentée et vérifiée (résultat biologiquement cohérent)
- [x] Croisement Compara implémenté pour les 3 miARN à locus unique
- [x] `run_all.sh` testé de bout en bout
- [ ] TargetScan (mentionné au plan initial, toujours non utilisé) — les deux autres
      sources de données du plan (miRBase, Ensembl Compara) couvrent déjà l'essentiel du
      besoin ; TargetScan resterait pertinent pour une analyse de conservation des sites
      cibles dans les 3'UTR, hors du périmètre actuel (conservation de la séquence du
      miARN lui-même)
