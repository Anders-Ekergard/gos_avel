"""
Gos_genomik_resurser: konkreta, publika genomikresurser for gos (Sander
lucioperca), samlade fran ett tyskt forskningsprogram (Rostock/Leibniz-
institutet) om "precision farming" for gos. Ren referensdata - inga
funktioner, inga API-anrop - tanken ar att andra moduler (candidate_gene_
dossier.py, framtida GRM/OCS-arbete) kan importera konstanterna nedan
istallet for att atersoka fakta som redan ar verifierade en gang.

Kalla: "Broodstock"-konceptdokumentet (2026-08-25/27, se minnesposten
reference_broodstock_artifact for URL), som i sin tur bygger pa de
publikationer som listas per konstant nedan.
"""

# GenBank-accession for referensgenomet, kromosomniva.
REFERENCE_GENOME_ACCESSION = "GCA_008315115.2"  # SLUC_FBN_1.2
REFERENCE_GENOME_CHROMOSOMES = 24
REFERENCE_GENOME_SIZE_MB = 900  # ungefarligt

# Ultratat SNP-lankningskarta - anvandbar som markorkalla for en riktad
# MAS-panel istallet for att kora full helgenomsregression.
SNP_LINKAGE_MAP_MARKER_COUNT = 723_360
SNP_LINKAGE_MAP_KALLA = (
    "Nguinkal et al. (2020), 'Investigation of the Genetic Diversity and "
    "Population Structure of Pikeperch (Sander lucioperca) Populations "
    "with a Novel Whole-Genome SNP Linkage Map', Scientific Reports"
)

# Rå helgenomsekvensering av halvsyskonfamiljer. OBS: detta ar rådata
# (FASTQ), inte en fardig genotyptabell - kraver en egen alignment/variant-
# calling-pipeline (t.ex. BWA + GATK/DeepVariant) innan den ar anvandbar
# har i pipelinen. Borja hellre med SNP_LINKAGE_MAP ovan.
SRA_PROJECT_HALVSYSKON = "PRJNA626522"

# Publicerat inavelsfynd fran just den datan - konkret sanity-check-mal
# for en framtida GRM/OCS-implementation: kan var kod rakna fram nagot i
# den har storleksordningen fran riktig data?
INAVELSKOEFFICIENT_F_AVKOMMA = 0.33
EFFEKTIV_POPULATIONSSTORLEK_NE = 12  # langt under rekommenderat minimum ~50
INAVELSFYND_KALLA = (
    "Autozygosity in full-sib pikeperch families, Frontiers in Genetics "
    "(2022) - SRA PRJNA626522"
)

# SNP:ar redan kopplade till tillvaxt HOS GOS SJALV (inte via ortologi-
# brygga fran en annan art som resten av CANDIDATE_GENE_SETS["tillvaxt"]).
TILLVAXT_SNP_GENER = ["igf1", "igf2", "ghr"]
TILLVAXT_SNP_KALLA = (
    "IGF-I/IGF-II/GHR growth SNPs in pikeperch, Aquaculture International "
    "(2020)"
)

# Populationsstrukturstudie - kontrastpunkt mot inavelsfyndet ovan: trots
# F~0.33/Ne~12 i halvsyskondatan har domesticerade gösbestånd i stort INTE
# lagre genetisk diversitet an vilda, vilket ar ovanligt for en
# domesticerad art.
POPULATIONSSTRUKTUR_ANTAL_POPULATIONER = 21
POPULATIONSSTRUKTUR_ANTAL_FISKAR = 958
POPULATIONSSTRUKTUR_KALLA = (
    "Genetic variation in wild vs. domesticated pikeperch, Animals (2022)"
)


if __name__ == "__main__":
    print(f"Referensgenom: {REFERENCE_GENOME_ACCESSION} "
          f"({REFERENCE_GENOME_CHROMOSOMES} kromosomer, ~{REFERENCE_GENOME_SIZE_MB} Mb)")
    print(f"SNP-lankningskarta: {SNP_LINKAGE_MAP_MARKER_COUNT:,} markorer")
    print(f"  kalla: {SNP_LINKAGE_MAP_KALLA}")
    print(f"SRA halvsyskonprojekt: {SRA_PROJECT_HALVSYSKON} (rådata, kraver egen pipeline)")
    print(f"  publicerat fynd: F={INAVELSKOEFFICIENT_F_AVKOMMA}, Ne={EFFEKTIV_POPULATIONSSTORLEK_NE}")
    print(f"  kalla: {INAVELSFYND_KALLA}")
    print(f"Tillvaxt-SNP-gener (direkt hos gos): {', '.join(TILLVAXT_SNP_GENER)}")
    print(f"  kalla: {TILLVAXT_SNP_KALLA}")
    print(f"Populationsstruktur: {POPULATIONSSTRUKTUR_ANTAL_POPULATIONER} populationer, "
          f"{POPULATIONSSTRUKTUR_ANTAL_FISKAR} fiskar")
    print(f"  kalla: {POPULATIONSSTRUKTUR_KALLA}")
