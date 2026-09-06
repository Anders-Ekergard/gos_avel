"""
Pikeperch_genomic_resources: concrete, public genomic resources for
pikeperch (Sander lucioperca), collected from a German research program
(Rostock/Leibniz Institute) on "precision farming" for pikeperch. Pure
reference data - no functions, no API calls - the idea is that other
modules (candidate_gene_dossier.py, future GRM/OCS work) can import the
constants below instead of re-researching facts that are already
verified once.

Source: the "Broodstock" concept document (2026-08-25/27, see the memory
entry reference_broodstock_artifact for the URL), which in turn is based
on the publications listed per constant below.
"""

# GenBank accession for the reference genome, chromosome level.
REFERENCE_GENOME_ACCESSION = "GCA_008315115.2"  # SLUC_FBN_1.2
REFERENCE_GENOME_CHROMOSOMES = 24
REFERENCE_GENOME_SIZE_MB = 900  # approximate

# Ultra-dense SNP linkage map - usable as a marker source for a targeted
# MAS panel instead of running full whole-genome regression.
SNP_LINKAGE_MAP_MARKER_COUNT = 723_360
SNP_LINKAGE_MAP_SOURCE = (
    "Nguinkal et al. (2020), 'Investigation of the Genetic Diversity and "
    "Population Structure of Pikeperch (Sander lucioperca) Populations "
    "with a Novel Whole-Genome SNP Linkage Map', Scientific Reports"
)

# Raw whole-genome sequencing of half-sib families. NOTE: this is raw data
# (FASTQ), not a finished genotype table - requires its own alignment/
# variant-calling pipeline (e.g. BWA + GATK/DeepVariant) before it's usable
# here in the pipeline. Better to start with SNP_LINKAGE_MAP above.
SRA_PROJECT_HALF_SIB = "PRJNA626522"

# Published inbreeding finding from that same data - a concrete sanity-check
# target for a future GRM/OCS implementation: can our code compute
# something in this ballpark from real data?
INBREEDING_COEFFICIENT_F_OFFSPRING = 0.33
EFFECTIVE_POPULATION_SIZE_NE = 12  # far below the recommended minimum of ~50
INBREEDING_FINDING_SOURCE = (
    "Autozygosity in full-sib pikeperch families, Frontiers in Genetics "
    "(2022) - SRA PRJNA626522"
)

# SNPs already linked to growth IN PIKEPERCH ITSELF (not via an ortholog
# bridge from another species like the rest of CANDIDATE_GENE_SETS["growth"]).
GROWTH_SNP_GENES = ["igf1", "igf2", "ghr"]
GROWTH_SNP_SOURCE = (
    "IGF-I/IGF-II/GHR growth SNPs in pikeperch, Aquaculture International "
    "(2020)"
)

# Population structure study - a contrast point against the inbreeding
# finding above: despite F~0.33/Ne~12 in the half-sib data, domesticated
# pikeperch stocks overall do NOT have lower genetic diversity than wild
# ones, which is unusual for a domesticated species.
POPULATION_STRUCTURE_NUM_POPULATIONS = 21
POPULATION_STRUCTURE_NUM_FISH = 958
POPULATION_STRUCTURE_SOURCE = (
    "Genetic variation in wild vs. domesticated pikeperch, Animals (2022)"
)


if __name__ == "__main__":
    print(f"Reference genome: {REFERENCE_GENOME_ACCESSION} "
          f"({REFERENCE_GENOME_CHROMOSOMES} chromosomes, ~{REFERENCE_GENOME_SIZE_MB} Mb)")
    print(f"SNP linkage map: {SNP_LINKAGE_MAP_MARKER_COUNT:,} markers")
    print(f"  source: {SNP_LINKAGE_MAP_SOURCE}")
    print(f"SRA half-sib project: {SRA_PROJECT_HALF_SIB} (raw data, requires its own pipeline)")
    print(f"  published finding: F={INBREEDING_COEFFICIENT_F_OFFSPRING}, Ne={EFFECTIVE_POPULATION_SIZE_NE}")
    print(f"  source: {INBREEDING_FINDING_SOURCE}")
    print(f"Growth SNP genes (directly in pikeperch): {', '.join(GROWTH_SNP_GENES)}")
    print(f"  source: {GROWTH_SNP_SOURCE}")
    print(f"Population structure: {POPULATION_STRUCTURE_NUM_POPULATIONS} populations, "
          f"{POPULATION_STRUCTURE_NUM_FISH} fish")
    print(f"  source: {POPULATION_STRUCTURE_SOURCE}")
