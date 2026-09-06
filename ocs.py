"""
OCS optimization (step N in the pipeline sketch) - simplified Optimal
Contribution Selection. Combines breeding value (step M) and genomic
similarity/GRM (step L) to propose breeding pairs: maximize genetic gain
while avoiding pairing the most closely related individuals.

NOTE: toy data (n=20). A real OCS (Meuwissen 1997) solves a continuous
optimization problem (quadratic programming) over how much each
individual should contribute to the next generation, under the
constraint that mean relatedness is kept under a target. Here a greedy
approximation is used that shows the same PRINCIPLE at a smaller scale:
pick the best individuals (lowest breeding value), pair them two by two
so that no pair exceeds a relatedness threshold.

Run: python ocs.py
"""

from itertools import combinations

from marker_phenotype_association import GENOTYPE_CSV, PHENOTYPE_CSV, read_genotypes, read_phenotypes
from breeding_values import marker_effects, calculate_breeding_values
from pikeperch_genomic_resources import (
    INBREEDING_COEFFICIENT_F_OFFSPRING as PUBLISHED_F,
    EFFECTIVE_POPULATION_SIZE_NE as PUBLISHED_NE,
    SRA_PROJECT_HALF_SIB,
)

# Future validation target for this module (not implemented yet): a
# published study on real pikeperch half-sib data (SRA project, see
# SRA_PROJECT_HALF_SIB above) found F=0.33, Ne=12 (see PUBLISHED_F/
# PUBLISHED_NE above) - seriously low, far below the recommended
# minimum of ~50. If/when genotypes/phenotypes from that data are made
# available here (requires its own alignment/variant-calling pipeline on
# raw FASTQ, see pikeperch_genomic_resources.py), a concrete sanity check
# is to run this module's f_projection/genomic_similarity logic against
# the real data and see if it reproduces the same order of magnitude -
# not just against the toy data (n=20) used now.

# How many of the best individuals (by breeding value) become breeding candidates
NUM_CANDIDATES = 8

# Fraction of shared allele dosages (identity-by-state) above which a pair
# is avoided, regardless of how good their breeding values are
RELATEDNESS_THRESHOLD = 0.7

# How many generations the F projection projects forward
GENERATIONS = 5


def genomic_similarity(genotypes: dict[str, dict[str, str]], ind1: str, ind2: str) -> float:
    """Fraction of markers where two individuals have identical genotype -
    the same principle as grm.py, but over all markers in the CSV data
    instead of the hardcoded A/B/C example."""
    markers = genotypes[ind1].keys()
    shared = sum(1 for m in markers if genotypes[ind1][m] == genotypes[ind2][m])
    return shared / len(markers)


def propose_pairs(
    candidates: list[str],
    genotypes: dict[str, dict[str, str]],
    breeding_values: dict[str, float],
    relatedness_threshold: float = RELATEDNESS_THRESHOLD,
) -> list[tuple[str, str, float]]:
    """Greedy pairing: take the individuals in breeding-value order, pair
    each unpaired individual with the best unpaired individual that is
    NOT above the relatedness threshold. If none exists, pick the least
    related."""
    remaining = sorted(candidates, key=lambda d: breeding_values[d])
    pairs: list[tuple[str, str, float]] = []

    while len(remaining) >= 2:
        first = remaining.pop(0)
        best_partner = None
        for candidate in remaining:
            if genomic_similarity(genotypes, first, candidate) <= relatedness_threshold:
                best_partner = candidate
                break
        if best_partner is None:
            best_partner = min(remaining, key=lambda k: genomic_similarity(genotypes, first, k))
        remaining.remove(best_partner)
        similarity = genomic_similarity(genotypes, first, best_partner)
        pairs.append((first, best_partner, similarity))

    return pairs


def naive_pairing(candidates: list[str], breeding_values: dict[str, float]) -> list[tuple[str, str]]:
    """Comparison baseline with no relatedness consideration: sort by
    breeding value and pair sequentially (best with second-best, etc.).
    Shows what the OCS step actually adds compared to just chasing
    breeding value."""
    ordered = sorted(candidates, key=lambda d: breeding_values[d])
    return [(ordered[i], ordered[i + 1]) for i in range(0, len(ordered) - 1, 2)]


def f_projection(mean_f_today: float, generations: int, ne: int) -> list[float]:
    """Standard approximation from quantitative genetics: rate of change
    in inbreeding per generation ~ 1/(2*Ne). F then grows according to
    F_t = F_(t-1) + (1 - F_(t-1)) * delta_F.

    NOTE: Ne is estimated crudely as the number of breeding animals
    actually used (no sex ratio in the toy data) - a real Ne is often
    lower than that (skewed sex ratio/family size lowers it). The
    projection also assumes that the same strategy and the same Ne
    repeat every generation - a trend, not a simulation of actual
    future pairs."""
    delta_f = 1 / (2 * ne)
    projection = []
    f = mean_f_today
    for _ in range(generations):
        f = f + (1 - f) * delta_f
        projection.append(f)
    return projection


def run_breeding_pipeline(
    genotypes: dict[str, dict[str, str]],
    phenotypes: dict[str, float],
    num_candidates: int = NUM_CANDIDATES,
    relatedness_threshold: float = RELATEDNESS_THRESHOLD,
    generations: int = GENERATIONS,
) -> dict:
    """Runs the full Part 2 chain (breeding value -> GRM/genomic similarity
    -> OCS pairing proposal -> F projection) and returns structured data,
    without printing anything.

    Pure data-gathering function - shared by the CLI (main) and the web
    API, the same pattern as candidate_gene_dossier.gather_dossier, so the
    two surfaces can never diverge in what they actually report."""
    effects = marker_effects(genotypes, phenotypes)
    values = calculate_breeding_values(genotypes, effects)

    candidates = sorted(values, key=lambda d: values[d])[:num_candidates]

    similarity_matrix = {
        f"{a}|{b}": genomic_similarity(genotypes, a, b)
        for a, b in combinations(candidates, 2)
    }

    pairs = propose_pairs(candidates, genotypes, values, relatedness_threshold)
    ocs_mean_f = sum(similarity for _, _, similarity in pairs) / len(pairs)

    naive_pairs = naive_pairing(candidates, values)
    naive_mean_f = sum(genomic_similarity(genotypes, a, b) for a, b in naive_pairs) / len(naive_pairs)

    ne = len(pairs) * 2
    projection = f_projection(ocs_mean_f, generations, ne)

    return {
        "marker_effects": {
            m: {"effect": e["effect"], "additive_deviation": e["additive_deviation"]}
            for m, e in effects.items()
        },
        "breeding_values": [
            {"animal": d, "breeding_value": values[d]} for d in sorted(values, key=lambda d: values[d])
        ],
        "candidates": candidates,
        "similarity_matrix": similarity_matrix,
        "relatedness_threshold": relatedness_threshold,
        "proposed_pairs": [
            {
                "a": a,
                "b": b,
                "f_proxy": similarity,
                "over_threshold": similarity > relatedness_threshold,
                "breeding_value_a": values[a],
                "breeding_value_b": values[b],
            }
            for a, b, similarity in pairs
        ],
        "ocs_mean_f": ocs_mean_f,
        "naive_mean_f": naive_mean_f,
        "ne": ne,
        "generations": generations,
        "f_projection": projection,
        "published_reference": {
            "f": PUBLISHED_F,
            "ne": PUBLISHED_NE,
            "sra_project": SRA_PROJECT_HALF_SIB,
        },
    }


def main():
    genotypes = read_genotypes(GENOTYPE_CSV)
    phenotypes = read_phenotypes(PHENOTYPE_CSV)

    result = run_breeding_pipeline(genotypes, phenotypes)

    print(f"--- {NUM_CANDIDATES} breeding candidates (lowest breeding value = best) ---")
    for animal in result["candidates"]:
        value = next(a["breeding_value"] for a in result["breeding_values"] if a["animal"] == animal)
        print(f"  {animal}: breeding value={value:+.3f}")
    print()

    print(f"--- Genomic similarity among the candidates (threshold: {RELATEDNESS_THRESHOLD}) ---")
    for pair_name, similarity in result["similarity_matrix"].items():
        a, b = pair_name.split("|")
        flag = "  <- above the threshold" if similarity > RELATEDNESS_THRESHOLD else ""
        print(f"  {a} & {b}: {similarity:.2f}{flag}")
    print()

    print("--- Proposed breeding pairs (OCS) ---")
    for p in result["proposed_pairs"]:
        warning = ""
        if p["over_threshold"]:
            warning = "  <- no partner under the threshold was available among the remaining candidates, least-related was chosen"
        print(
            f"  {p['a']} x {p['b']}  (F-proxy offspring: {p['f_proxy']:.2f}, "
            f"breeding values: {p['breeding_value_a']:+.3f} / {p['breeding_value_b']:+.3f}){warning}"
        )
    print()

    print("--- F comparison: OCS vs naive pairing (breeding value only) ---")
    print(f"  OCS mean F this generation:  {result['ocs_mean_f']:.3f}")
    print(f"  Naive mean F this generation: {result['naive_mean_f']:.3f}")
    if result["naive_mean_f"] > 0:
        reduction = (result["naive_mean_f"] - result["ocs_mean_f"]) / result["naive_mean_f"]
        print(f"  -> OCS lowers the mean by {reduction:.0%}")
    print()

    print(f"--- F projection, OCS strategy repeated (Ne~{result['ne']}) ---")
    for generation, f in enumerate(result["f_projection"], start=1):
        print(f"  Generation +{generation}: F~{f:.3f}")
    print()
    print("NOTE: rough trend (constant Ne, no actual simulation of future")
    print("genotypes) - see the f_projection docstring for assumptions.")
    print()
    ref = result["published_reference"]
    print(f"Reference point (real data, not this toy data): a published study on")
    print(f"pikeperch half-sibs (SRA {ref['sra_project']}) found F={ref['f']}, Ne={ref['ne']}.")


if __name__ == "__main__":
    main()
