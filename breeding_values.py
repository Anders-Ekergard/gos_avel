"""
Breeding value per individual (step M in the pipeline sketch) - combines
the markers' additive effects (step K, marker_phenotype_association.py)
into one aggregate value per animal.

NOTE: toy data (n=20) - shows the PRINCIPLE for how a breeding value is
calculated from marker effects, no conclusion about pikeperch.

Run: python breeding_values.py
"""

from marker_phenotype_association import GENOTYPE_CSV, PHENOTYPE_CSV, read_genotypes, read_phenotypes, mean

# Allele dosage: number of "B" alleles at the marker (0 = AA, 1 = AB, 2 = BB)
ALLELE_DOSAGE = {"AA": 0, "AB": 1, "BB": 2}


def marker_effects(genotypes: dict[str, dict[str, str]], phenotypes: dict[str, float]) -> dict[str, dict]:
    """For each marker: additive effect per "B" allele (half the difference
    between the homozygote groups), plus how well the heterozygote follows an
    additive model (distance from the expected midpoint (AA+BB)/2)."""
    animal_ids = [d for d in genotypes if d in phenotypes]
    marker_names = next(iter(genotypes.values())).keys()

    effects = {}
    for marker in marker_names:
        groups: dict[str, list[float]] = {"AA": [], "AB": [], "BB": []}
        for animal in animal_ids:
            groups[genotypes[animal][marker]].append(phenotypes[animal])

        if not groups["AA"] or not groups["BB"]:
            effects[marker] = {"effect": 0.0, "additive_deviation": None}
            continue

        mean_aa = mean(groups["AA"])
        mean_bb = mean(groups["BB"])
        effect = (mean_bb - mean_aa) / 2  # change in FCR per "B" allele

        additive_deviation = None
        if groups["AB"]:
            expected_ab = (mean_aa + mean_bb) / 2
            additive_deviation = mean(groups["AB"]) - expected_ab

        effects[marker] = {"effect": effect, "additive_deviation": additive_deviation}
    return effects


def calculate_breeding_values(
    genotypes: dict[str, dict[str, str]], effects: dict[str, dict]
) -> dict[str, float]:
    """Per individual: the sum of allele dosage * marker effect over all
    markers. Lower value = lower expected FCR = better, since feed
    conversion ratio should be minimized."""
    values = {}
    for animal, animal_genotype in genotypes.items():
        value = 0.0
        for marker, genotype in animal_genotype.items():
            dosage = ALLELE_DOSAGE[genotype]
            value += dosage * effects[marker]["effect"]
        values[animal] = value
    return values


def main():
    genotypes = read_genotypes(GENOTYPE_CSV)
    phenotypes = read_phenotypes(PHENOTYPE_CSV)

    effects = marker_effects(genotypes, phenotypes)

    print("--- Marker effects (additive model: AA vs BB) ---")
    for marker, e in effects.items():
        warning = ""
        if e["additive_deviation"] is not None and abs(e["additive_deviation"]) > 0.15:
            warning = "  <- heterozygote does NOT follow the additive model well, use with caution"
        print(f"  {marker}: effect/B-allele={e['effect']:+.3f}{warning}")
    print()

    values = calculate_breeding_values(genotypes, effects)

    print("--- Breeding value per animal (lower = better, expected FCR contribution) ---")
    for animal, value in sorted(values.items(), key=lambda x: x[1]):
        print(f"  {animal}: breeding value={value:+.3f}  (actual FCR: {phenotypes[animal]:.2f})")

    print()
    print("NOTE: all markers are counted, even those whose heterozygote doesn't follow")
    print("an additive model well (flagged above) - their contribution is therefore noise,")
    print("not signal. In a real selection, such markers should be down-weighted or")
    print("excluded until more data confirms them.")


if __name__ == "__main__":
    main()
