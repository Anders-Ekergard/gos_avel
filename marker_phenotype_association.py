"""
Exercise: link marker genotype to phenotype (FCR) in the toy dataset
(genotypes_example.csv / phenotypes_example.csv).

Shows the principle behind a simple marker-phenotype association - the
basic idea behind QTL mapping/GWAS: group individuals by genotype at a
marker, compare mean phenotype between the groups. A marker whose groups
differ clearly is a QTL candidate.

NOTE: toy data (n=20), just to show the principle - no conclusion about
pikeperch. See also the Part 2 discussion on GRM/OCS.

Run: python marker_phenotype_association.py
"""

import csv
import io
from pathlib import Path

GENOTYPE_CSV = Path(__file__).parent / "genotypes_example.csv"
PHENOTYPE_CSV = Path(__file__).parent / "phenotypes_example.csv"


def parse_genotypes(csv_text: str) -> dict[str, dict[str, str]]:
    """Like read_genotypes, but from a CSV string instead of a file on
    disk - used by the web API for uploaded/pasted CSV data."""
    genotypes = {}
    reader = csv.DictReader(io.StringIO(csv_text))
    marker_columns = [k for k in (reader.fieldnames or []) if k != "animal_id"]
    for row in reader:
        genotypes[row["animal_id"]] = {m: row[m] for m in marker_columns}
    return genotypes


def read_genotypes(path: Path) -> dict[str, dict[str, str]]:
    """animal_id -> {marker_name: genotype (AA/AB/BB)}"""
    with open(path, newline="", encoding="utf-8") as f:
        return parse_genotypes(f.read())


def parse_phenotypes(csv_text: str, column: str = "fcr") -> dict[str, float]:
    """Like read_phenotypes, but from a CSV string instead of a file on
    disk - used by the web API for uploaded/pasted CSV data."""
    phenotypes = {}
    reader = csv.DictReader(io.StringIO(csv_text))
    if column not in (reader.fieldnames or []):
        raise ValueError(f"Column '{column}' not found - available columns: {reader.fieldnames}")
    for row in reader:
        phenotypes[row["animal_id"]] = float(row[column])
    return phenotypes


def read_phenotypes(path: Path, column: str = "fcr") -> dict[str, float]:
    """animal_id -> value in the given phenotype column (default "fcr", for
    backward compatibility with the toy data). Change column to link to
    a different phenotype CSV without changing the read function itself."""
    with open(path, newline="", encoding="utf-8") as f:
        return parse_phenotypes(f.read(), column)


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def association_per_marker(genotypes: dict[str, dict[str, str]], phenotypes: dict[str, float]) -> None:
    animal_ids = [d for d in genotypes if d in phenotypes]
    marker_names = next(iter(genotypes.values())).keys()

    print(f"{len(animal_ids)} animals with both genotype and phenotype\n")

    for marker in marker_names:
        groups: dict[str, list[float]] = {"AA": [], "AB": [], "BB": []}
        for animal in animal_ids:
            genotype = genotypes[animal][marker]
            groups[genotype].append(phenotypes[animal])

        print(f"--- {marker} ---")
        mean_values = {}
        for genotype, values in groups.items():
            if values:
                m = mean(values)
                mean_values[genotype] = m
                print(f"  {genotype}: n={len(values):<2} mean-FCR={m:.3f}")
            else:
                print(f"  {genotype}: n=0")

        if "AA" in mean_values and "BB" in mean_values:
            difference = mean_values["BB"] - mean_values["AA"]
            if difference < 0:
                interpretation = "BB has lower (better) FCR than AA"
            elif difference > 0:
                interpretation = "AA has lower (better) FCR than BB"
            else:
                interpretation = "no difference"
            print(f"  Difference BB - AA: {difference:+.3f}  ({interpretation})")
        print()

    print("NOTE: toy data, n=20 split across up to 3 groups - too little for")
    print("statistical confidence (no significance test is run here, by design).")
    print("The point is to show the PRINCIPLE: a marker whose genotype groups have")
    print("clearly different mean FCR is a QTL candidate, whose link to")
    print("physiology can then be interpreted with the Part 1 methodology (the candidate-gene dossier).")


def main():
    genotypes = read_genotypes(GENOTYPE_CSV)
    phenotypes = read_phenotypes(PHENOTYPE_CSV)
    association_per_marker(genotypes, phenotypes)


if __name__ == "__main__":
    main()
