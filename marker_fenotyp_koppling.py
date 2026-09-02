"""
Ovning: koppla markorgenotyp till fenotyp (FCR) i leksaksdatan
(genotyper_exempel.csv / fenotyper_exempel.csv).

Visar principen bakom en enkel markor-fenotyp-koppling - grundiden i
QTL-kartlaggning/GWAS: gruppera individer efter genotyp vid en markor,
jamfor medelfenotyp mellan grupperna. En markor vars grupper skiljer sig
tydligt at ar en QTL-kandidat.

OBS: leksaksdata (n=20), bara for att visa principen - ingen slutsats om
gos. Se aven Del 2-diskussionen om GRM/OCS.

Kor: python marker_fenotyp_koppling.py
"""

import csv
import io
from pathlib import Path

GENOTYP_CSV = Path(__file__).parent / "genotyper_exempel.csv"
FENOTYP_CSV = Path(__file__).parent / "fenotyper_exempel.csv"


def parse_genotyper(csv_text: str) -> dict[str, dict[str, str]]:
    """Som las_genotyper, men fran en CSV-strang istallet for en fil pa
    disk - anvands av webb-API:t for uppladdad/inklistrad CSV-data."""
    genotyper = {}
    reader = csv.DictReader(io.StringIO(csv_text))
    markor_kolumner = [k for k in (reader.fieldnames or []) if k != "djur_id"]
    for rad in reader:
        genotyper[rad["djur_id"]] = {m: rad[m] for m in markor_kolumner}
    return genotyper


def las_genotyper(path: Path) -> dict[str, dict[str, str]]:
    """djur_id -> {markornamn: genotyp (AA/AB/BB)}"""
    with open(path, newline="", encoding="utf-8") as f:
        return parse_genotyper(f.read())


def parse_fenotyper(csv_text: str, kolumn: str = "fcr") -> dict[str, float]:
    """Som las_fenotyper, men fran en CSV-strang istallet for en fil pa
    disk - anvands av webb-API:t for uppladdad/inklistrad CSV-data."""
    fenotyper = {}
    reader = csv.DictReader(io.StringIO(csv_text))
    if kolumn not in (reader.fieldnames or []):
        raise ValueError(f"Kolumnen '{kolumn}' finns inte - tillgangliga kolumner: {reader.fieldnames}")
    for rad in reader:
        fenotyper[rad["djur_id"]] = float(rad[kolumn])
    return fenotyper


def las_fenotyper(path: Path, kolumn: str = "fcr") -> dict[str, float]:
    """djur_id -> varde i angiven fenotyp-kolumn (standard "fcr", for
    bakatkompatibilitet med leksaksdatan). Byt kolumn for att lanka mot en
    annan fenotyp-CSV utan att andra sjalva lasfunktionen."""
    with open(path, newline="", encoding="utf-8") as f:
        return parse_fenotyper(f.read(), kolumn)


def medel(varden: list[float]) -> float:
    return sum(varden) / len(varden)


def koppling_per_markor(genotyper: dict[str, dict[str, str]], fenotyper: dict[str, float]) -> None:
    djur_ids = [d for d in genotyper if d in fenotyper]
    markor_namn = next(iter(genotyper.values())).keys()

    print(f"{len(djur_ids)} djur med både genotyp och fenotyp\n")

    for markor in markor_namn:
        grupper: dict[str, list[float]] = {"AA": [], "AB": [], "BB": []}
        for djur in djur_ids:
            genotyp = genotyper[djur][markor]
            grupper[genotyp].append(fenotyper[djur])

        print(f"--- {markor} ---")
        medelvarden = {}
        for genotyp, varden in grupper.items():
            if varden:
                m = medel(varden)
                medelvarden[genotyp] = m
                print(f"  {genotyp}: n={len(varden):<2} medel-FCR={m:.3f}")
            else:
                print(f"  {genotyp}: n=0")

        if "AA" in medelvarden and "BB" in medelvarden:
            skillnad = medelvarden["BB"] - medelvarden["AA"]
            if skillnad < 0:
                tolkning = "BB har lägre (bättre) FCR än AA"
            elif skillnad > 0:
                tolkning = "AA har lägre (bättre) FCR än BB"
            else:
                tolkning = "ingen skillnad"
            print(f"  Skillnad BB - AA: {skillnad:+.3f}  ({tolkning})")
        print()

    print("OBS: leksaksdata, n=20 uppdelat på upp till 3 grupper - för lite för")
    print("statistisk säkerhet (ingen signifikanstest görs här med avsikt).")
    print("Poängen är att visa PRINCIPEN: en markör vars genotypgrupper har")
    print("tydligt olika medel-FCR är en QTL-kandidat, vars koppling till")
    print("fysiologi sedan kan tolkas med Del 1-metodiken (kandidatgen-dossiern).")


def main():
    genotyper = las_genotyper(GENOTYP_CSV)
    fenotyper = las_fenotyper(FENOTYP_CSV)
    koppling_per_markor(genotyper, fenotyper)


if __name__ == "__main__":
    main()
