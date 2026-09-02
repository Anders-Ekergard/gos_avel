"""
Testskript: finns räkor (t.ex. Litopenaeus vannamei) i EVA (European Variation Archive)?
Kör: python eva_rakor_test.py
"""

import requests

BASE = "https://www.ebi.ac.uk/eva/webservices/rest/v1"

# Sökord vi letar efter i art- och studienamn
SOKORD = list(input("Ange sökord (kommaseparerade): ").split(","))  #["shrimp", "prawn", "penaeus", "litopenaeus", "vannamei", "monodon"]


def matchar(text):
    text = (text or "").lower()
    return any(ord_ in text for ord_ in SOKORD)


def sok_arter():
    print("== Söker i artlistan ==")
    r = requests.get(f"{BASE}/meta/species/list", timeout=30)
    r.raise_for_status()
    arter = r.json().get("response", [{}])[0].get("result", [])
    traff = [a for a in arter if matchar(a.get("taxonomyScientificName", "")) or matchar(a.get("taxonomyCommonName", ""))]
    if traff:
        for a in traff:
            print(" -", a.get("taxonomyScientificName"), "|", a.get("taxonomyCode"), a.get("assemblyCode"))
    else:
        print(" Inga räk-arter hittade i artlistan.")
    return traff


def sok_studier():
    print("\n== Söker bland alla registrerade studier ==")
    r = requests.get(f"{BASE}/meta/studies/all", timeout=30)
    r.raise_for_status()
    studier = r.json().get("response", [{}])[0].get("result", [])
    traff = [s for s in studier if matchar(s.get("name", "")) or matchar(s.get("description", ""))]
    if traff:
        for s in traff:
            print(" -", s.get("id"), "|", s.get("name"))
    else:
        print(" Inga räk-relaterade studier hittade bland studienamn/beskrivningar.")
    return traff


if __name__ == "__main__":
    arter = sok_arter()
    studier = sok_studier()

    print("\n== Resultat ==")
    print(f"Räk-arter hittade: {len(arter)}")
    print(f"Räk-studier hittade: {len(studier)}")