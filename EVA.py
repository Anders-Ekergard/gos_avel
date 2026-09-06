"""
Test script: is shrimp (e.g. Litopenaeus vannamei) registered in EVA (European Variation Archive)?
Run: python EVA.py
"""

import requests

BASE = "https://www.ebi.ac.uk/eva/webservices/rest/v1"

# Search terms we're looking for in species and study names
SEARCH_TERMS = list(input("Enter search terms (comma-separated): ").split(","))  #["shrimp", "prawn", "penaeus", "litopenaeus", "vannamei", "monodon"]


def matches(text):
    text = (text or "").lower()
    return any(term in text for term in SEARCH_TERMS)


def search_species():
    print("== Searching the species list ==")
    r = requests.get(f"{BASE}/meta/species/list", timeout=30)
    r.raise_for_status()
    species = r.json().get("response", [{}])[0].get("result", [])
    hits = [a for a in species if matches(a.get("taxonomyScientificName", "")) or matches(a.get("taxonomyCommonName", ""))]
    if hits:
        for a in hits:
            print(" -", a.get("taxonomyScientificName"), "|", a.get("taxonomyCode"), a.get("assemblyCode"))
    else:
        print(" No shrimp species found in the species list.")
    return hits


def search_studies():
    print("\n== Searching all registered studies ==")
    r = requests.get(f"{BASE}/meta/studies/all", timeout=30)
    r.raise_for_status()
    studies = r.json().get("response", [{}])[0].get("result", [])
    hits = [s for s in studies if matches(s.get("name", "")) or matches(s.get("description", ""))]
    if hits:
        for s in hits:
            print(" -", s.get("id"), "|", s.get("name"))
    else:
        print(" No shrimp-related studies found among study names/descriptions.")
    return hits


if __name__ == "__main__":
    species = search_species()
    studies = search_studies()

    print("\n== Results ==")
    print(f"Shrimp species found: {len(species)}")
    print(f"Shrimp studies found: {len(studies)}")
