"""
Candidate-gene dossier: looks up a gene for one or more species against
UniProt, Ensembl, NCBI SRA/GEO and EVA - to see what's already publicly
known before planning any wet-lab work or synteny analysis.

Run: python candidate_gene_dossier.py
"""

from concurrent.futures import ThreadPoolExecutor

import requests


def clean_input(prompt: str) -> str:
    """input() but strips any UTF-8 BOM that some terminals/pipes
    prepend to input."""
    return input(prompt).replace("﻿", "").strip()

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
ENSEMBL_BASE = "https://rest.ensembl.org"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EVA_BASE = "https://www.ebi.ac.uk/eva/webservices/rest/v1"

# Candidate genes named in Sun & Zhu (2019), "Designing future farmed fishes
# using genome editing", split by phenotype - one category at a time, in
# line with the project's priority order (growth/feed conversion first).
#
# fat1/fat2 (meat quality, omega-3) are deliberately excluded: they are
# "humanized" C. elegans genes (fat-1/fat-2 desaturases) introduced via
# transgenesis in the article's study - they cannot arise through breeding
# and don't belong in a MAS/GRM/OCS candidate list.
#
# igf1/igf2/ghr - added 2026-08-27: unlike the rest of the list (carried
# over from other fish species via orthology/bridge species) these are
# already directly linked to growth IN PIKEPERCH ITSELF, via a SNP study
# on Sander lucioperca (see GROWTH_SNP_SOURCE in
# pikeperch_genomic_resources.py) - the strongest evidence level of the
# whole list, no species bridge needed for these three.
CANDIDATE_GENE_SETS: dict[str, list[str]] = {
    "growth": ["socs1a", "mstnb", "igf1", "igf2", "ghr"],
    "sex": ["foxl2", "cyp19a", "cyp17a", "dnd"],
    "meat_quality": ["sp7", "stat3"],
    "disease_stress": ["fih", "vhl"],
}

# Standing bridge species: much of the functional evidence (knockout ->
# phenotype) for fish genes comes from zebrafish, but zebrafish
# (Cypriniformes) sits ~230+ million years from pikeperch (Percomorpha).
# Tilapia (Cichliformes) sits within the same larger group as pikeperch,
# ~100-120 million years away - better preserved synteny/regulation, a
# better bridge for judging whether a zebrafish finding likely also
# applies to percids.
BRIDGE_SPECIES = "Oreochromis niloticus"


def uniprot_gene_hits(gene_symbol: str, organism: str) -> list[dict]:
    """Search UniProtKB for a gene symbol in a species.

    Tries an exact match first. If that gives zero hits, tries the most
    common teleost paralog suffixes a/b (same pattern as gh1a/gh1b,
    mstna/mstnb - and now igf2b: pikeperch has no plain "igf2" but does
    have "igf2b", discovered 2026-08-27 via a manual UniProt search that
    the exact match missed entirely).

    An open suffix wildcard (gene:X*) was tried first but produced too
    much noise - e.g. "ghr" deliberately also matched ghrl (ghrelin) and
    ghrhr (the GHRH receptor), entirely different hormone genes. Hence
    only a/b, never a free "*", and the actually matched gene name is
    always returned so a broadened hit is never silent."""
    hits = _uniprot_query(gene_symbol, organism)
    if hits:
        return hits

    broadened_hits: list[dict] = []
    for suffix in ("a", "b"):
        broadened_hits += _uniprot_query(f"{gene_symbol}{suffix}", organism)
    return broadened_hits


def _uniprot_query(gene_query: str, organism: str) -> list[dict]:
    params = {
        "query": f'(gene:{gene_query}) AND (organism_name:"{organism}")',
        "format": "json",
        "fields": "accession,gene_names,protein_name,xref_ensembl",
    }
    r = requests.get(f"{UNIPROT_BASE}/search", params=params, timeout=30)
    r.raise_for_status()

    hits = []
    for entry in r.json().get("results", []):
        ensembl_gene = None
        for xref in entry.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "Ensembl":
                for prop in xref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        ensembl_gene = prop.get("value")

        protein_desc = entry.get("proteinDescription", {})
        name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")
        if not name:
            submission_names = protein_desc.get("submissionNames", [])
            if submission_names:
                name = submission_names[0].get("fullName", {}).get("value")

        matched_gene_name = None
        genes_field = entry.get("genes", [])
        if genes_field:
            matched_gene_name = genes_field[0].get("geneName", {}).get("value")

        hits.append({
            "accession": entry.get("primaryAccession"),
            "protein_name": name,
            "ensembl_gene": ensembl_gene,
            "matched_gene_name": matched_gene_name,
        })
    return hits


def unique_ensembl_genes(hits: list[dict]) -> set[str]:
    """If several UniProt hits share the same Ensembl gene, they are
    likely isoforms of the same gene, not separate paralogs."""
    return {h["ensembl_gene"] for h in hits if h["ensembl_gene"]}


def ensembl_paralog_count(ensembl_gene: str, ensembl_species: str, timeout: int = 60) -> int:
    """Number of paralogs according to Ensembl Compara.
    NOTE: may include the whole gene family (ancient duplications), not
    just recently arisen species-specific copies.

    timeout is configurable (default 60s, same as before) because the
    web API in api/gene_dossier.py runs against a Vercel function with a
    30s maxDuration for the WHOLE request - there a shorter timeout is
    needed so a slow Ensembl response fails fast (already caught by the
    caller's except RequestException) instead of hitting the ceiling
    itself."""
    url = f"{ENSEMBL_BASE}/homology/id/{ensembl_species}/{ensembl_gene}"
    params = {"content-type": "application/json", "type": "paralogues", "format": "condensed"}
    r = requests.get(url, params=params, timeout=timeout)
    if r.status_code == 404:
        return 0
    r.raise_for_status()
    data = r.json().get("data", [])
    if not data:
        return 0
    return len(data[0].get("homologies", []))


def sra_rna_seq_count(species_list: list[str], tissue_terms: list[str] | None = None) -> int:
    """Number of SRA records (RNA-seq) for a list of species, optionally filtered by tissue."""
    organism_query = " OR ".join(f"{s}[Organism]" for s in species_list)
    term = f"({organism_query}) AND RNA-Seq[Strategy]"
    if tissue_terms:
        tissue_query = " OR ".join(f"{t}[All Fields]" for t in tissue_terms)
        term += f" AND ({tissue_query})"
    params = {"db": "sra", "term": term, "retmode": "json", "retmax": 0}
    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    return int(r.json()["esearchresult"]["count"])


def geo_dataset_count(species_list: list[str], tissue_terms: list[str] | None = None) -> int:
    """Number of GEO records for a list of species, optionally filtered by tissue."""
    organism_query = " OR ".join(f"{s}[Organism]" for s in species_list)
    term = f"({organism_query})"
    if tissue_terms:
        tissue_query = " OR ".join(f"{t}[All Fields]" for t in tissue_terms)
        term += f" AND ({tissue_query})"
    params = {"db": "gds", "term": term, "retmode": "json", "retmax": 0}
    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    return int(r.json()["esearchresult"]["count"])


def eva_species_registered(scientific_name: str) -> bool:
    """Is the species registered in EVA (the variant archive)?"""
    r = requests.get(f"{EVA_BASE}/meta/species/list", timeout=30)
    r.raise_for_status()
    species = r.json().get("response", [{}])[0].get("result", [])
    return any(
        scientific_name.lower() in (a.get("taxonomyScientificName") or "").lower()
        for a in species
    )


def _gather_one_species(gene_symbol: str, organism: str, paralog_timeout: int) -> dict:
    hits = uniprot_gene_hits(gene_symbol, organism)
    genes = unique_ensembl_genes(hits)
    matched_names = sorted({h["matched_gene_name"] for h in hits if h["matched_gene_name"]})
    broadened = bool({n.lower() for n in matched_names} - {gene_symbol.lower()})

    paralogs: dict[str, int | str] = {}
    ensembl_species = organism.lower().replace(" ", "_")
    for gene in genes:
        try:
            paralogs[gene] = ensembl_paralog_count(gene, ensembl_species, timeout=paralog_timeout)
        except requests.exceptions.RequestException as e:
            paralogs[gene] = f"error: {e.__class__.__name__}"

    return {
        "organism": organism,
        "hit_count": len(hits),
        "hits": hits,
        "unique_ensembl_genes": sorted(genes),
        "isoform_warning": len(hits) > len(genes) and bool(genes),
        "matched_gene_names": matched_names,
        "broadened_search": broadened,
        "paralogs": paralogs,
    }


def gather_dossier(gene_symbol: str, species_list: list[str], paralog_timeout: int = 60) -> dict:
    """Collects all dossier data for a gene, without printing anything.

    Pure data-gathering function - used by both print_dossier (CLI) and
    the web API (api/gene_dossier.py) so the two surfaces can never
    diverge in what they actually report.

    paralog_timeout: see ensembl_paralog_count - the CLI uses the default
    60s, the web API passes in a shorter value to fit within Vercel's 30s
    maxDuration for the whole request.

    Everything (the per-species lookups including the paralog calls, plus
    SRA/GEO/EVA) runs IN PARALLEL across threads, not sequentially - an
    earlier sequential version with an 8s paralog timeout took ~8-16s
    total locally, but timed out repeatedly in production on Vercel
    (a different network path to Ensembl, 2026-09-02, discovered by the
    user after igf1 gave a ReadTimeout on O. niloticus twice in a row).
    Parallelizing lets us raise paralog_timeout without the sum of all
    calls growing linearly - one slow service only blocks its own part,
    not the whole request."""
    with ThreadPoolExecutor(max_workers=len(species_list) * 2 + 2) as executor:
        species_futures = {
            organism: executor.submit(_gather_one_species, gene_symbol, organism, paralog_timeout)
            for organism in species_list
        }
        sra_future = executor.submit(sra_rna_seq_count, species_list)
        geo_future = executor.submit(geo_dataset_count, species_list)
        eva_futures = {organism: executor.submit(eva_species_registered, organism) for organism in species_list}

        per_species = [species_futures[organism].result() for organism in species_list]
        eva_status = {organism: eva_futures[organism].result() for organism in species_list}

    return {
        "gene_symbol": gene_symbol,
        "species": per_species,
        "sra_rna_seq_count": sra_future.result(),
        "geo_dataset_count": geo_future.result(),
        "eva_status": eva_status,
    }


def print_dossier(gene_symbol: str, species_list: list[str]) -> None:
    dossier = gather_dossier(gene_symbol, species_list)
    print(f"\n== Candidate-gene dossier: {dossier['gene_symbol']} ==\n")

    for sp in dossier["species"]:
        print(f"--- {sp['organism']} ---")
        if not sp["hit_count"]:
            print("  UniProt: no hits")
            print()
            continue

        print(f"  UniProt: {sp['hit_count']} record(s), {len(sp['unique_ensembl_genes'])} unique Ensembl gene(s)")
        if sp["isoform_warning"]:
            print("    -> several UniProt records share an Ensembl gene: likely isoforms, not paralogs")
        if sp["broadened_search"]:
            print(f"    Actually matched gene name(s): {', '.join(sp['matched_gene_names'])}"
                  f" (exact match gave zero, a/b suffix tried instead - verify the right gene was found)")

        for gene, count in sp["paralogs"].items():
            if isinstance(count, str):
                print(f"    {gene}: could not fetch paralogs ({count})")
            else:
                print(f"    {gene}: {count} Ensembl paralog(s) (may be the whole gene family - check before drawing conclusions)")
        print()

    print("--- Public sequence data (all species in the list, all tissues) ---")
    print(f"  SRA RNA-seq records: {dossier['sra_rna_seq_count']}")
    print(f"  GEO datasets: {dossier['geo_dataset_count']}")
    print()

    print("--- EVA (variant archive) ---")
    for organism, status in dossier["eva_status"].items():
        print(f"  {organism}: {'registered' if status else 'not found'}")


def batch_screen(gene_symbols: list[str], species_list: list[str]) -> None:
    """Quick triage across many candidate genes: UniProt presence +
    isoform/paralog signal per gene. SRA/GEO/EVA are only counted once
    for the whole species list, since those lookups aren't gene-specific."""
    print(f"\n== Batch screening: {len(gene_symbols)} candidate genes ==\n")

    print("--- Species level (applies to all genes below - not gene-specific) ---")
    print(f"  SRA RNA-seq records: {sra_rna_seq_count(species_list)}")
    print(f"  GEO datasets: {geo_dataset_count(species_list)}")
    for organism in species_list:
        status = "registered" if eva_species_registered(organism) else "not found"
        print(f"  EVA ({organism}): {status}")
    print()

    col_labels = [o.split()[-1] for o in species_list]
    header = f"{'Gene':<10} | " + " | ".join(f"{c:<12}" for c in col_labels) + " | Isoform? | Paralogs*"
    print(header)
    print("-" * len(header))

    broadened_hits: dict[str, set[str]] = {}
    for gene in gene_symbols:
        counts = []
        isoform_flag = False
        primary_paralogs = None
        matched_names: set[str] = set()
        for i, organism in enumerate(species_list):
            hits = uniprot_gene_hits(gene, organism)
            genes = unique_ensembl_genes(hits)
            counts.append(len(hits))
            matched_names |= {h["matched_gene_name"] for h in hits if h["matched_gene_name"]}
            if len(hits) > len(genes) and genes:
                isoform_flag = True
            if i == 0 and genes:
                ensembl_species = organism.lower().replace(" ", "_")
                try:
                    primary_paralogs = ensembl_paralog_count(next(iter(genes)), ensembl_species)
                except requests.exceptions.RequestException:
                    primary_paralogs = "timeout"

        if {n.lower() for n in matched_names} - {gene.lower()}:
            broadened_hits[gene] = matched_names

        counts_str = " | ".join(f"{c:<12}" for c in counts)
        isoform_str = "yes" if isoform_flag else ""
        paralog_str = str(primary_paralogs) if primary_paralogs is not None else "-"
        print(f"{gene:<10} | {counts_str} | {isoform_str:<8} | {paralog_str}")

    print(f"\n* Paralogs are only counted for {species_list[0]} (the first species in the list),"
          " and may be the whole gene family - not just recently arisen copies.")

    if broadened_hits:
        print("\n--- Exact match gave zero, a/b suffix found other gene names - verify ---")
        for gene, names in broadened_hits.items():
            print(f"  {gene} -> {', '.join(sorted(names))}")


def main():
    mode = clean_input("Mode - (s)ingle gene or (b)atch over the candidate-gene list? [s/b]: ").lower()
    species_input = clean_input(
        "Species, comma-separated (Enter for default: "
        "Sander lucioperca, Sander vitreus, Perca fluviatilis): "
    )
    if species_input:
        species_list = [s.strip() for s in species_input.split(",")]
    else:
        species_list = ["Sander lucioperca", "Sander vitreus", "Perca fluviatilis"]

    if BRIDGE_SPECIES not in species_list:
        species_list.append(BRIDGE_SPECIES)

    if mode == "b":
        categories = list(CANDIDATE_GENE_SETS)
        category = clean_input(
            f"Phenotype category ({'/'.join(categories)}) [Enter = growth]: "
        ).lower()
        if not category:
            category = "growth"
        if category not in CANDIDATE_GENE_SETS:
            print(f"Unknown category '{category}', using 'growth'.")
            category = "growth"
        batch_screen(CANDIDATE_GENE_SETS[category], species_list)
    else:
        gene_symbol = clean_input("Gene symbol to look up (e.g. mstnb): ")
        print_dossier(gene_symbol, species_list)


if __name__ == "__main__":
    main()
