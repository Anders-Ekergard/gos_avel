"""
Kandidatgen-dossier: slar upp en gen hos en eller flera arter mot UniProt,
Ensembl, NCBI SRA/GEO och EVA - for att se vad som redan finns publikt
tillgangligt innan egen provtagning eller synteny-analys planeras.

Kor: python candidate_gene_dossier.py
"""

import requests


def clean_input(prompt: str) -> str:
    """input() men stripper bort ev. UTF-8 BOM som vissa terminaler/pipes
    lagger till forst i indata."""
    return input(prompt).replace("﻿", "").strip()

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
ENSEMBL_BASE = "https://rest.ensembl.org"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EVA_BASE = "https://www.ebi.ac.uk/eva/webservices/rest/v1"

# Kandidatgener namnda i Sun & Zhu (2019), "Designing future farmed fishes
# using genome editing", uppdelade per fenotyp - en kategori i taget, i
# linje med projektets prioritetsordning (tillvaxt/foderkvot forst).
#
# fat1/fat2 (kottkvalitet, omega-3) ar medvetet uteslutna: de ar
# "humaniserade" C. elegans-gener (fat-1/fat-2-desaturaser) inforda via
# transgenes i artikelns studie - de kan inte uppsta genom avel och hor
# inte hemma i en MAS/GRM/OCS-kandidatlista.
#
# igf1/igf2/ghr - tillagda 2026-08-27: till skillnad fran resten av listan
# (overford fran andra fiskarter via ortologi/bryggart) ar dessa redan
# direkt kopplade till tillvaxt HOS GOS SJALV, via en SNP-studie pa
# Sander lucioperca (se GOS_GENOMIK_RESURSER.TILLVAXT_SNP_KALLA i
# gos_genomik_resurser.py) - starkast bevisniva av hela listan, ingen
# artbrygga kravs for just dessa tre.
CANDIDATE_GENE_SETS: dict[str, list[str]] = {
    "tillvaxt": ["socs1a", "mstnb", "igf1", "igf2", "ghr"],
    "kon": ["foxl2", "cyp19a", "cyp17a", "dnd"],
    "kottkvalitet": ["sp7", "stat3"],
    "sjukdom_stress": ["fih", "vhl"],
}

# Staende bryggart: mycket av det funktionella beviset (knockout -> fenotyp)
# for fiskgener kommer fran zebrafisk, men zebrafisk (Cypriniformes) ligger
# ~230+ miljoner ar fran gos (Percomorpha). Tilapia (Cichliformes) ligger
# inom samma stora grupp som gos, ~100-120 miljoner ar bort - battre bevarad
# synteny/reglering, battre brygga for att bedoma om ett zebrafisk-fynd
# sannolikt galler percider ocksa.
BRIDGE_SPECIES = "Oreochromis niloticus"


def uniprot_gene_hits(gene_symbol: str, organism: str) -> list[dict]:
    """Sok UniProtKB efter en gensymbol hos en art.

    Provar forst exakt match. Ger den noll traffar, provas de vanligaste
    teleost-paralogsuffixen a/b (samma monster som gh1a/gh1b, mstna/mstnb -
    och nu igf2b: gos saknar en ren "igf2" men har "igf2b", upptackt
    2026-08-27 via manuell UniProt-sokning som exakt match missade helt).

    Ett oppet suffix-wildcard (gene:X*) provades forst men gav for mycket
    brus - t.ex. matchade "ghr" avsiktligt in bade ghrl (ghrelin) och ghrhr
    (GHRH-receptorn), helt andra hormongener. Darfor bara a/b, aldrig ett
    fritt "*", och det faktiska matchade gennamnet returneras alltid sa en
    breddad trafff aldrig ar tyst."""
    hits = _uniprot_query(gene_symbol, organism)
    if hits:
        return hits

    breddade_hits: list[dict] = []
    for suffix in ("a", "b"):
        breddade_hits += _uniprot_query(f"{gene_symbol}{suffix}", organism)
    return breddade_hits


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
    """Om flera UniProt-traffar delar samma Ensembl-gen ar de troligen
    isoformer av samma gen, inte skilda paraloger."""
    return {h["ensembl_gene"] for h in hits if h["ensembl_gene"]}


def ensembl_paralog_count(ensembl_gene: str, ensembl_species: str, timeout: int = 60) -> int:
    """Antal paraloger enligt Ensembl Compara.
    OBS: kan inkludera hela genfamiljen (avlagsna dupliceringar), inte bara
    nyligen uppkomna artspecifika kopior.

    timeout ar konfigurerbar (standard 60s, samma som tidigare) eftersom
    webb-API:t i api/pipeline.py kor mot en Vercel-funktion med 30s
    maxDuration for HELA anropet - dar behovs en kortare timeout sa ett
    trogt Ensembl-svar felar snabbt (fangas redan av anroparens except
    RequestException) istallet for att sjalvt sla i taket."""
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
    """Antal SRA-poster (RNA-seq) for en lista arter, ev. filtrerat pa vavnad."""
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
    """Antal GEO-poster for en lista arter, ev. filtrerat pa vavnad."""
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
    """Ar arten registrerad i EVA (variant-arkiv)?"""
    r = requests.get(f"{EVA_BASE}/meta/species/list", timeout=30)
    r.raise_for_status()
    arter = r.json().get("response", [{}])[0].get("result", [])
    return any(
        scientific_name.lower() in (a.get("taxonomyScientificName") or "").lower()
        for a in arter
    )


def gather_dossier(gene_symbol: str, species_list: list[str], paralog_timeout: int = 60) -> dict:
    """Samlar all dossier-data for en gen, utan att skriva ut nagot.

    Ren datainsamlingsfunktion - anvands bade av print_dossier (CLI) och
    av webb-API:t (api/pipeline.py) sa de tva ytorna aldrig kan divergera
    i vad de faktiskt rapporterar.

    paralog_timeout: se ensembl_paralog_count - CLI:t anvander standard 60s,
    webb-API:t skickar in ett kortare varde for att hinna inom Vercels
    30s maxDuration for hela requesten."""
    per_species = []
    for organism in species_list:
        hits = uniprot_gene_hits(gene_symbol, organism)
        genes = unique_ensembl_genes(hits)
        matchade_namn = sorted({h["matched_gene_name"] for h in hits if h["matched_gene_name"]})
        broadened = bool({n.lower() for n in matchade_namn} - {gene_symbol.lower()})

        paralogs: dict[str, int | str] = {}
        ensembl_species = organism.lower().replace(" ", "_")
        for gene in genes:
            try:
                paralogs[gene] = ensembl_paralog_count(gene, ensembl_species, timeout=paralog_timeout)
            except requests.exceptions.RequestException as e:
                paralogs[gene] = f"fel: {e.__class__.__name__}"

        per_species.append({
            "organism": organism,
            "hit_count": len(hits),
            "hits": hits,
            "unique_ensembl_genes": sorted(genes),
            "isoform_warning": len(hits) > len(genes) and bool(genes),
            "matched_gene_names": matchade_namn,
            "broadened_search": broadened,
            "paralogs": paralogs,
        })

    return {
        "gene_symbol": gene_symbol,
        "species": per_species,
        "sra_rna_seq_count": sra_rna_seq_count(species_list),
        "geo_dataset_count": geo_dataset_count(species_list),
        "eva_status": {organism: eva_species_registered(organism) for organism in species_list},
    }


def print_dossier(gene_symbol: str, species_list: list[str]) -> None:
    dossier = gather_dossier(gene_symbol, species_list)
    print(f"\n== Kandidatgen-dossier: {dossier['gene_symbol']} ==\n")

    for sp in dossier["species"]:
        print(f"--- {sp['organism']} ---")
        if not sp["hit_count"]:
            print("  UniProt: inga traffar")
            print()
            continue

        print(f"  UniProt: {sp['hit_count']} post(er), {len(sp['unique_ensembl_genes'])} unik(a) Ensembl-gen(er)")
        if sp["isoform_warning"]:
            print("    -> flera UniProt-poster delar Ensembl-gen: sannolikt isoformer, inte paraloger")
        if sp["broadened_search"]:
            print(f"    Faktiskt matchade gennamn: {', '.join(sp['matched_gene_names'])}"
                  f" (exakt match gav noll, a/b-suffix testades istället - kontrollera att rätt gen hittades)")

        for gene, count in sp["paralogs"].items():
            if isinstance(count, str):
                print(f"    {gene}: kunde inte hämta paraloger ({count})")
            else:
                print(f"    {gene}: {count} Ensembl-paralog(er) (kan vara hela genfamiljen - kolla innan slutsats dras)")
        print()

    print("--- Publik sekvensdata (alla arter i listan, alla vavnader) ---")
    print(f"  SRA RNA-seq-poster: {dossier['sra_rna_seq_count']}")
    print(f"  GEO-dataset: {dossier['geo_dataset_count']}")
    print()

    print("--- EVA (variant-arkiv) ---")
    for organism, status in dossier["eva_status"].items():
        print(f"  {organism}: {'registrerad' if status else 'ej hittad'}")


def batch_screen(gene_symbols: list[str], species_list: list[str]) -> None:
    """Snabb triage over manga kandidatgener: UniProt-narvaro + isoform/paralog-
    signal per gen. SRA/GEO/EVA racknas bara en gang for hela artlistan,
    eftersom de sokningarna inte ar gen-specifika."""
    print(f"\n== Batch-screening: {len(gene_symbols)} kandidatgener ==\n")

    print("--- Artnivå (gäller alla gener nedan - inte gen-specifikt) ---")
    print(f"  SRA RNA-seq-poster: {sra_rna_seq_count(species_list)}")
    print(f"  GEO-dataset: {geo_dataset_count(species_list)}")
    for organism in species_list:
        status = "registrerad" if eva_species_registered(organism) else "ej hittad"
        print(f"  EVA ({organism}): {status}")
    print()

    col_labels = [o.split()[-1] for o in species_list]
    header = f"{'Gen':<10} | " + " | ".join(f"{c:<12}" for c in col_labels) + " | Isoform? | Paraloger*"
    print(header)
    print("-" * len(header))

    breddade_traffar: dict[str, set[str]] = {}
    for gene in gene_symbols:
        counts = []
        isoform_flag = False
        primary_paralogs = None
        matchade_namn: set[str] = set()
        for i, organism in enumerate(species_list):
            hits = uniprot_gene_hits(gene, organism)
            genes = unique_ensembl_genes(hits)
            counts.append(len(hits))
            matchade_namn |= {h["matched_gene_name"] for h in hits if h["matched_gene_name"]}
            if len(hits) > len(genes) and genes:
                isoform_flag = True
            if i == 0 and genes:
                ensembl_species = organism.lower().replace(" ", "_")
                try:
                    primary_paralogs = ensembl_paralog_count(next(iter(genes)), ensembl_species)
                except requests.exceptions.RequestException:
                    primary_paralogs = "timeout"

        if {n.lower() for n in matchade_namn} - {gene.lower()}:
            breddade_traffar[gene] = matchade_namn

        counts_str = " | ".join(f"{c:<12}" for c in counts)
        isoform_str = "ja" if isoform_flag else ""
        paralog_str = str(primary_paralogs) if primary_paralogs is not None else "-"
        print(f"{gene:<10} | {counts_str} | {isoform_str:<8} | {paralog_str}")

    print(f"\n* Paraloger raknas bara for {species_list[0]} (forsta arten i listan),"
          " och kan vara hela genfamiljen - inte bara nyligen uppkomna kopior.")

    if breddade_traffar:
        print("\n--- Exakt match gav noll, a/b-suffix hittade andra gennamn - kontrollera ---")
        for gene, namn in breddade_traffar.items():
            print(f"  {gene} -> {', '.join(sorted(namn))}")


def main():
    mode = clean_input("Läge - (e)nskild gen eller (b)atch över kandidatgenlistan? [e/b]: ").lower()
    species_input = clean_input(
        "Arter, kommaseparerade (Enter för standard: "
        "Sander lucioperca, Sander vitreus, Perca fluviatilis): "
    )
    if species_input:
        species_list = [s.strip() for s in species_input.split(",")]
    else:
        species_list = ["Sander lucioperca", "Sander vitreus", "Perca fluviatilis"]

    if BRIDGE_SPECIES not in species_list:
        species_list.append(BRIDGE_SPECIES)

    if mode == "b":
        kategorier = list(CANDIDATE_GENE_SETS)
        kategori = clean_input(
            f"Fenotypkategori ({'/'.join(kategorier)}) [Enter = tillvaxt]: "
        ).lower()
        if not kategori:
            kategori = "tillvaxt"
        if kategori not in CANDIDATE_GENE_SETS:
            print(f"Okänd kategori '{kategori}', använder 'tillvaxt'.")
            kategori = "tillvaxt"
        batch_screen(CANDIDATE_GENE_SETS[kategori], species_list)
    else:
        gene_symbol = clean_input("Gensymbol att slå upp (t.ex. mstnb): ")
        print_dossier(gene_symbol, species_list)


if __name__ == "__main__":
    main()
