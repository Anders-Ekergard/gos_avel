"""
Text-mining triage: searches PubMed/PubTator3 for genes co-mentioned with
a phenotype in zebrafish literature (Danio rerio - the fish species with
the richest functional literature), and can optionally expand hits via
STRING's text-mining channel plus KEGG pathways (the same two sources
already used for human genes in pathway_mvp.py, rewritten here for
zebrafish organism codes instead of hsa/Homo sapiens).

This does NOT replace manual review - the hits are a starting list to
verify (UniProt/Ensembl via candidate_gene_dossier.py) and paste into
CANDIDATE_GENE_SETS by hand, one category at a time.

Run: python text_mining_triage.py
"""

import functools
import time

import requests

from candidate_gene_dossier import CANDIDATE_GENE_SETS, EUTILS_BASE, clean_input

PUBTATOR_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
STRING_BASE = "https://string-db.org/api"
KEGG_BASE = "https://rest.kegg.jp"

# Choice of "bridge species" for text mining: zebrafish has by far the most
# functional literature for fish (see rationale in candidate_gene_dossier.py).
DISCOVERY_SPECIES = "Danio rerio"
DISCOVERY_TAXON_ID = 7955       # NCBI taxon id, required by STRING
DISCOVERY_KEGG_ORGANISM = "dre"  # KEGG's organism code for Danio rerio

# Search terms per phenotype category - English, since PubMed/PubTator is
# English-language literature. The keys must match CANDIDATE_GENE_SETS in
# candidate_gene_dossier.py so hits can be cross-referenced against
# already-known candidates.
PHENOTYPE_QUERY_TERMS: dict[str, list[str]] = {
    "growth": ["growth rate", "feed conversion ratio", "muscle growth", "body weight"],
    "sex": ["sex determination", "sex differentiation", "gonad development"],
    "meat_quality": ["muscle fiber", "fillet quality", "flesh quality"],
    "disease_stress": ["hypoxia tolerance", "stress response", "disease resistance", "immune response"],
}

MAX_PMIDS_PER_TERM = 50       # cap per search term - triage, not a full literature review
BIOCJSON_CHUNK_SIZE = 50      # number of PMIDs per annotations call
REQUEST_DELAY_SECONDS = 0.4   # politeness pause between calls (~3/sec, same norm as NCBI E-utilities)


def _polite_pause() -> None:
    """Simple politeness pause between calls to PubTator3/STRING/KEGG - all
    three are smaller, publicly funded services without a documented hard
    rate limit for this kind of call; better to be polite."""
    time.sleep(REQUEST_DELAY_SECONDS)


# --- PubTator3 (discovery) ---------------------------------------------------

def pubtator_search_pmids(query: str, max_results: int = MAX_PMIDS_PER_TERM) -> list[str]:
    """Search PubTator3 for articles matching a free-text query
    (e.g. '"feed conversion ratio" Danio rerio'). Returns PMIDs.

    NOTE: PubTator3 hasn't been called from this repo before - the field
    names below (`results`, `pmid`) are a best guess from public
    documentation. Verify empirically (print(r.json())) before trusting
    the hits blindly."""
    pmids: list[str] = []
    page = 1
    while len(pmids) < max_results:
        params = {"text": query, "page": page}
        r = requests.get(f"{PUBTATOR_BASE}/search/", params=params, timeout=30)
        r.raise_for_status()
        data = r.json()

        results = data.get("results", [])
        if not results:
            break
        for item in results:
            pmid = str(item.get("pmid") or item.get("_id", "")).replace("PMID:", "")
            if pmid and pmid not in pmids:
                pmids.append(pmid)
        page += 1
        _polite_pause()
    return pmids[:max_results]


def pubtator_fetch_annotations(pmids: list[str]) -> list[dict]:
    """Fetches BioC JSON with normalized gene/species annotations for a
    list of PMIDs, in chunks of BIOCJSON_CHUNK_SIZE."""
    documents: list[dict] = []
    for i in range(0, len(pmids), BIOCJSON_CHUNK_SIZE):
        chunk = pmids[i:i + BIOCJSON_CHUNK_SIZE]
        params = {"pmids": ",".join(chunk)}
        r = requests.get(f"{PUBTATOR_BASE}/publications/export/biocjson", params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        documents.extend(data.get("PubTator3", data.get("documents", [])))
        _polite_pause()
    return documents


def extract_gene_mentions(documents: list[dict], taxon_id: int = DISCOVERY_TAXON_ID) -> list[dict]:
    """Picks out Gene annotations from BioC documents, filtered to
    documents that also have a Species annotation matching taxon_id."""
    mentions: list[dict] = []
    for doc in documents:
        pmid = doc.get("pmid") or doc.get("id")
        species_ids = set()
        gene_annotations = []
        for passage in doc.get("passages", []):
            for ann in passage.get("annotations", []):
                infons = ann.get("infons", {})
                if infons.get("type") == "Species":
                    species_ids.add(infons.get("identifier"))
                elif infons.get("type") == "Gene":
                    gene_annotations.append({
                        "text": ann.get("text"),
                        "ncbi_gene_id": infons.get("identifier"),
                    })
        if str(taxon_id) not in species_ids:
            continue
        for g in gene_annotations:
            mentions.append({"pmid": pmid, **g})
    return mentions


def rank_gene_candidates(mentions: list[dict]) -> list[dict]:
    """Aggregates gene mentions into a ranked candidate list: most unique
    PMIDs first (broader support in the literature), then most total
    mentions.

    NOTE: co-occurrence in literature is correlation, not causation - see
    print_caveats()."""
    by_gene: dict[str, dict] = {}
    for m in mentions:
        key = m["ncbi_gene_id"] or m["text"]
        entry = by_gene.setdefault(key, {
            "ncbi_gene_id": m["ncbi_gene_id"],
            "symbols_seen": set(),
            "pmids": set(),
            "mention_count": 0,
        })
        entry["symbols_seen"].add(m["text"])
        entry["pmids"].add(m["pmid"])
        entry["mention_count"] += 1

    ranked = []
    for entry in by_gene.values():
        ranked.append({
            "ncbi_gene_id": entry["ncbi_gene_id"],
            "symbols_seen": sorted(entry["symbols_seen"]),
            "pmid_count": len(entry["pmids"]),
            "mention_count": entry["mention_count"],
            "example_pmids": sorted(entry["pmids"])[:5],
        })
    ranked.sort(key=lambda e: (e["pmid_count"], e["mention_count"]), reverse=True)
    return ranked


def build_queries_for_category(category: str, extra_terms: list[str] | None = None) -> list[str]:
    """One search query per synonym term (not a single combined boolean
    query) - easier to get right without knowing PubTator3's exact query
    syntax, and easier to debug term by term."""
    terms = list(PHENOTYPE_QUERY_TERMS.get(category, []))
    if extra_terms:
        terms += extra_terms
    return [f'"{term}" {DISCOVERY_SPECIES}' for term in terms]


def discover_candidates_for_category(category: str, extra_terms: list[str] | None = None) -> list[dict]:
    """Full discovery pipeline for a phenotype category: search -> fetch
    annotations -> filter -> rank."""
    queries = build_queries_for_category(category, extra_terms)

    all_pmids: list[str] = []
    for q in queries:
        try:
            pmids = pubtator_search_pmids(q)
            print(f"  '{q}': {len(pmids)} PMIDs")
            for p in pmids:
                if p not in all_pmids:
                    all_pmids.append(p)
        except requests.exceptions.RequestException as e:
            print(f"  '{q}': search failed ({e.__class__.__name__})")

    documents = pubtator_fetch_annotations(all_pmids)
    mentions = extract_gene_mentions(documents)
    return rank_gene_candidates(mentions)


# --- STRING (expansion from a known gene, text-mining channel) --------------

def string_get_string_id(gene_symbol: str, species_taxon_id: int = DISCOVERY_TAXON_ID) -> str | None:
    """Looks up STRING's internal ID for a gene symbol - required by
    interaction_partners (doesn't accept free symbols directly, to avoid
    ambiguous matching)."""
    params = {"identifiers": gene_symbol, "species": species_taxon_id, "limit": 1}
    r = requests.get(f"{STRING_BASE}/json/get_string_ids", params=params, timeout=30)
    r.raise_for_status()
    hits = r.json()
    if not hits:
        return None
    return hits[0].get("stringId")


def string_textmining_partners(
    gene_symbol: str,
    species_taxon_id: int = DISCOVERY_TAXON_ID,
    min_tscore: float = 0.4,
    limit: int = 20,
) -> list[dict]:
    """Functional partners of a known gene, filtered to STRING's
    text-mining channel (tscore) - i.e. partners supported by
    co-occurrence in the literature, not just database/experimental
    evidence."""
    string_id = string_get_string_id(gene_symbol, species_taxon_id)
    if not string_id:
        return []

    params = {"identifiers": string_id, "species": species_taxon_id, "limit": limit}
    r = requests.get(f"{STRING_BASE}/json/interaction_partners", params=params, timeout=30)
    r.raise_for_status()
    partners = r.json()

    result = []
    for p in partners:
        tscore = p.get("tscore")
        if tscore is not None and tscore >= min_tscore:
            result.append({
                "partner_symbol": p.get("preferredName_B"),
                "tscore": tscore,
                "combined_score": p.get("score"),
            })
    result.sort(key=lambda e: e["tscore"], reverse=True)
    return result


# --- KEGG (pathways for a known gene) ---------------------------------------
# Same two-step method (Entrez -> gene id, KEGG link/list) already used
# for human genes in pathway_mvp.py - here with "Danio rerio[orgn]"/"dre"
# instead of "Homo sapiens[orgn]"/"hsa".

@functools.lru_cache(maxsize=64)
def resolve_kegg_gene_id(symbol: str, organism: str = DISCOVERY_SPECIES) -> str | None:
    """Looks up a KEGG gene id (e.g. 'dre:napp') via NCBI Entrez esearch,
    the same eutils service already used in candidate_gene_dossier.py."""
    params = {"db": "gene", "term": f"{symbol}[sym] AND {organism}[orgn]", "retmode": "json"}
    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return f"{DISCOVERY_KEGG_ORGANISM}:{ids[0]}" if ids else None


def parse_kegg_tsv(text: str) -> list[tuple[str, str]]:
    """Parses KEGG's flat, tab-separated response format ('key\\tvalue'
    per line), shared by both the link and list operations."""
    pairs = []
    for line in text.strip().splitlines():
        if not line:
            continue
        key, value = line.split("\t", 1)
        pairs.append((key, value))
    return pairs


def fetch_kegg_pathway_ids(gene_id: str) -> list[str]:
    """Which KEGG pathways a gene belongs to."""
    r = requests.get(f"{KEGG_BASE}/link/pathway/{gene_id}", timeout=30)
    r.raise_for_status()
    return [pathway_id for _gene, pathway_id in parse_kegg_tsv(r.text)]


@functools.lru_cache(maxsize=8)
def _fetch_kegg_pathway_list(organism: str) -> tuple[tuple[str, str], ...]:
    """Fetches and caches an organism's full pathway id -> name list
    (the same response regardless of which gene triggered the lookup)."""
    r = requests.get(f"{KEGG_BASE}/list/pathway/{organism}", timeout=30)
    r.raise_for_status()
    return tuple(parse_kegg_tsv(r.text))


def fetch_kegg_pathway_names(pathway_ids: list[str], organism: str = DISCOVERY_KEGG_ORGANISM) -> dict[str, str]:
    """Names for KEGG pathway ids (the link response only gives ids, the
    list response is needed for the names)."""
    if not pathway_ids:
        return {}
    names_by_bare_id = dict(_fetch_kegg_pathway_list(organism))
    return {
        pid: names_by_bare_id[pid.removeprefix("path:")]
        for pid in pathway_ids
        if pid.removeprefix("path:") in names_by_bare_id
    }


def kegg_pathways_for_gene(gene_symbol: str) -> list[dict]:
    """High-level function: gene symbol -> list of KEGG pathways
    ({'id', 'name'}) for zebrafish. Empty list if Entrez doesn't find the
    gene (e.g. a misspelling/synonym that doesn't match the Entrez index)
    - doesn't crash."""
    gene_id = resolve_kegg_gene_id(gene_symbol)
    if gene_id is None:
        return []
    pathway_ids = fetch_kegg_pathway_ids(gene_id)
    pathway_names = fetch_kegg_pathway_names(pathway_ids)
    return [{"id": pid, "name": pathway_names.get(pid, pid)} for pid in pathway_ids]


def enrich_with_string_partners(ranked_candidates: list[dict], top_n: int = 10) -> list[dict]:
    """Adds STRING text-mining partner signal for the top_n highest-ranked
    candidates (not all of them - keeps the number of external calls
    reasonable)."""
    for entry in ranked_candidates[:top_n]:
        symbol = entry["symbols_seen"][0] if entry["symbols_seen"] else None
        if not symbol:
            entry["string_partners"] = []
            continue
        try:
            entry["string_partners"] = string_textmining_partners(symbol)
        except requests.exceptions.RequestException as e:
            entry["string_partners"] = []
            print(f"    STRING lookup failed for {symbol} ({e.__class__.__name__})")
        _polite_pause()
    return ranked_candidates


def enrich_with_kegg_pathways(ranked_candidates: list[dict], top_n: int = 10) -> list[dict]:
    """Adds KEGG pathways for the top_n highest-ranked candidates."""
    for entry in ranked_candidates[:top_n]:
        symbol = entry["symbols_seen"][0] if entry["symbols_seen"] else None
        if not symbol:
            entry["kegg_pathways"] = []
            continue
        try:
            entry["kegg_pathways"] = kegg_pathways_for_gene(symbol)
        except requests.exceptions.RequestException as e:
            entry["kegg_pathways"] = []
            print(f"    KEGG lookup failed for {symbol} ({e.__class__.__name__})")
        _polite_pause()
    return ranked_candidates


# --- Display -----------------------------------------------------------------

def print_candidate_table(ranked_candidates: list[dict], category: str) -> None:
    """Prints a ranked candidate table for manual review. Flags genes that
    are already in CANDIDATE_GENE_SETS[category]."""
    known = set(s.lower() for s in CANDIDATE_GENE_SETS.get(category, []))

    print(f"\n== Text-mining candidates: {category} (zebrafish literature) ==\n")
    header = (
        f"{'Symbol(s)':<25} | {'PMIDs':<6} | {'Mentions':<9} | {'Already listed?':<16} "
        "| STRING text-mining partners | KEGG pathways"
    )
    print(header)
    print("-" * len(header))

    for entry in ranked_candidates:
        symbols = "/".join(entry["symbols_seen"])
        already = "yes" if known & set(s.lower() for s in entry["symbols_seen"]) else ""
        partners = entry.get("string_partners")
        partner_str = ", ".join(p["partner_symbol"] for p in partners[:3]) if partners else "-"
        pathways = entry.get("kegg_pathways")
        pathway_str = ", ".join(p["name"] for p in pathways[:2]) if pathways else "-"
        print(
            f"{symbols:<25} | {entry['pmid_count']:<6} | {entry['mention_count']:<9} "
            f"| {already:<16} | {partner_str:<27} | {pathway_str}"
        )

    print_caveats()


def print_caveats() -> None:
    print(
        "\nNOTE - read before pasting anything into CANDIDATE_GENE_SETS:\n"
        "  1) Co-occurrence in literature = correlation, not causation. A high PMID rank\n"
        "     means 'often named together with the search term', not 'proven cause of the phenotype'.\n"
        "  2) Symbol normalization: PubTator returns the text form the article used, which\n"
        "     may differ from ZFIN's lowercase zebrafish convention (e.g. 'mstnb' vs 'MSTNB' vs\n"
        "     'myostatin b'). Verify the correct canonical symbol and species ortholog via UniProt\n"
        "     (candidate_gene_dossier.uniprot_gene_hits) before writing it in.\n"
        "  3) STRING text-mining score (tscore) and KEGG pathway membership are clues about\n"
        "     further candidates/mechanisms to look at, not confirmed findings for pikeperch.\n"
        "  4) This script does NOT write to CANDIDATE_GENE_SETS. Copy the symbols you judge\n"
        "     reasonable into candidate_gene_dossier.py by hand, then run that\n"
        "     dossier (single/batch mode) to verify against the pikeperch species."
    )


def main():
    categories = list(PHENOTYPE_QUERY_TERMS)
    category = clean_input(
        f"Phenotype category ({'/'.join(categories)}) [Enter = growth]: "
    ).lower()
    if not category:
        category = "growth"
    if category not in PHENOTYPE_QUERY_TERMS:
        print(f"Unknown category '{category}', using 'growth'.")
        category = "growth"

    extra = clean_input("Extra search term (Enter for none): ")
    extra_terms = [extra] if extra else None

    ranked = discover_candidates_for_category(category, extra_terms)

    if not ranked:
        print("No gene candidates found - try a different category or extra search term.")
        return

    enrich = clean_input("Enrich top candidates with STRING + KEGG? [y/N]: ").lower()
    if enrich == "y":
        ranked = enrich_with_string_partners(ranked)
        ranked = enrich_with_kegg_pathways(ranked)

    print_candidate_table(ranked, category)


if __name__ == "__main__":
    main()
