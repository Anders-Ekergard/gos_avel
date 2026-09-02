"""
Text-mining-triage: soker PubMed/PubTator3 efter gener som samnamns med
en fenotyp hos zebrafisk (Danio rerio - den fiskart med rikast funktionell
litteratur), och kan valfritt utoka traffarna via STRINGs textmining-kanal
samt KEGG-pathways (samma tva kallor som redan anvands for humana gener i
pathway_mvp.py, har omskrivna for zebrafisk-organismkoder istallet for
hsa/Homo sapiens).

Detta ersatter INTE manuell granskning - traffarna ar en startlista att
verifiera (UniProt/Ensembl via candidate_gene_dossier.py) och klippa in
i CANDIDATE_GENE_SETS for hand, en kategori i taget.

Kor: python text_mining_triage.py
"""

import functools
import time

import requests

from candidate_gene_dossier import CANDIDATE_GENE_SETS, EUTILS_BASE, clean_input

PUBTATOR_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
STRING_BASE = "https://string-db.org/api"
KEGG_BASE = "https://rest.kegg.jp"

# Val av "brygg-art" for textmining: zebrafisk har overlagset mest
# funktionell litteratur for fisk (se motivering i candidate_gene_dossier.py).
DISCOVERY_SPECIES = "Danio rerio"
DISCOVERY_TAXON_ID = 7955       # NCBI taxon id, kravs av STRING
DISCOVERY_KEGG_ORGANISM = "dre"  # KEGG:s organismkod for Danio rerio

# Sokord per fenotypkategori - engelska, eftersom PubMed/PubTator ar
# engelskspraktig litteratur. Nycklarna maste matcha CANDIDATE_GENE_SETS
# i candidate_gene_dossier.py sa att traffar kan korsrefereras mot
# redan kanda kandidater.
PHENOTYPE_QUERY_TERMS: dict[str, list[str]] = {
    "tillvaxt": ["growth rate", "feed conversion ratio", "muscle growth", "body weight"],
    "kon": ["sex determination", "sex differentiation", "gonad development"],
    "kottkvalitet": ["muscle fiber", "fillet quality", "flesh quality"],
    "sjukdom_stress": ["hypoxia tolerance", "stress response", "disease resistance", "immune response"],
}

MAX_PMIDS_PER_TERM = 50       # tak per sokterm - triage, inte fullstandig litteraturgenomgang
BIOCJSON_CHUNK_SIZE = 50      # antal PMIDs per annotations-anrop
REQUEST_DELAY_SECONDS = 0.4   # artighetspaus mellan anrop (~3/sek, samma norm som NCBI E-utilities)


def _polite_pause() -> None:
    """Enkel artighetspaus mellan anrop till PubTator3/STRING/KEGG - alla
    tre ar mindre, publikt finansierade tjanster utan dokumenterad hard
    rate-limit for den har typen av anrop; battre att vara snall."""
    time.sleep(REQUEST_DELAY_SECONDS)


# --- PubTator3 (upptackt) ---------------------------------------------------

def pubtator_search_pmids(query: str, max_results: int = MAX_PMIDS_PER_TERM) -> list[str]:
    """Sok PubTator3 efter artiklar som matchar en fritextfraga
    (t.ex. '"feed conversion ratio" Danio rerio'). Returnerar PMIDs.

    OBS: PubTator3 har inte anropats fran det har repot tidigare - fältnamnen
    nedan (`results`, `pmid`) ar basta-gissning fran offentlig dokumentation.
    Verifiera empiriskt (print(r.json())) innan traffarna litas pa blint."""
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
    """Hamtar BioC-JSON med normaliserade gen/art-annoteringar for en
    lista PMIDs, i klumpar om BIOCJSON_CHUNK_SIZE."""
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
    """Plockar ut Gene-annoteringar fran BioC-dokument, filtrerat till
    dokument som ocksa har en Species-annotering som matchar taxon_id."""
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
    """Aggregerar genmentions till en rankad kandidatlista: flest unika
    PMIDs forst (bredare stod i litteraturen), sedan flest mentions totalt.

    OBS: samforekomst i litteratur ar korrelation, inte kausalitet - se
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


def build_queries_for_category(kategori: str, extra_terms: list[str] | None = None) -> list[str]:
    """En sokfraga per synonymterm (inte en enda kombinerad boolesk
    fraga) - enklare att fa ratt utan att veta PubTator3s exakta
    frag-syntax, och lattare att debugga term for term."""
    terms = list(PHENOTYPE_QUERY_TERMS.get(kategori, []))
    if extra_terms:
        terms += extra_terms
    return [f'"{term}" {DISCOVERY_SPECIES}' for term in terms]


def discover_candidates_for_category(kategori: str, extra_terms: list[str] | None = None) -> list[dict]:
    """Full discovery-pipeline for en fenotypkategori: sok -> hamta
    annoteringar -> filtrera -> rangordna."""
    queries = build_queries_for_category(kategori, extra_terms)

    all_pmids: list[str] = []
    for q in queries:
        try:
            pmids = pubtator_search_pmids(q)
            print(f"  '{q}': {len(pmids)} PMIDs")
            for p in pmids:
                if p not in all_pmids:
                    all_pmids.append(p)
        except requests.exceptions.RequestException as e:
            print(f"  '{q}': sokning misslyckades ({e.__class__.__name__})")

    documents = pubtator_fetch_annotations(all_pmids)
    mentions = extract_gene_mentions(documents)
    return rank_gene_candidates(mentions)


# --- STRING (utokning fran kand gen, textmining-kanalen) --------------------

def string_get_string_id(gene_symbol: str, species_taxon_id: int = DISCOVERY_TAXON_ID) -> str | None:
    """Slar upp STRINGs interna ID for en gensymbol - kravs av
    interaction_partners (accepterar inte fria symboler direkt utan
    risk for tvetydig matchning)."""
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
    """Funktionella partners till en kand gen, filtrerat till STRINGs
    textmining-kanal (tscore) - dvs partners med stod i samforekomst-
    litteratur, inte bara databas-/experimentbevis."""
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


# --- KEGG (pathways for kand gen) -------------------------------------------
# Samma tvasteg-metodik (Entrez -> gen-id, KEGG link/list) som redan anvands
# for humana gener i pathway_mvp.py - har med "Danio rerio[orgn]"/"dre"
# istallet for "Homo sapiens[orgn]"/"hsa".

@functools.lru_cache(maxsize=64)
def resolve_kegg_gene_id(symbol: str, organism: str = DISCOVERY_SPECIES) -> str | None:
    """Slar upp ett KEGG-gen-id (t.ex. 'dre:napp') via NCBI Entrez esearch,
    samma eutils-tjanst som redan anvands i candidate_gene_dossier.py."""
    params = {"db": "gene", "term": f"{symbol}[sym] AND {organism}[orgn]", "retmode": "json"}
    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return f"{DISCOVERY_KEGG_ORGANISM}:{ids[0]}" if ids else None


def parse_kegg_tsv(text: str) -> list[tuple[str, str]]:
    """Tolkar KEGGs platta, tabbseparerade svarsformat ('nyckel\\tvarde'
    per rad), delat av bade link- och list-operationerna."""
    pairs = []
    for line in text.strip().splitlines():
        if not line:
            continue
        key, value = line.split("\t", 1)
        pairs.append((key, value))
    return pairs


def fetch_kegg_pathway_ids(gene_id: str) -> list[str]:
    """Vilka KEGG-pathways en gen tillhor."""
    r = requests.get(f"{KEGG_BASE}/link/pathway/{gene_id}", timeout=30)
    r.raise_for_status()
    return [pathway_id for _gene, pathway_id in parse_kegg_tsv(r.text)]


@functools.lru_cache(maxsize=8)
def _fetch_kegg_pathway_list(organism: str) -> tuple[tuple[str, str], ...]:
    """Hamtar och cachar en organisms fullstandiga pathway-id -> namn-lista
    (samma svar oavsett vilken gen som utlost uppslaget)."""
    r = requests.get(f"{KEGG_BASE}/list/pathway/{organism}", timeout=30)
    r.raise_for_status()
    return tuple(parse_kegg_tsv(r.text))


def fetch_kegg_pathway_names(pathway_ids: list[str], organism: str = DISCOVERY_KEGG_ORGANISM) -> dict[str, str]:
    """Namn for KEGG-pathway-id:n (link-svaret ger bara id:n, list-svaret
    behovs for namnen)."""
    if not pathway_ids:
        return {}
    names_by_bare_id = dict(_fetch_kegg_pathway_list(organism))
    return {
        pid: names_by_bare_id[pid.removeprefix("path:")]
        for pid in pathway_ids
        if pid.removeprefix("path:") in names_by_bare_id
    }


def kegg_pathways_for_gene(gene_symbol: str) -> list[dict]:
    """Hognivafunktion: gensymbol -> lista KEGG-pathways ({'id', 'name'})
    for zebrafisk. Tom lista om Entrez inte hittar genen (t.ex. fel
    stavning/synonym som inte matchar Entrez-index) - kraschar inte."""
    gene_id = resolve_kegg_gene_id(gene_symbol)
    if gene_id is None:
        return []
    pathway_ids = fetch_kegg_pathway_ids(gene_id)
    pathway_names = fetch_kegg_pathway_names(pathway_ids)
    return [{"id": pid, "name": pathway_names.get(pid, pid)} for pid in pathway_ids]


def enrich_with_string_partners(ranked_candidates: list[dict], top_n: int = 10) -> list[dict]:
    """Lagger till STRING-textmining-partnersignal for de top_n hogst
    rankade kandidaterna (inte alla - hall antalet externa anrop rimligt)."""
    for entry in ranked_candidates[:top_n]:
        symbol = entry["symbols_seen"][0] if entry["symbols_seen"] else None
        if not symbol:
            entry["string_partners"] = []
            continue
        try:
            entry["string_partners"] = string_textmining_partners(symbol)
        except requests.exceptions.RequestException as e:
            entry["string_partners"] = []
            print(f"    STRING-uppslag misslyckades for {symbol} ({e.__class__.__name__})")
        _polite_pause()
    return ranked_candidates


def enrich_with_kegg_pathways(ranked_candidates: list[dict], top_n: int = 10) -> list[dict]:
    """Lagger till KEGG-pathways for de top_n hogst rankade kandidaterna."""
    for entry in ranked_candidates[:top_n]:
        symbol = entry["symbols_seen"][0] if entry["symbols_seen"] else None
        if not symbol:
            entry["kegg_pathways"] = []
            continue
        try:
            entry["kegg_pathways"] = kegg_pathways_for_gene(symbol)
        except requests.exceptions.RequestException as e:
            entry["kegg_pathways"] = []
            print(f"    KEGG-uppslag misslyckades for {symbol} ({e.__class__.__name__})")
        _polite_pause()
    return ranked_candidates


# --- Visning -----------------------------------------------------------------

def print_candidate_table(ranked_candidates: list[dict], kategori: str) -> None:
    """Skriver ut rankad kandidattabell for manuell granskning. Markerar
    gener som redan finns i CANDIDATE_GENE_SETS[kategori]."""
    kanda = set(s.lower() for s in CANDIDATE_GENE_SETS.get(kategori, []))

    print(f"\n== Text-mining-kandidater: {kategori} (zebrafisk-litteratur) ==\n")
    header = (
        f"{'Symbol(er)':<25} | {'PMIDs':<6} | {'Mentions':<9} | {'Redan i listan?':<16} "
        "| STRING textmining-partners | KEGG-pathways"
    )
    print(header)
    print("-" * len(header))

    for entry in ranked_candidates:
        symbols = "/".join(entry["symbols_seen"])
        redan = "ja" if kanda & set(s.lower() for s in entry["symbols_seen"]) else ""
        partners = entry.get("string_partners")
        partner_str = ", ".join(p["partner_symbol"] for p in partners[:3]) if partners else "-"
        pathways = entry.get("kegg_pathways")
        pathway_str = ", ".join(p["name"] for p in pathways[:2]) if pathways else "-"
        print(
            f"{symbols:<25} | {entry['pmid_count']:<6} | {entry['mention_count']:<9} "
            f"| {redan:<16} | {partner_str:<27} | {pathway_str}"
        )

    print_caveats()


def print_caveats() -> None:
    print(
        "\nOBS - las innan nagot klipps in i CANDIDATE_GENE_SETS:\n"
        "  1) Samforekomst i litteratur = korrelation, inte kausalitet. Hog PMID-rankning\n"
        "     betyder 'namns ofta tillsammans med sokordet', inte 'bevisad orsak till fenotypen'.\n"
        "  2) Symbol-normalisering: PubTator returnerar den textform artikeln raknat med, som\n"
        "     kan skilja fran ZFIN:s gemena zebrafisk-konvention (t.ex. 'mstnb' vs 'MSTNB' vs\n"
        "     'myostatin b'). Verifiera ratt kanonisk symbol och art-ortholog via UniProt\n"
        "     (candidate_gene_dossier.uniprot_gene_hits) innan den skrivs in.\n"
        "  3) STRING textmining-score (tscore) och KEGG-pathwaymedlemskap ar ledtradar om\n"
        "     ytterligare kandidater/mekanismer att titta pa, inte bekraftade fynd for gos.\n"
        "  4) Detta skript skriver INTE till CANDIDATE_GENE_SETS. Kopiera de symboler du\n"
        "     bedomer rimliga in i candidate_gene_dossier.py for hand, och kor sedan den\n"
        "     dossiern (enskild/batch-lage) for att verifiera mot gos-arterna."
    )


def main():
    kategorier = list(PHENOTYPE_QUERY_TERMS)
    kategori = clean_input(
        f"Fenotypkategori ({'/'.join(kategorier)}) [Enter = tillvaxt]: "
    ).lower()
    if not kategori:
        kategori = "tillvaxt"
    if kategori not in PHENOTYPE_QUERY_TERMS:
        print(f"Okand kategori '{kategori}', anvander 'tillvaxt'.")
        kategori = "tillvaxt"

    extra = clean_input("Extra sokterm (Enter for ingen): ")
    extra_terms = [extra] if extra else None

    ranked = discover_candidates_for_category(kategori, extra_terms)

    if not ranked:
        print("Inga genkandidater hittade - prova en annan kategori eller extra sokterm.")
        return

    berika = clean_input("Berika topp-kandidater med STRING + KEGG? [j/N]: ").lower()
    if berika == "j":
        ranked = enrich_with_string_partners(ranked)
        ranked = enrich_with_kegg_pathways(ranked)

    print_candidate_table(ranked, kategori)


if __name__ == "__main__":
    main()
