"""
Vercel serverless function exposing the pikeperch in-silico breeding pipeline as
a web API for the interactive demo (see index.html): candidate_gene_dossier
.gather_dossier (Part 1) and ocs.run_breeding_pipeline (Part 2 - marker/phenotype
association -> breeding value -> GRM -> OCS pairing -> F-projection).

/api/gene_dossier and /api/breeding are both served from this single
handler class, not two separate files, because this project's Vercel
Python build (pyproject.toml-driven, declares its one entrypoint via
[tool.vercel]) only supports a single Python entrypoint per deployment -
same constraint documented in the kmer-assembler repo this project was
split out from (see project memory, project_gos_avelspipeline.md). The two
routes stay independent HTTP requests either way (the browser fires
separate fetch calls, so a slow breeding-pipeline run can't block a gene
lookup in flight).

Standalone project - split out from the "kmer-assembler" repo 2026-09-02
because that repo's Vercel deployment forced this gene-dossier lookup and
the unrelated Trinity-like assembly pipeline demo onto the same page.
"""
import json
import sys
from http.server import BaseHTTPRequestHandler
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from candidate_gene_dossier import BRIDGE_SPECIES, gather_dossier
from marker_phenotype_association import parse_phenotypes, parse_genotypes
from ocs import NUM_CANDIDATES, GENERATIONS, RELATEDNESS_THRESHOLD, run_breeding_pipeline


class handler(BaseHTTPRequestHandler):
    def do_GET(self):
        index_path = Path(__file__).resolve().parent.parent / "index.html"
        body = index_path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        if self.path.startswith("/api/breeding"):
            self._handle_breeding()
        else:
            self._handle_gene_dossier()

    def _handle_gene_dossier(self):
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length)

        try:
            payload = json.loads(body or b"{}")
            gene_symbol = str(payload.get("gene_symbol", "")).strip()
            if not gene_symbol:
                raise ValueError("gene_symbol is required")

            species_list = [
                str(s).strip() for s in payload.get("species", []) if str(s).strip()
            ]
            if not species_list:
                # Kept short (1 species + bridge species = 2) by design:
                # each extra species multiplies the number of network
                # calls, and the whole request must fit within Vercel's
                # 30s maxDuration (see vercel.json).
                species_list = ["Sander lucioperca"]
            if BRIDGE_SPECIES not in species_list:
                species_list.append(BRIDGE_SPECIES)

            # Shorter than the CLI's default 60s for the same reason - see
            # gather_dossier/ensembl_paralog_count docstrings in
            # candidate_gene_dossier.py. gather_dossier runs everything in
            # parallel across threads, so a 20s timeout doesn't raise the
            # total response time linearly - an earlier sequential 8s
            # variant timed out repeatedly in production (a different
            # network path to Ensembl than tested locally), see history in
            # candidate_gene_dossier.py.
            result = gather_dossier(gene_symbol, species_list, paralog_timeout=20)
            status = 200
        except Exception as exc:
            result = {"error": str(exc)}
            status = 400

        self._send_json(result, status)

    def _handle_breeding(self):
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length)

        try:
            payload = json.loads(body or b"{}")
            genotype_csv = payload.get("genotype_csv")
            phenotype_csv = payload.get("phenotype_csv")
            if not genotype_csv or not phenotype_csv:
                raise ValueError("genotype_csv and phenotype_csv (CSV text) are required")

            phenotype_column = str(payload.get("phenotype_column") or "fcr").strip()
            num_candidates = int(payload.get("num_candidates", NUM_CANDIDATES))
            relatedness_threshold = float(payload.get("relatedness_threshold", RELATEDNESS_THRESHOLD))
            generations = int(payload.get("generations", GENERATIONS))

            genotypes = parse_genotypes(genotype_csv)
            phenotypes = parse_phenotypes(phenotype_csv, phenotype_column)
            if not genotypes:
                raise ValueError("genotype_csv contained no rows")
            if len(genotypes) < 2:
                raise ValueError("need at least 2 animals with genotype to form pairs")
            num_candidates = min(num_candidates, len(genotypes))

            result = run_breeding_pipeline(
                genotypes,
                phenotypes,
                num_candidates=num_candidates,
                relatedness_threshold=relatedness_threshold,
                generations=generations,
            )
            status = 200
        except Exception as exc:
            result = {"error": str(exc)}
            status = 400

        self._send_json(result, status)

    def _send_json(self, result: dict, status: int) -> None:
        response = json.dumps(result).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(response)))
        self.end_headers()
        self.wfile.write(response)
