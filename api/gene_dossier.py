"""
Vercel serverless function exposing the gös in-silico breeding pipeline as
a web API for the interactive demo (see index.html): candidate_gene_dossier
.gather_dossier (Del 1) and ocs.kor_avelspipeline (Del 2 - marker/phenotype
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
from marker_fenotyp_koppling import parse_fenotyper, parse_genotyper
from ocs import ANTAL_KANDIDATER, GENERATIONER, SLAKTSKAP_TROSKEL, kor_avelspipeline


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
                # Kept short (1 art + bryggart = 2 st) med avsikt: varje extra
                # art multiplicerar antalet nätverksanrop, och hela requesten
                # måste rymmas inom Vercels 30s maxDuration (se vercel.json).
                species_list = ["Sander lucioperca"]
            if BRIDGE_SPECIES not in species_list:
                species_list.append(BRIDGE_SPECIES)

            # Kortare an CLI:ts standard-60s av samma skal - se
            # gather_dossier/ensembl_paralog_count docstrings i
            # candidate_gene_dossier.py. gather_dossier kor allt parallellt
            # i trådar, sa 20s timeout hojer inte totala svarstiden linjart -
            # en tidigare sekventiell 8s-variant timeoutade upprepat i
            # produktion (annan natverksväg till Ensembl an lokalt testat),
            # se historik i candidate_gene_dossier.py.
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
            genotyp_csv = payload.get("genotyp_csv")
            fenotyp_csv = payload.get("fenotyp_csv")
            if not genotyp_csv or not fenotyp_csv:
                raise ValueError("genotyp_csv and fenotyp_csv (CSV text) are required")

            fenotyp_kolumn = str(payload.get("fenotyp_kolumn") or "fcr").strip()
            antal_kandidater = int(payload.get("antal_kandidater", ANTAL_KANDIDATER))
            slaktskap_troskel = float(payload.get("slaktskap_troskel", SLAKTSKAP_TROSKEL))
            generationer = int(payload.get("generationer", GENERATIONER))

            genotyper = parse_genotyper(genotyp_csv)
            fenotyper = parse_fenotyper(fenotyp_csv, fenotyp_kolumn)
            if not genotyper:
                raise ValueError("genotyp_csv innehöll inga rader")
            if len(genotyper) < 2:
                raise ValueError("behöver minst 2 djur med genotyp för att bilda par")
            antal_kandidater = min(antal_kandidater, len(genotyper))

            result = kor_avelspipeline(
                genotyper,
                fenotyper,
                antal_kandidater=antal_kandidater,
                slaktskap_troskel=slaktskap_troskel,
                generationer=generationer,
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
