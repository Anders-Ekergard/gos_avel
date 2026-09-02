"""
Vercel serverless function exposing candidate_gene_dossier.gather_dossier as
a web API for the interactive demo (see index.html).

Standalone project - split out from the "kmer-assembler" repo 2026-09-02
because that repo's Vercel deployment only supports one Python entrypoint,
which forced this gene-dossier lookup and the unrelated Trinity-like
assembly pipeline demo onto the same page. See project memory
(project_gos_avelspipeline.md) for why.
"""
import json
import sys
from http.server import BaseHTTPRequestHandler
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from candidate_gene_dossier import BRIDGE_SPECIES, gather_dossier


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

        response = json.dumps(result).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(response)))
        self.end_headers()
        self.wfile.write(response)
