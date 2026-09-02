import sys
from http.server import ThreadingHTTPServer
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "api"))
from gene_dossier import handler

# ThreadingHTTPServer (not plain HTTPServer): a real browser opens more than
# one connection per page (keep-alive / speculative preconnects). See the
# kmer-assembler repo's _dev_server.py, where this same fix is documented -
# this file was copied from there when the gos project was split out.
#
# `vercel dev` doesn't work on Windows checkouts of the parent repo (a
# Vercel CLI path-escaping bug, see kmer-assembler's _dev_server.py) -
# running the real `handler` class directly sidesteps that, same code
# Vercel would run, just served locally.
server = ThreadingHTTPServer(("127.0.0.1", 3112), handler)
print("Serving gos-avel on http://127.0.0.1:3112")
server.serve_forever()
