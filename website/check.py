"""Check the published entry points and local notebook links."""
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = set()

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.add(href)


site = Path(sys.argv[1])
for entry in ("index.html", "examples.html", "notebooks/notebook.css", ".nojekyll"):
    if not (site / entry).is_file():
        sys.exit(f"Missing published file: {entry}")

index = Links()
index.feed((site / "examples.html").read_text(encoding="utf-8"))
targets = set()
for href in index.hrefs:
    url = urlsplit(href)
    if url.scheme or url.netloc:
        continue
    path = unquote(url.path).removeprefix("/Tcc2DDpi.jl/").removeprefix("./")
    if path.startswith("notebooks/") and path.endswith(".html"):
        targets.add(path)
expected = {p.relative_to(site).as_posix() for p in (site / "notebooks").glob("*.html")}
if not expected or targets != expected:
    sys.exit(f"Notebook index mismatch: missing files {targets - expected}, missing links {expected - targets}")
print(f"Verified {len(expected)} notebook links")
