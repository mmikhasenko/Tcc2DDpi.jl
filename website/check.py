"""Check the published entry points and local notebook links."""
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = set()
        self.assets = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag in {"link", "script"}:
            asset = attrs.get("href") if tag == "link" else attrs.get("src")
            if asset and asset.startswith("/Tcc2DDpi.jl/"):
                self.assets.add(asset)
        if tag == "a":
            href = attrs.get("href")
            if href:
                self.hrefs.add(href)


site = Path(sys.argv[1])
for entry in ("index.html", "examples.html", "notebooks/notebook.css", ".nojekyll"):
    if not (site / entry).is_file():
        sys.exit(f"Missing published file: {entry}")

for page in ("index.html", "examples.html"):
    html = Links()
    html.feed((site / page).read_text(encoding="utf-8"))
    if not html.assets:
        sys.exit(f"{page} has no assets under the GitHub Pages base")
    for asset in html.assets:
        target = unquote(urlsplit(asset).path).removeprefix("/Tcc2DDpi.jl/")
        if not (site / target).is_file():
            sys.exit(f"{page} links to missing asset: {asset}")

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
if "--complete" in sys.argv[2:]:
    sources = Path(__file__).resolve().parents[1] / "notebooks"
    required = set()
    for source in sources.glob("*.jl"):
        with source.open(encoding="utf-8") as notebook:
            if notebook.readline().strip() == "### A Pluto.jl notebook ###":
                required.add(f"notebooks/{source.stem}.html")
    if required != expected:
        sys.exit(f"Notebook coverage mismatch: missing exports {required - expected}, obsolete exports {expected - required}")
print(f"Verified {len(expected)} notebook links")
