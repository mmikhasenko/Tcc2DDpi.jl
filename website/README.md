# Refresh the documentation site

Use Julia 1.12, Node.js with npm, Python 3, and the sibling
`xDDPhaseSpace.jl` checkout at commit
`3590e02cda1d2366b9d2e32fc9342fde4ef7d44c` (the package CI pin).
The notebooks share `notebooks/Project.toml` and its lockfile. The
exporter and documentation tools use the separate `website/` environment.

From the repository root, set up PyCall's Python and instantiate both
Julia environments once:

```sh
python3 -m venv /tmp/tcc-notebook-python
/tmp/tcc-notebook-python/bin/python -m pip install -r notebooks/requirements.txt
PYTHON=/tmp/tcc-notebook-python/bin/python julia +1.12 --project=notebooks -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
julia +1.12 --project=website -e 'using Pkg; Pkg.instantiate()'
julia +1.12 --project=website website/build.jl N-100-resolution-animation.jl
python3 website/check.py public-site --complete
mkdir -p /tmp/tcc-pages-preview
ln -s "$PWD/public-site" /tmp/tcc-pages-preview/Tcc2DDpi.jl
python3 -m http.server --directory /tmp/tcc-pages-preview 8000
```

Open `http://localhost:8000/Tcc2DDpi.jl/` for a preview.

Pass one or more notebook filenames to refresh those snapshots and retain the
others. Pass `--all` for a full rebuild before replacing the old site, or omit arguments to rebuild just the
documentation around the existing snapshots. Failed notebook cells stop the
build before the published directory is replaced. Review the exports, then
commit `public-site/` with the notebook changes. Pages CI validates and deploys
that directory. Notebook execution and npm are local build steps.

The notebook HTML preserves Pluto's own renderer and adds a small stylesheet
and a link back to the DocumenterVitepress examples index. Interactive controls
that need Julia require running the source notebook in Pluto.
