# Refresh the documentation site

Use Julia 1.12, Node.js with npm, Python 3, and the package's normal sibling
`xDDPhaseSpace.jl` checkout. Notebook environments can have additional
requirements such as Python for SymPy.

From the repository root:

```sh
julia +1.12 --project=website -e 'using Pkg; Pkg.instantiate()'
julia +1.12 --project=website website/build.jl N-100-resolution-animation.jl
python3 website/check.py public-site
mkdir -p /tmp/tcc-pages-preview
ln -s "$PWD/public-site" /tmp/tcc-pages-preview/Tcc2DDpi.jl
python3 -m http.server --directory /tmp/tcc-pages-preview 8000
```

Open `http://localhost:8000/Tcc2DDpi.jl/` for a preview.

Pass one or more notebook filenames to refresh those snapshots and retain the
others. Pass `--all` for a full rebuild, or omit arguments to rebuild just the
documentation around the existing snapshots. Failed notebook cells stop the
build before the published directory is replaced. Review the exports, then
commit `public-site/` with the notebook changes. Pages CI validates and deploys
that directory. Notebook execution and npm are local build steps.

The notebook HTML preserves Pluto's own renderer and adds a small stylesheet
and a link back to the DocumenterVitepress examples index. Interactive controls
that need Julia require running the source notebook in Pluto.
