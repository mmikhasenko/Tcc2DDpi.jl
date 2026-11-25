using Documenter, DocumenterVitepress, Pluto

const root = dirname(@__DIR__)
const published = joinpath(root, "public-site")
ENV["GKSwstype"] = "100"

function export_notebook(session, path, output)
    notebook = Pluto.SessionActions.open(session, path; run_async=false)
    try
        failed = [cell.cell_id for cell in notebook.cells if cell.errored]
        isempty(failed) || error("$(basename(path)) has failed cells: $failed")
        html = Pluto.generate_html(notebook)
        nav = """<nav class="notebook-nav"><a href="../examples.html">← Notebook examples</a></nav>"""
        html = replace(html, "</head>" => """<link rel="stylesheet" href="notebook.css"></head>"""; count=1)
        write(output, replace(html, r"<body[^>]*>" => SubstitutionString("\\0" * nav); count=1))
    finally
        Pluto.SessionActions.shutdown(session, notebook)
    end
end

mktempdir() do tmp
    source = joinpath(tmp, "src")
    cp(joinpath(@__DIR__, "src"), source)
    exports = mkpath(joinpath(tmp, "notebooks"))
    if isdir(joinpath(published, "notebooks")) && ARGS != ["--all"]
        cp(joinpath(published, "notebooks"), exports; force=true)
    end
    cp(joinpath(@__DIR__, "notebook.css"), joinpath(exports, "notebook.css"); force=true)

    notebooks = sort(filter(readdir(joinpath(root, "notebooks"))) do file
        endswith(file, ".jl") && Pluto.is_pluto_notebook(joinpath(root, "notebooks", file))
    end)
    selected = ARGS == ["--all"] ? notebooks : ARGS
    all(in(notebooks), selected) || error("Pass notebook filenames from notebooks/ or --all")
    session = Pluto.ServerSession(; options=Pluto.Configuration.from_flat_kwargs(
        disable_writing_notebook_files=true))
    for file in selected
        export_notebook(session, joinpath(root, "notebooks", file),
            joinpath(exports, splitext(file)[1] * ".html"))
    end

    open(joinpath(source, "examples.md"), "w") do io
        println(io, "# Notebook examples\n\nRead the saved outputs or open the Julia source in Pluto.\n")
        for file in notebooks
            name = splitext(file)[1]
            isfile(joinpath(exports, name * ".html")) || continue
            println(io, "- <a href=\"./notebooks/$name.html\">$name</a> · [Julia source](https://github.com/mmikhasenko/Tcc2DDpi.jl/blob/master/notebooks/$file)")
        end
    end

    withenv("CI" => "true") do
        makedocs(; root=tmp, sitename="X2DDpi", pages=["Overview" => "index.md", "Examples" => "examples.md"],
            format=MarkdownVitepress(repo="github.com/mmikhasenko/Tcc2DDpi.jl", devbranch="master", write_inventory=false,
                deploy_url="https://mmikhasenko.github.io/Tcc2DDpi.jl/"))
    end
    site = joinpath(tmp, "build", "1")
    isfile(joinpath(site, "examples.html")) || error("Missing built examples page")
    cp(exports, joinpath(site, "notebooks"))
    write(joinpath(site, ".nojekyll"), "")
    run(`python3 $(joinpath(@__DIR__, "check.py")) $site`)
    cp(site, published; force=true)
end
