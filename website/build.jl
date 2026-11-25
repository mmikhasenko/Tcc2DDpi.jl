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

mktempdir() do tempdir
    tmp = realpath(tempdir)
    source = joinpath(tmp, "src")
    cp(joinpath(@__DIR__, "src"), source)
    exports = mkpath(joinpath(tmp, "notebooks"))
    if isdir(joinpath(published, "notebooks")) && ARGS != ["--all"]
        cp(joinpath(published, "notebooks"), exports; force=true)
    end

    notebooks = sort(filter(readdir(joinpath(root, "notebooks"))) do file
        endswith(file, ".jl") && Pluto.is_pluto_notebook(joinpath(root, "notebooks", file))
    end)
    expected = Set(splitext(file)[1] * ".html" for file in notebooks)
    for file in readdir(exports)
        endswith(file, ".html") && !(file in expected) && rm(joinpath(exports, file))
    end
    cp(joinpath(@__DIR__, "notebook.css"), joinpath(exports, "notebook.css"); force=true)
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
        println(io, "```@raw html\n<ul>")
        for file in notebooks
            name = splitext(file)[1]
            isfile(joinpath(exports, name * ".html")) || continue
            heading = match(r"md\x22\x22\x22\r?\n# ([^\r\n]+)", read(joinpath(root, "notebooks", file), String))
            title = isnothing(heading) ? name : heading.captures[1]
            title = replace(title, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
            println(io, """<li><a href="./notebooks/$name.html">$title</a> · <a href="https://github.com/mmikhasenko/Tcc2DDpi.jl/blob/master/notebooks/$file">Julia source</a></li>""")
        end
        println(io, "</ul>\n```")
    end

    withenv("CI" => "true") do
        makedocs(; root=tmp, remotes=nothing, sitename="X2DDpi", pages=["Overview" => "index.md", "Examples" => "examples.md"],
            format=MarkdownVitepress(repo="github.com/mmikhasenko/Tcc2DDpi.jl", devbranch="master", write_inventory=false,
                deploy_url="https://mmikhasenko.github.io/Tcc2DDpi.jl/"))
    end
    site = joinpath(tmp, "build", "1")
    isfile(joinpath(site, "examples.html")) || error("Missing built examples page")
    cp(exports, joinpath(site, "notebooks"))
    write(joinpath(site, ".nojekyll"), "")
    # DocumenterVitepress references this even for a single unversioned site.
    write(joinpath(site, "versions.js"), "var DOC_VERSIONS = [];\n")
    run(`python3 $(joinpath(@__DIR__, "check.py")) $site`)
    cp(site, published; force=true)
end
