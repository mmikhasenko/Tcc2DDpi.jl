using Documenter, DocumenterVitepress, TOML

const root = dirname(@__DIR__)
const titles = TOML.parsefile(joinpath(@__DIR__, "notebook-titles.toml"))
const published = joinpath(root, "public-site")
const snapshots = joinpath(@__DIR__, "snapshots")

mktempdir() do tempdir
    tmp = realpath(tempdir)
    source = joinpath(tmp, "src")
    cp(joinpath(@__DIR__, "src"), source)
    exports = joinpath(tmp, "notebooks")
    cp(snapshots, exports)
    cp(joinpath(@__DIR__, "notebook.css"), joinpath(exports, "notebook.css"))
    notebooks = sort(filter(readdir(joinpath(root, "notebooks"))) do file
        endswith(file, ".jl") && open(readline, joinpath(root, "notebooks", file)) == "### A Pluto.jl notebook ###"
    end)
    open(joinpath(source, "examples.md"), "w") do io
        println(io, "# Notebook examples\n\nRead the saved outputs or open the Julia source in Pluto.\n")
        println(io, "```@raw html\n<ul>")
        for file in notebooks
            name = splitext(file)[1]
            isfile(joinpath(exports, name * ".html")) || continue
            title = get(titles, name, name)
            title = replace(title, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
            println(io, """<li><a href="./notebooks/$name.html">$title</a> · <a href="https://github.com/mmikhasenko/Tcc2DDpi.jl/blob/master/notebooks/$file">Julia source</a></li>""")
        end
        println(io, "</ul>\n```")
    end

    withenv("CI" => "true") do
        makedocs(; root=tmp, remotes=nothing, sitename="X2DDpi", pages=["Overview" => "index.md", "Examples" => "examples.md"],
            format=MarkdownVitepress(repo="github.com/mmikhasenko/Tcc2DDpi.jl", devbranch="master", write_inventory=false,
                deploy_url="https://mmikhasenko.github.io/Tcc2DDpi.jl/",
                deploy_decision=Documenter.DeployDecision(all_ok=false)))
    end
    site = joinpath(tmp, "build", "1")
    isfile(joinpath(site, "examples.html")) || error("Missing built examples page")
    cp(exports, joinpath(site, "notebooks"))
    write(joinpath(site, ".nojekyll"), "")
    # DocumenterVitepress references this even for a single unversioned site.
    write(joinpath(site, "versions.js"), "var DOC_VERSIONS = [];\n")
    run(`python3 $(joinpath(@__DIR__, "check.py")) $site --complete`)
    cp(site, published; force=true)
end
