using Pluto

const root = dirname(@__DIR__)
const snapshots = joinpath(@__DIR__, "snapshots")
ENV["GKSwstype"] = "100"
isempty(ARGS) && error("Pass notebook filenames or --all")

function export_notebook(session, path, output)
    notebook = Pluto.SessionActions.open(session, path; run_async=false)
    try
        failed = filter(cell -> cell.errored, notebook.cells)
        for cell in failed
            @error "Notebook cell failed" notebook=basename(path) cell=cell.cell_id output=cell.output.body
        end
        isempty(failed) || error("$(basename(path)) has $(length(failed)) failed cells")
        html = Pluto.generate_html(notebook)
        nav = """<nav class="notebook-nav"><a href="../examples.html">← Notebook examples</a></nav>"""
        html = replace(html, "</head>" => """<link rel="stylesheet" href="notebook.css"></head>"""; count=1)
        write(output, replace(html, r"<body[^>]*>" => SubstitutionString("\\0" * nav); count=1))
    finally
        Pluto.SessionActions.shutdown(session, notebook)
    end
end

mktempdir() do tmp
    exports = joinpath(tmp, "snapshots")
    if ARGS == ["--all"]
        mkdir(exports)
    else
        cp(snapshots, exports)
    end
    notebooks = sort(filter(readdir(joinpath(root, "notebooks"))) do file
        endswith(file, ".jl") && open(readline, joinpath(root, "notebooks", file)) == "### A Pluto.jl notebook ###"
    end)
    expected = Set(splitext(file)[1] * ".html" for file in notebooks)
    for file in readdir(exports)
        endswith(file, ".html") && !(file in expected) && rm(joinpath(exports, file))
    end
    selected = ARGS == ["--all"] ? notebooks : ARGS
    all(in(notebooks), selected) || error("Pass notebook filenames from notebooks/ or --all")
    session = Pluto.ServerSession(; options=Pluto.Configuration.from_flat_kwargs(
        disable_writing_notebook_files=true))
    failures = String[]
    for file in selected
        @info "Exporting notebook" file
        try
            elapsed = @elapsed withenv("TCC_EXPORT_OUTPUTS" => tmp) do
                export_notebook(session, joinpath(root, "notebooks", file),
                    joinpath(exports, splitext(file)[1] * ".html"))
            end
            @info "Exported notebook" file seconds=round(elapsed; digits=1)
        catch err
            err isa InterruptException && rethrow()
            @error "Notebook export failed" file exception=(err, catch_backtrace())
            push!(failures, file)
        end
    end
    isempty(failures) || error("Notebook exports failed: $(join(failures, ", "))")

    cp(exports, snapshots; force=true)
end
