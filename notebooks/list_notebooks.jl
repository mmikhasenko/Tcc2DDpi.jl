using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
# 
using Pluto

# make sure the script is located in the notebooks folder
# the script assumes it is a subfolder of the project
notebooks_dir = "notebooks"
@assert splitpath(@__DIR__)[end] == notebooks_dir

files = readdir(@__DIR__)
notebooks = filter(files) do f
    Pluto.is_pluto_notebook(joinpath(@__DIR__, f))
end

function extract_md_blocks(content::AbstractString)
    r = r"md\"\"\"\r?\n(.*?)\"\"\""s
    return [m.captures[1] for m in eachmatch(r, content)]
end

function parse_heading(heading::AbstractString)
    firstheading = split(heading, '\n') |> first
    m = match(r"^(#+) (.*)\r?$", firstheading)
    level = length(m.captures[1])
    title = strip(m.captures[2])
    return level, title
end

function build_index_structure(content::AbstractString)
    # Extract all markdown blocks
    blocks = extract_md_blocks(content)
    filter!(blocks) do block
        occursin(r"^(#+) \r?.*$"s, block)
    end
    # Extract headings from blocks and their levels
    headings = parse_heading.(blocks)

    # Initialize root of the structure
    root = Pair("root", [])
    current_struct = [root]

    for (level, title) in headings
        while length(current_struct) < level
            # Fill missing levels with "Untitled"
            new_untitled = Pair("Untitled", [])
            push!(last(current_struct).second, new_untitled)
            push!(current_struct, new_untitled)
        end
        while length(current_struct) > level
            pop!(current_struct)
        end

        # Add current title to the structure
        new_entry = Pair(title, [])
        push!(last(current_struct).second, new_entry)

        if level < length(current_struct)
            current_struct[level] = new_entry
        else
            push!(current_struct, new_entry)
        end
    end

    return root.second
end



foldername = last(splitpath(pwd()))
# 
html_string = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(foldername)</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <div class="container">
        <!-- Index Panel on the left -->
        <div class="index-panel">
            <h3>Index</h3>
            <ul>
                <li><a href="README.html" target="notebook-frame">README</a></li>
            </ul>
        </div>
        <!-- Notebook Display on the right -->
        <div class="iframe-container">
            <iframe name="notebook-frame" class="notebook-display" src="README.html"></iframe>
        </div>
    </div>
</body>
</html>
"""


function extend_itemized_list(notebooks, html_string)
    index_end = first(findfirst("</ul>", html_string)) - 1
    # 
    list_items = ""
    for notebook in notebooks
        notebook_html = replace(notebook, ".jl" => ".html")
        index = build_index_structure(read(joinpath(@__DIR__, notebook), String))

        # Extract titles and subtitles from the index
        length(index) != 1 && error("Fishy index, not `Title=>[Subs...]`\n\nindex = $(index)\n")
        title = first(index[1])
        subtitles = last(index[1])

        subtitle_list_items = ""
        for subtitle in subtitles
            subtitle_name = first(subtitle)
            subtitle_list_items *= """<li>$(subtitle_name)</li>"""
        end

        # incorporate to html
        list_items *= """
        <li>
            <a href="notebooks/$(notebook_html)" target="notebook-frame">$(notebook)</a>
            <div class="notebook-details">
                <span class="title">$(title)</span>
                <ul class="subtitles">$(subtitle_list_items)</ul>
            </div>
        </li>
        """
    end

    extended_html_string = string(html_string[1:index_end-1], list_items, html_string[index_end:end])

    return extended_html_string
end

extended_html_string = extend_itemized_list(notebooks, html_string)
write(joinpath(@__DIR__, "..", "index.html"), extended_html_string)
