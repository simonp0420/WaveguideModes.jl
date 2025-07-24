using Pkg: Pkg
Pkg.activate(@__DIR__)
using Literate, Documenter, DocumenterCitations, WaveguideModes, DocumenterMermaid, Unitful
olddir = pwd()
cd(@__DIR__)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

literate_list = ["Contents", "index", "Tutorial", "autodocs"]
for file in literate_list
    Literate.markdown(joinpath("literate", file*".jl"), "src")
end

#Literate.notebook(joinpath("literate", "Tutorial.jl"), "notebooks"; execute=false)

makedocs(;
    modules = [WaveguideModes],
    repo = "github.com/simonp0420/WaveguideModes.jl.git",
    format = Documenter.HTML(
        assets=String["assets/citations.css"],
    ),
    pages = [
        "Contents" => "Contents.md",
        "Introduction" => "index.md",
        "Tutorial" => "Tutorial.md",
        "API Reference" => "autodocs.md",
        "References" => "references.md",
    ],
    sitename = "WaveguideModes",
    plugins = [bib,],
)

deploydocs(
    repo = "github.com/simonp0420/WaveguideModes.jl.git",
)

cd(olddir)