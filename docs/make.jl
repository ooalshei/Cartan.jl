using Documenter, RedCarD

DocMeta.setdocmeta!(RedCarD, :DocTestSetup, :(using RedCarD); recursive=true)

makedocs(
    sitename="RedCarD.jl",
    modules=[RedCarD],
    authors="Omar Alsheikh",
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Manual" => [
            "Algebras and decompositions" => "decompositions.md",
            "Finding the circuit" => "optimization.md",
        ],
        "API reference" => "api.md",
    ],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://ooalshei.github.io/RedCarD.jl",
        repolink="https://github.com/ooalshei/RedCarD.jl",
    ),
)

deploydocs(
    repo = "https://github.com/ooalshei/RedCarD.jl",
    devbranch = "main"
)
