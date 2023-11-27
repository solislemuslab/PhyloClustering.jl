using Documenter, PhyloClustering

makedocs(
    sitename="PhyloClustering.jl",
    authors="Yibo Kong, Claudia Solís-Lemus, and contributors",
    modules=[PhyloClustering],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/installation.md",
            "Implementation" => "man/implementation.md",
        ],
        "Library" => [
            "Public Methods" => "lib/public_methods.md",
        ]
    ]
)

deploydocs(
    # repo = "github.com/solislemuslab/PhyloDiamond.jl.git",
)