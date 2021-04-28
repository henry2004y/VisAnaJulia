using Documenter, SpaceAnalysis, PyPlot

makedocs(;
    modules=[SpaceAnalysis],
    authors="Hongyang Zhou <hyzhou@umich.edu> and contributors",
    sitename="SpaceAnalysis",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/VisAnaJulia",
        assets=String[],
    ),
    pages = ["Home"          => "index.md",
             "Analysis"      => "man/analysis.md",
             "Example"       => "man/examples.md",
             "Log"           => "man/log.md",
             "Parker Spiral" => "man/parker_spiral.md",
             "Internal"      => "man/internal.md"
    ]
)

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia.git",
    target = "build",
    branch = "gh-pages"
)
