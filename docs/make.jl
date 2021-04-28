using Documenter, VisAna, PyPlot

makedocs(
    sitename="VisAna",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = ["Home"           => "index.md",
             "Analysis"   => "analysis.md",
             "Example"   => "examples.md",
             "Log"       => "log.md",
             "Parker Spiral" => "parker_spiral.md"
    ])
)

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia.git",
    target = "build",
    branch = "gh-pages"
)
