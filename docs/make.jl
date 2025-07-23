using Documenter
using HAMMERHEAD

makedocs(
    modules = [HAMMERHEAD],
    format=Documenter.HTML(; prettyurls=!("local" in ARGS), highlights=["yaml"], ansicolor=true), 
    sitename = "HAMMERHEAD.jl",
    authors = "Jordan Murphy",
    pages = [
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "Propagators" => Any[
            "propagators/cowell.md",
            "propagators/gauss_ve.md",
            "propagators/edromo.md",
            "propagators/kustaanheimo_stiefel.md",
            "propagators/stiefel_scheifele.md",
            "propagators/milankovich.md",
            "propagators/usm.md",
        ],
        "API Reference" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/jmurphy6895/HAMMERHEAD.jl.git",
    target = "build",
)
