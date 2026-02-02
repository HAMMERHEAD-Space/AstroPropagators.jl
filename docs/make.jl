using Documenter
using AstroPropagators

makedocs(;
    modules=[AstroPropagators],
    format=Documenter.HTML(;
        prettyurls=(!("local" in ARGS)), highlights=["yaml"], ansicolor=true
    ),
    sitename="AstroPropagators.jl",
    authors="Jordan Murphy",
    pages=[
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "Propagators" => Any[
            "propagators/index.md",
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

deploydocs(; repo="github.com/HAMMERHEAD-Space/AstroPropagators.jl.git", target="build")
