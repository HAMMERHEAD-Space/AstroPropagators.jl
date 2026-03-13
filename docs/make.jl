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
            "propagators/geqoe.md",
            "propagators/modified_equinoctial.md",
        ],
        "Events" => Any[
            "events/index.md",
            "events/orbital_detectors.md",
            "events/eclipse_detectors.md",
            "events/geometric_detectors.md",
            "events/utility_detectors.md",
            "events/maneuvers.md",
        ],
        "API Reference" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo="github.com/HAMMERHEAD-Space/AstroPropagators.jl.git", target="build")
