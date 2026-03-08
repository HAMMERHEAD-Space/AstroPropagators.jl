using Documenter
using AstroPropagators
using AstroMeasurements

const _MeasExt = Base.get_extension(AstroPropagators, :AstroPropagatorsMeasurementsExt)

makedocs(;
    modules=[AstroPropagators, _MeasExt],
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
        ],
        "Events" => Any[
            "events/index.md",
            "events/orbital_detectors.md",
            "events/eclipse_detectors.md",
            "events/geometric_detectors.md",
            "events/utility_detectors.md",
            "events/maneuvers.md",
            "events/access.md",
        ],
        "API Reference" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo="github.com/HAMMERHEAD-Space/AstroPropagators.jl.git", target="build")
