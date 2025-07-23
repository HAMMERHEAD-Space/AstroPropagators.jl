# ==============
# benchmark body
# ==============

const SUITE = BenchmarkGroup()
function tune_benchmarks!(g::BenchmarkGroup; seconds=nothing, gcsample=true)
    for v in values(g)
        if seconds !== nothing
            v.params.seconds = seconds
        end
        v.params.gcsample = gcsample
        v.params.evals = 1 # `setup` must be functional
    end
end

# warm up

# analysis performance
# --------------------

#TODO: SEE JET.jl
