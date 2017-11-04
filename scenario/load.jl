if nprocs() < Sys.CPU_CORES
    addprocs(Sys.CPU_CORES - nprocs())
end

push!(LOAD_PATH, @__DIR__)

using Measurements, AlphaBeta

Plots.plotlyjs()
grant_meas()
