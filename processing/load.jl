if nprocs() < Sys.CPU_CORES
    addprocs(Sys.CPU_CORES - nprocs())
end

push!(LOAD_PATH, @__DIR__)
push!(LOAD_PATH, joinpath(@__DIR__, "..\\scenario\\"))
push!(LOAD_PATH, joinpath(@__DIR__, ".\\algos\\"))
push!(LOAD_PATH, joinpath(@__DIR__, ".\\utils\\"))

using Proc

Plots.gr()
run_processing()
