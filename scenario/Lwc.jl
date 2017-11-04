__precompile__(false)

module Lwc

using FileIO, JLD2, NLopt, PDMats, Plots, Distributions
export grant_lwc

function construct{T <: Real, W <: Real}(focuses_pos::Vector{Pair{T, T}},
                                          mixture_components::Integer, w_max::W;
                                          debug::Bool=false)
  focuses = length(focuses_pos)

  # sort points by distance from zero
  sort!(focuses_pos, by = p -> first(p).^2 + last(p).^2)

  mixture_components_array::Vector{DiagNormal} = []
  mixture_priority_array::Vector{Float64} = []
  centers_array::Vector{Vector{Float64}} = []

  for i in 1:focuses
    prop = 1 / (4^(i - 1))
    pos_std = 3000 / i
    pos_var_mean = 3000^2 / i
    pos = focuses_pos[i]

    for c in 1:mixture_components
      mu = [first(pos), last(pos)] + randn(2) * pos_std
      corr_diag_elem = randexp() * pos_var_mean
      push!(mixture_components_array, DiagNormal(mu, PDiagMat(fill(corr_diag_elem, 2))))
      push!(mixture_priority_array, sqrt(corr_diag_elem))
      push!(centers_array, mu)
    end
  end

  mixture_priority_array ./= sum(mixture_priority_array)
  mixture = MixtureModel(mixture_components_array, Categorical(mixture_priority_array))

  opt = Opt(:LD_LBFGS, 2)
  global_max = -Inf
  obj(x, grad) = pdf(mixture, x)
  max_objective!(opt, obj)
  for center in centers_array
    (f_max, x_max, ret) = optimize(opt, center)
    if f_max > global_max
      global_max = f_max
    end
  end

  if debug
    x = linspace(-8000, 8000, 100)
    display(surface(x, x, pdf.(mixture, broadcast(vcat, x, x'))' / global_max * w_max))
  end

  p::Pair{Float64, Float64} -> pdf(mixture, [first(p), last(p)]) / global_max * w_max
end

function grant_lwc(force_refresh::Bool=false)
  file_path = joinpath(@__DIR__, "..\\res\\w_fun.jld2")
  if !isfile(file_path) || force_refresh
    w_fun = construct([Pair(-3000, -3000), Pair(0, 0), Pair(5000, 5000)], 4, 5, debug=true)
    save(file_path, "w_fun", w_fun)
  else
    w_fun = load(file_path, "w_fun")
  end

  w_fun
end

end
