module MinRMS

using Basic, Descriptors, AlphaBeta, Projections
export apply_minrms

function apply_minrms(Zm::Tuple{Vararg{Matrix{Float64}}}, W::Vector{Array{Float64, 3}},
                M::Vector{Matrix{Bool}}, alpha_beta::KtoZ, cell_geometry::CellGeometry,
                radar_params::RadarParams, iterations::Integer,
                pts::Vector{Pair{Float64, Float64}})
  Zest = Zm
  X = Vector{Float64}()

  # run simple algorithm on polar grid
  for iter in 1:iterations
    Zest = map(z -> slow_nonconverging_algorithm(z[1], z[2], alpha_beta, radar_params), Tuple(zip(Zm, Zest)))
    X = combine_projections_dist(Zest, pts, cell_geometry, radar_params, alpha_beta)

    Zest_proj = Vector{Matrix{Float64}}(cell_geometry.radars)
    for i in 1:cell_geometry.radars
      Zest_proj[i] = squeeze(sum(X .* permutedims(W[i], [3, 1, 2]), 1), 1)
      Zest_proj[i] .*= (Zest_proj[i] .> 0)
    end
    Zest = Tuple(Zest_proj)
  end

  X
end

end
