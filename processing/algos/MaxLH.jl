module MaxLH

using Basic, Descriptors, AlphaBeta, Projections, NLopt, Plots
export apply_maxlh

function radar_objective(X0::Vector{Float64}, grad::Vector{Float64},
                          Zm::Matrix{Float64}, W::Array{Float64, 3}, M::Matrix{Bool},
                          alpha_beta::KtoZ, radar_params::RadarParams)
  # interpolate prevoius estimate on polar grids
  Z_interp = squeeze(sum(X0 .* permutedims(W, [3, 1, 2]), 1), 1)

  # imitate detection
  Z_interp .*= (Z_interp .> radar_params.det_thr)

  # masks of "good" polar grid points
  M .&= (Zm .!= 0.0) .& (Z_interp .!= 0.0)

  # useful const
  gam = 0.46 * alpha_beta.alpha * radar_params.range_elem * 1e-3

  # using eps to avoid NaNs
  Z_interp += eps(Float64)

  Z_pred = Z_interp .* exp.(-gam * cumsum(Z_interp.^alpha_beta.beta, 2))

  # predicted/true relations
  Z_mult = Zm ./ Z_pred

  acc = radar_params.accum_pulses
  F0 = sum((Z_mult + log.(Z_pred)) .* M)

  if length(grad) > 0
    # calculate one_sided braces masking bad grid points
    Z_braces = W .* (Z_interp.^(alpha_beta.beta - 1)  .* M)
    Z_braces = W ./ Z_interp - gam * alpha_beta.beta * cumsum(Z_braces, 2)
    Z_braces .*= M

    J = squeeze(sum(sum((1 - Z_mult) .* Z_braces, 1), 2), (1, 2))
    grad[:] += log(10) / 10 * J .* X0
  end

  F0
end

function objective(X0::Vector{Float64}, grad::Vector{Float64},
                  Zm::Tuple{Vararg{Matrix{Float64}}}, weights::Vector{Array{Float64, 3}},
                  masks::Vector{Matrix{Bool}}, alpha_beta::KtoZ,
                  radar_params::RadarParams, cell_geometry::CellGeometry)
  X0 = 10 .^ (X0 / 10)
  F0 = 0

  for i in 1:cell_geometry.radars
    F0 += radar_objective(X0, grad, Zm[i], weights[i], masks[i], alpha_beta, radar_params)
  end

  println(F0)
  F0
end

function apply_maxlh(X0::Vector{Float64}, Zm::Tuple{Vararg{Matrix{Float64}}},
                rweights::Vector{Array{Float64, 3}}, rmasks::Vector{Matrix{Bool}},
                alpha_beta::KtoZ, radar_params::RadarParams,
                cell_geometry::CellGeometry)
  to_minimize(X, grad) = objective(X, grad, Zm, rweights, rmasks, alpha_beta, radar_params, cell_geometry)

  opt = Opt(:LD_LBFGS, length(X0))
  min_objective!(opt, to_minimize)
  maxeval!(opt, 20)
  lower_bounds!(opt, 10log10(radar_params.min_z) * ones(length(X0)))
  upper_bounds!(opt, 10log10.(X0) + 10)

  (f_min, X, ret) = try
    optimize(opt, 10log10.(X0))
  catch x
    if isa(x, InterruptException)
      interrupt(workers())
      throw(ErrorException("execution was interrupted"))
    else
      throw(x)
    end
  end

  println(f_min, " ", ret)
  X = 10 .^ (X / 10)
end

end
