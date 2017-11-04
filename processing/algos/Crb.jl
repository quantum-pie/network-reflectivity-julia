module Crb

using Descriptors, AlphaBeta
export calc_crb

function radar_crb(proj_z::Matrix{Float64}, weights::Array{Float64, 3},
                    masks::Matrix{Bool}, alpha_beta::KtoZ,
                    radar_params::RadarParams)
  # useful const
  gam = 0.46 * alpha_beta.alpha * radar_params.range_elem * 1e-3

  # using eps to avoid NaNs
  proj_z += eps(Float64)

  # calculate one sided braces masking bad grid points
  Z_braces = weights .* (proj_z.^(alpha_beta.beta - 1)  .* masks)
  Z_braces = weights ./ proj_z - gam * alpha_beta.beta * cumsum(Z_braces, 2)
  Z_braces .*= masks

  pts = size(weights, 3)
  J = zeros(pts, pts)
  for k in 1:pts
      # calculate jacobi elements on k-th row
      J[k, k:end] = squeeze(sum(sum(Z_braces[:,:,k] .* Z_braces[:,:,k:end], 1), 2), (1, 2))
  end
  J
end

function calc_crb(proj_z::Vector{Matrix{Float64}}, weights::Vector{Array{Float64, 3}},
                  masks::Vector{Matrix{Bool}}, alpha_beta::KtoZ,
                  radar_params::RadarParams, cell_geometry::CellGeometry)
  J = 0
  for i in 1:cell_geometry.radars
      J += radar_crb(proj_z[i], weights[i], masks[i], alpha_beta, radar_params)
  end

  crb = diag(inv(radar_params.accum_pulses * (J + triu(J, 1)')))
end

end
