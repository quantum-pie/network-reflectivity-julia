module Basic

using AlphaBeta, Descriptors
export nfir, stupid_slow_nonconverging_algorithm

function nfir(Zm::Matrix{Float64}, alpha_beta::KtoZ, radar_params::RadarParams)
  to_power = 1 - log(10) / 5 * alpha_beta.alpha * alpha_beta.beta *
        radar_params.range_elem * 1e-3 * cumsum(Zm.^alpha_beta.beta, 2)

  powered = map(v -> v < 0 ? Inf : v^(-1/alpha_beta.beta), to_power)
  Z = Zm .* powered
  Z[isnan.(Z)] = 0
  Z
end

function stupid_slow_nonconverging_algorithm(Zm::Matrix{Float64}, Zest::Matrix{Float64},
                                alpha_beta::KtoZ, radar_params::RadarParams)
  gam = log(10) / 5 * alpha_beta.alpha * radar_params.range_elem * 1e-3

  Z = Zm .* exp.(gam * cumsum(Zest.^alpha_beta.beta, 2))
  Z
end

end
