module Measurements

using FileIO, JLD2, Distributions, QuadGK, Lwc, Wxr, Descriptors, Plots
export grant_meas

function generic_worker(x::Float64, angle::Float64, pos::Pair{Float64, Float64},
                        w_fun, K_hlpr)
  integrand(r) = K_hlpr(w_fun(Pair(pos[1] + cos(angle)*r, pos[2] + sin(angle)*r)))
  QuadGK.quadgk(integrand, 0, x)[1]
end

function noisify(Z::Float64, min_z::Float64, pulses::Integer, thr::Float64)
  res = Z/pulses * rand(Gamma(pulses, 1))
  res > thr ? res : 0
end

function generate_measurements(alpha::Vector{Float64}, r::Vector{Float64},
                                rad_desc::RadarParams, geom_desc::RadarGeometry,
                                w_fun, K_hlpr, Z_hlpr, realizations::Integer = 200;
                                debug::Bool = false)
  # calculate cartesian coordinates of polar grid points
  x = geom_desc.position[1] + cos.(alpha[geom_desc.beams_range])*r'
  y = geom_desc.position[2] + sin.(alpha[geom_desc.beams_range])*r'

  water_content = w_fun.(Pair.(x, y))

  # expensive shit
  Zt = pmap(Z_hlpr, water_content)

  Ka = zeros(size(Zt))
  for b in 1:geom_desc.beams
    angle = alpha[geom_desc.beams_range[b]]
    pos = geom_desc.position
    @time Ka[b, :] = pmap(x -> generic_worker(x, angle, pos, w_fun, K_hlpr), r)
  end

  Zt .*= exp.(-log(10)/5*Ka*1e-3) .* geom_desc.triangle_filter.(Pair.(x, y))

  Z = Vector{Matrix{Float64}}(realizations)
  for i in 1:realizations
    Z[i] = noisify.(Zt, rad_desc.min_z, rad_desc.accum_pulses, rad_desc.det_thr)
  end

  if debug
    reflec = 10log10.(Z[1])
    map!(x -> if isinf(x) NaN else x end, reflec, reflec)
    display(surface(reflec))
  end

  Z
end

function main()
  w_fun = grant_lwc()
  (radar_params, cell_geometry) = grant_desc()

  alpha = collect(0:(radar_params.beams - 1)) * radar_params.beam_width
  r = collect(1:radar_params.range_elems) * radar_params.range_elem

  wave_len = radar_params.wave_len
  K_closure(w::Float64) = K(w, wave_len)
  Z_closure(w::Float64) = Z(w, wave_len)

  Zm = Vector{Vector{Matrix{Float64}}}(cell_geometry.radars)
  for i in 1:cell_geometry.radars
    Zm[i] = generate_measurements(alpha, r, radar_params, cell_geometry.radars_geom[i], w_fun, K_closure, Z_closure, 2)
  end

  (collect(zip(Zm...)), Z_closure)
end

function grant_meas(force_refresh::Bool=false)
  meas_file_path = joinpath(@__DIR__, "..\\res\\meas.jld2")
  hlpr_file_path = joinpath(@__DIR__, "..\\res\\z_hlpr.jld2")
  if !isfile(meas_file_path) || force_refresh
    (Zm, Z_hlpr) = main()
    save(meas_file_path, "Zm", Zm)
    save(hlpr_file_path, "Z_hlpr", Z_hlpr)
  else
    Zm = load(meas_file_path, "Zm")
    Z_hlpr = load(hlpr_file_path, "Z_hlpr")
  end

  (Zm, Z_hlpr)
end

end
