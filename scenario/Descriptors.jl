module Descriptors

using FileIO, JLD2, Constants
export RadarParams, RadarGeometry, CellGeometry, grant_desc

immutable RadarParams
  wave_len::Float64
  tr_power::Float64
  pulse_len::Float64
  ant_gain::Float64
  beam_width::Float64
  noise_pwr::Float64
  max_dist::Float64
  min_z::Float64
  det_thr::Float64
  range_elem::Float64
  beams::Integer
  range_elems::Integer
  accum_pulses::Integer
end

immutable RadarGeometry
  position::Pair{Float64, Float64}
  beams_range::UnitRange{Integer}
  beams::Integer
  density_filter::Function
  triangle_filter::Function
end

immutable CellGeometry
  radars::Integer
  radars_geom::Vector{RadarGeometry}
  xlim::Pair{Float64, Float64}
  ylim::Pair{Float64, Float64}
  elem::Pair{Float64, Float64}
end

function generate_descriptors()
  rad_wl = 0.032
  rad_tp = 30
  rad_pl = 19.5e-6
  rad_ag = 2700
  rad_bw = 3pi/180
  rad_np = Constants.N0 / rad_pl
  rad_md = 30e3
  rad_zm = 1e18 * rad_np * 2^10 * log(2) * rad_wl^2 * rad_md^2 /
            (pi^3 * rad_tp * rad_ag^2 * rad_bw^2 * Constants.c * rad_pl *
            Constants.Km_sq_abs)
  rad_dt = 10
  rad_re = 150
  rad_b = 2pi / rad_bw
  rad_res = rad_md / rad_re - 1
  rad_ap = 64

  radar_params = RadarParams(rad_wl, rad_tp, rad_pl, rad_ag, rad_bw, rad_np,
                            rad_md, rad_zm, rad_dt, rad_re, rad_b, rad_res,
                            rad_ap)

  h = sqrt(rad_md^2 - (rad_md / 2)^2)
  low_y = -h / 3
  high_y = h + low_y
  k1 = -sqrt(3) / 3;
  k2 = -k1;
  k3 = -sqrt(3);
  k4 = -k3;

  geom1 = RadarGeometry(Pair(-15e3, low_y), 1:21, 21,
                      (p::Pair{Float64, Float64}) -> ((p[1] <= 0.0) && (p[2] <= k1 * p[1])),
                      (p::Pair{Float64, Float64}) -> (p[2] < k3 * p[1] + high_y))

  geom2 = RadarGeometry(Pair(15e3, low_y), 41:61, 21,
                      (p::Pair{Float64, Float64}) -> ((p[1] > 0.0) && (p[2] <= k2 * p[1])),
                      (p::Pair{Float64, Float64}) -> (p[2] < k4 * p[1] + high_y))

  geom3 = RadarGeometry(Pair(0, high_y), 81:101, 21,
                      (p::Pair{Float64, Float64}) -> ( ((p[1] <= 0.0) && (p[2] >  k1 * p[1] )) ||
                                                    ((p[1] > 0.0) && (p[2] >  k2 * p[1] )) ),
                      (p::Pair{Float64, Float64}) -> (p[2] > low_y))

  cell_geometry = CellGeometry(3, [geom1, geom2, geom3],
                              Pair(-15e3, 15e3), Pair(-12e3, 18e3), Pair(rad_re, rad_re))

  (radar_params, cell_geometry)
end

function grant_desc(force_refresh::Bool=false)
  radar_file_path = joinpath(@__DIR__, "..\\res\\radar_params.jld2")
  geom_file_path = joinpath(@__DIR__, "..\\res\\cell_geometry.jld2")
  if !isfile(radar_file_path) || !isfile(geom_file_path) || force_refresh
    radar_params, cell_geometry = generate_descriptors()
    save(radar_file_path, "radar_params", radar_params)
    save(geom_file_path, "cell_geometry", cell_geometry)
  else
    radar_params = load(radar_file_path, "radar_params")
    cell_geometry = load(geom_file_path, "cell_geometry")
  end

  (radar_params, cell_geometry)
end

end
