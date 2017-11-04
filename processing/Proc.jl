module Proc

using FileIO, JLD2, AlphaBeta, Measurements, Descriptors, Lwc, Wxr, Drawing,
        Projections, Basic, Chandr, MaxLH, MinRMS, Plots, Crb

export run_processing

immutable Performance
  std_abs::Float64
  bias_abs::Float64
  std_db::Float64
  bias_db::Float64
end

function adaptive_discretize(cell_geometry::CellGeometry, radar_params::RadarParams)
  work_pts = Pair{Float64, Float64}[]
  alpha = collect(0:(radar_params.beams - 1)) * radar_params.beam_width
  r = collect(0:(radar_params.range_elems - 1)) * radar_params.range_elem + radar_params.range_elem/2

  for i in 1:cell_geometry.radars
    geom_desc = cell_geometry.radars_geom[i]
    x = geom_desc.position[1] + cos.(alpha[geom_desc.beams_range])*r'
    y = geom_desc.position[2] + sin.(alpha[geom_desc.beams_range])*r'

    good_pts = Iterators.filter(geom_desc.density_filter, Pair.(x, y))
    append!(work_pts, good_pts)
  end

  work_pts
end

function evaluate(Z_est::Vector{Float64}, Z_t::Vector{Float64},
                  cell_geometry::CellGeometry, radar_params::RadarParams,
                  pts::Vector{Pair{Float64, Float64}}, to_draw::Bool=false)
  std_abs = std(Z_t - Z_est)
  bias_abs = mean(Z_t - Z_est)
  std_db = std(10log10.(Z_t) - 10log10.(Z_est))
  bias_db = mean(10log10.(Z_t) - 10log10.(Z_est))

  Z_grid = project_to_grid(Z_est, cell_geometry, pts)
  Z_grid .* (Z_grid .>= radar_params.det_thr)

  if to_draw
      draw_z(Z_grid, cell_geometry)
  end

  perf = Performance(std_abs, bias_abs, std_db, bias_db)
  (Z_grid, perf)
end

function run_processing()
  # load various shit
  alpha_beta = grant_ab()
  (Zm, Z_hlpr) = grant_meas()
  w_fun = grant_lwc()
  (radar_params, cell_geometry) = grant_desc()

  # create main image cartesian grid
  X = collect( (cell_geometry.xlim[1]):(cell_geometry.elem[1]):(cell_geometry.xlim[2]) )
  Y = collect( (cell_geometry.ylim[1]):(cell_geometry.elem[2]):(cell_geometry.ylim[2]) )
  XY_meshgrid = Pair.(X, Y')

  # true z generation section
  tz_file_path = joinpath(@__DIR__, "..\\res\\true_z.jld2")
  if !isfile(tz_file_path)
    true_z = pmap(Z_hlpr, w_fun.(XY_meshgrid))
    save(tz_file_path, "true_z", true_z)
  else
    true_z = load(tz_file_path, "true_z")
  end

  # work points generation section
  pts_file_path = joinpath(@__DIR__, "..\\res\\work_pts.jld2")
  wz_file_path = joinpath(@__DIR__, "..\\res\\work_z.jld2")
  if !isfile(pts_file_path)
    work_pts = adaptive_discretize(cell_geometry, radar_params)
    work_z = pmap(Z_hlpr, w_fun.(work_pts))
    save(pts_file_path, "work_pts", work_pts)
    save(wz_file_path, "work_z", work_z)
  else
    work_pts = load(pts_file_path, "work_pts")
    work_z = load(wz_file_path, "work_z")
  end

  # calc true Z on local grids
  pz_file_path = joinpath(@__DIR__, "..\\res\\proj_z.jld2")
  if !isfile(pz_file_path)
    proj_z = projected_z(cell_geometry, radar_params, Z_hlpr, w_fun)
    save(pz_file_path, "proj_z", proj_z)
  else
    proj_z = load(pz_file_path, "proj_z")
  end

  # calc Cramer-Rao bound
  crb_file_path = joinpath(@__DIR__, "..\\res\\crb_abs.jld2")
  crb_pts_file_path = joinpath(@__DIR__, "..\\res\\crb_pts.jld2")
  if !isfile(crb_file_path)
    good = work_z .> radar_params.det_thr
    crb_pts = work_pts[good]
    crb_z = work_z[good]
    (crweights, crmasks) = interp_mtx(crb_pts, cell_geometry, radar_params)

    # filter abnormal points (not affecting any measurement)
    good = squeeze(sum(sum(sum(crweights), 1), 2), (1, 2)) .!= 0
    for i in 1:cell_geometry.radars
      crweights[i] = crweights[i][:, :, good]
    end

    crb_z = crb_z[good]
    crb_pts = crb_pts[good]

    @time crb_abs = calc_crb(proj_z, crweights, crmasks, alpha_beta, radar_params,
                        cell_geometry)

    save(crb_file_path, "crb_abs", crb_abs)
    save(crb_pts_file_path, "crb_pts", crb_pts)
  else
    crb_abs = load(crb_file_path, "crb_abs")
    crb_pts = load(crb_pts_file_path, "crb_pts")
  end

  draw_crb(project_to_grid(crb_abs, cell_geometry, crb_pts), maximum(true_z), cell_geometry)
  draw_z(true_z, cell_geometry)

  for meas in Zm
    @time X0c = apply_chandr(meas, cell_geometry, radar_params, alpha_beta, work_pts)

    # filter nonzero points based on initial guess and update weights
    good = X0c .> radar_params.det_thr
    X0 = X0c[good]
    pts = work_pts[good]
    rz = work_z[good]
    (rweights, rmasks) = interp_mtx(pts, cell_geometry, radar_params)

    (chandr_grid_Z, chandr_performance) = evaluate(X0, rz, cell_geometry, radar_params, pts, true)
    println("Chandr std: ", chandr_performance.std_db, " Chandr bias: ", chandr_performance.bias_db)

    @time X0s = apply_minrms(meas, rweights, rmasks, alpha_beta, cell_geometry, radar_params, 10, pts)
    (simple_grid_Z, simple_performance) = evaluate(X0s, rz, cell_geometry, radar_params, pts, true)
    println("Simple std: ", simple_performance.std_db, " Simple bias: ", simple_performance.bias_db)

    @time X0c = apply_maxlh(X0, meas, rweights, rmasks, alpha_beta, radar_params, cell_geometry)
    (complex_grid_Z, complex_performance) = evaluate(X0c, rz, cell_geometry, radar_params, pts, true)
    println("Complex std: ", complex_performance.std_db, " Complex bias: ", complex_performance.bias_db)
  end
end

end
