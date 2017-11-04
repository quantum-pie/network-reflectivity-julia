module Projections

using Descriptors, GeometricalPredicates, VoronoiDelaunay, Plots,
        Wxr, AlphaBeta, SpecialFunctions

import Base.convert

export grid_interp, interp_mtx, projected_z, project_to_grid, project_to_space,
        combine_projections_dist, combine_projections_var

# Convert a Point2D to a vector
convert(::Type{Vector{Float64}}, p::Point2D) =
        [getx(p), gety(p)]

convert(::Type{Vector{Float64}}, p::Pair{Float64, Float64}) =
        [first(p), last(p)]

promote_rule(::Type{Vector{Float64}}, ::Type{Pair{Float64, Float64}}) = Vector{Float64}

function scale_value(x0::Float64, lim_old::Pair{Float64, Float64},
                      lim_new::Pair{Float64, Float64})
  x_new = (x0 - lim_old[1]) / (lim_old[2] - lim_old[1]) *
          (lim_new[2] - lim_new[1]) + lim_new[1]
end

function generate_scalers(cell_geometry::CellGeometry)
  # project points to tesselation working space
  predicate_lims = Pair(VoronoiDelaunay.min_coord, VoronoiDelaunay.max_coord)
  xscaler(x) = scale_value(x, cell_geometry.xlim, predicate_lims)

  # magic to ensure that all local grid point are inside good box
  # (basic lims in cell_geometry are about triangle limits)
  y_insurance_addition = 6000
  fixed_ylim = Pair(cell_geometry.ylim[1] - y_insurance_addition, cell_geometry.ylim[2])
  yscaler(y) = scale_value(y, fixed_ylim, predicate_lims)
  (xscaler, yscaler)
end

function barycentric(a::Point2D, b::Point2D, c::Point2D, p::Point2D)
  av = convert(Vector{Float64}, a)
  bv = convert(Vector{Float64}, b)
  cv = convert(Vector{Float64}, c)
  pv = convert(Vector{Float64}, p)

  v0 = bv - av
  v1 = cv - av
  v2 = pv - av

  d00 = dot(v0, v0)
  d01 = dot(v0, v1)
  d11 = dot(v1, v1)
  d20 = dot(v2, v0)
  d21 = dot(v2, v1)

  denom = d00 * d11 - d01 * d01
  v = (d11 * d20 - d01 * d21) / denom
  w = (d00 * d21 - d01 * d20) / denom
  u = 1.0 - v - w

  [u, v, w]
end

function grid_interp(x_from::Vector{Float64}, y_from::Vector{Float64},
                      z_from::Vector{Float64}, x_to::Vector{Float64},
                      y_to::Vector{Float64}, xscaler, yscaler)
  tess = DelaunayTessellation(length(x_from))

  pt_to_z = Dict{Point2D, Float64}()
  for (idx, (x, y)) in enumerate(zip(x_from, y_from))
    x_new, y_new = xscaler(x), yscaler(y)
    normalized_point = Point(x_new, y_new)
    pt_to_z[normalized_point] = z_from[idx]
    push!(tess, normalized_point)
  end

  out_z = fill(NaN, length(x_to))
  for (idx, (x, y)) in enumerate(zip(x_to, y_to))
    xn, yn = xscaler(x), yscaler(y)
    pt = Point(xn, yn)
    t = locate(tess, pt)
    if !isexternal(t)
      a = geta(t)
      b = getb(t)
      c = getc(t)

      bc = barycentric(a, b, c, pt)

      out_z[idx] = bc[1] * pt_to_z[a] + bc[2] * pt_to_z[b] +
                    bc[3] * pt_to_z[c]
    end
  end

  out_z
end

function interp_mtx_dbg(tess, cell_geometry::CellGeometry, radar_params::RadarParams,
                        dbg_radar::Integer, xscaler, yscaler)
  x, y = getplotxy(delaunayedges(tess))
  plot(x, y)

  geom_desc = cell_geometry.radars_geom[dbg_radar]
  alpha = collect((geom_desc.beams_range - 1) * radar_params.beam_width)
  r = collect(0:(radar_params.range_elems - 1)) * radar_params.range_elem + radar_params.range_elem/2
  x = geom_desc.position[1] + cos.(alpha)*r'
  y = geom_desc.position[2] + sin.(alpha)*r'

  good_pts = geom_desc.triangle_filter.(Pair.(x, y))
  x .*= good_pts
  y .*= good_pts

  display(scatter!(xscaler.(x), yscaler.(y), markercolor = :red,
          label = "", markersize = 1))
end

function interp_mtx(work_pts::Vector{Pair{Float64, Float64}},
                    cell_geometry::CellGeometry, radar_params::RadarParams;
                    debug::Bool=false)
  # initialize Delaunau tesselation object
  tess = DelaunayTessellation(length(work_pts))
  (xscaler, yscaler) = generate_scalers(cell_geometry)

  pt_to_idx = Dict{Point2D, Integer}()
  for (idx, pt) in enumerate(work_pts)
    x, y = pt[1], pt[2]
    x_new, y_new = xscaler(x), yscaler(y)
    normalized_point = Point(x_new, y_new)
    pt_to_idx[normalized_point] = idx
    push!(tess, normalized_point)
  end

  if debug
    interp_mtx_dbg(tess, cell_geometry, radar_params, 3, xscaler, yscaler)
  end

  weights = Vector{Array{Float64, 3}}(cell_geometry.radars)
  masks = Vector{Matrix{Bool}}(cell_geometry.radars)
  for i in 1:cell_geometry.radars
    geom_desc = cell_geometry.radars_geom[i]
    weights[i] = fill(0., (geom_desc.beams, radar_params.range_elems, length(work_pts)))
    masks[i] = fill(false, (geom_desc.beams, radar_params.range_elems))

    for (id, beam) in enumerate(geom_desc.beams_range)
      alpha = (beam - 1) * radar_params.beam_width
      r = collect(0:(radar_params.range_elems - 1)) * radar_params.range_elem + radar_params.range_elem/2
      x = geom_desc.position[1] + cos(alpha)*r
      y = geom_desc.position[2] + sin(alpha)*r

      for n in 1:radar_params.range_elems
        xn, yn = xscaler(x[n]), yscaler(y[n])
        pt = Point(xn, yn)
        t = locate(tess, pt)
        if !isexternal(t)
          a = geta(t)
          b = getb(t)
          c = getc(t)

          bc = barycentric(a, b, c, pt)

          weights[i][id, n, [pt_to_idx[a], pt_to_idx[b], pt_to_idx[c]]] = bc
          masks[i][id, n] = true
        end
      end
    end
  end

  (weights, masks)
end

function projected_z(cell_geometry::CellGeometry, radar_params::RadarParams,
                        Z_hlpr, w_fun)
  Z = Vector{Matrix{Float64}}(cell_geometry.radars)

  alpha = collect(0:(radar_params.beams - 1)) * radar_params.beam_width
  r = collect(1:radar_params.range_elems) * radar_params.range_elem

  for i in 1:cell_geometry.radars
    geom_desc = cell_geometry.radars_geom[i]
    x = geom_desc.position[1] + cos.(alpha[geom_desc.beams_range])*r'
    y = geom_desc.position[2] + sin.(alpha[geom_desc.beams_range])*r'

    Z[i] = pmap(Z_hlpr, w_fun.(Pair.(x, y)))
  end

  Z
end

function project_to_grid(Z::Vector{Float64}, cell_geometry::CellGeometry,
                          work_pts::Vector{Pair{Float64, Float64}})
  X = first.(work_pts)
  Y = last.(work_pts)

  X_ = collect( (cell_geometry.xlim[1]):(cell_geometry.elem[1]):(cell_geometry.xlim[2]) )
  Y_ = collect( (cell_geometry.ylim[1]):(cell_geometry.elem[2]):(cell_geometry.ylim[2]) )

  CV = collect(Iterators.product(X_, Y_))[:]

  (xscaler, yscaler) = generate_scalers(cell_geometry)

  Z_ = grid_interp(X, Y, Z, first.(CV), last.(CV), xscaler, yscaler)

  Z_[isnan.(Z_)] = 0
  Z_ = reshape(Z_, (length(X_), length(Y_)))
end

function project_to_space(Z::Matrix{Float64}, cell_geometry::CellGeometry,
                          geom_desc::RadarGeometry, radar_params::RadarParams,
                          pts::Vector{Pair{Float64, Float64}})
  alpha = collect(0:(radar_params.beams - 1)) * radar_params.beam_width
  r = collect(1:radar_params.range_elems) * radar_params.range_elem

  X = geom_desc.position[1] + cos.(alpha[geom_desc.beams_range])*r'
  Y = geom_desc.position[2] + sin.(alpha[geom_desc.beams_range])*r'

  X_ = first.(pts)
  Y_ = last.(pts)

  (xscaler, yscaler) = generate_scalers(cell_geometry)

  Z_ = grid_interp(X[:], Y[:], Z[:], X_, Y_, xscaler, yscaler)
  Z_[isnan.(Z_)] = 0
  Z_
end

function combine_projections_var(Z::Tuple{Vararg{Matrix{Float64}}}, pts::Vector{Pair{Float64, Float64}},
          cell_geometry::CellGeometry, radar_params::RadarParams, alpha_beta::KtoZ)
  X = zeros(length(pts))
  Wc = zeros(length(pts))

  for i in 1:cell_geometry.radars
      Zi = Z[i]
      Zi_db = 10log10.(Zi)
      Zi_db[Zi_db .== -Inf] = 0
      vars = zeros(size(Zi))
      vars[:, 1] = polygamma(1, radar_params.accum_pulses) * (10 / log(10))^2

      # calculate HB estimate variance
      for n in 2:radar_params.range_elems
          vars[:, n] = vars[:, n-1] .*
              (1 + 4 * (radar_params.range_elem * 1e-3)^2 * alpha_beta.alpha_db^2 *
              alpha_beta.beta_db^2 * Zi_db[:, n-1].^(2*alpha_beta.beta_db - 2))
      end

      X_proj = project_to_space(Zi, cell_geometry, cell_geometry.radars_geom[i], radar_params, pts)
      X_proj[X_proj .<= 0] = 1
      X_proj_db = 10log10.(X_proj)
      X_proj_db[isinf.(X_proj_db)] = 0
      vars_proj = project_to_space(vars, cell_geometry, cell_geometry.radars_geom[i], radar_params, pts)
      vars_proj[vars_proj .== 0] = polygamma(1, radar_params.accum_pulses) * (10 / log(10))^2
      Wnew = (X_proj_db .!= 0) ./ vars_proj
      X += X_proj_db .* Wnew
      Wc += Wnew
  end

  X = 10 .^ (X ./ Wc / 10)
end

function combine_projections_dist(Z::Tuple{Vararg{Matrix{Float64}}}, pts::Vector{Pair{Float64, Float64}},
          cell_geometry::CellGeometry, radar_params::RadarParams, alpha_beta::KtoZ)
  X = zeros(length(pts))
  Wc = zeros(length(pts))

  for i in 1:cell_geometry.radars
    geom_desc = cell_geometry.radars_geom[i]
    X_proj = project_to_space(Z[i], cell_geometry, geom_desc, radar_params, pts)
    dists = map(v -> sum(v.^2), convert.(Vector{Float64}, pts) .- [convert(Vector{Float64}, geom_desc.position)])
    Wnew = (X_proj .!= 0) ./ dists
    X += X_proj .* Wnew
    Wc += Wnew
  end

  X ./= Wc
end

end
