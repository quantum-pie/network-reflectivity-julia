module Chandr

using Basic, Descriptors, AlphaBeta, Projections
export apply_chandr

function apply_chandr(Zm::Tuple{Vararg{Matrix{Float64}}}, cell_geometry::CellGeometry,
                radar_params::RadarParams, alpha_beta::KtoZ,
                work_pts::Vector{Pair{Float64, Float64}})
  fir_filtered = map(v -> Basic.nfir(v, alpha_beta, radar_params), Zm)
  X = combine_projections_var(fir_filtered, work_pts, cell_geometry, radar_params, alpha_beta)
end

end
