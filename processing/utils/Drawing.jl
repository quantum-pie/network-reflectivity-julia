module Drawing

using Plots, Descriptors, ImageFiltering, FileIO, Images
export draw_crb, draw_z

function draw_crb(crb::Matrix{Float64}, z_max::Float64,
                    cell_geometry::CellGeometry)
  X_ = collect( (cell_geometry.xlim[1]):(cell_geometry.elem[1]):(cell_geometry.xlim[2]) )
  Y_ = collect( (cell_geometry.ylim[1]):(cell_geometry.elem[2]):(cell_geometry.ylim[2]) )

  crb = imfilter(crb, Kernel.gaussian(2)) .* (crb .!= 0)
  crb = 10log10.(sqrt.(crb) / z_max)
  crb[isinf.(crb)] = -60.0
  heatmap(X_ * 1e-3, Y_ * 1e-3, crb', c = :inferno, clim = (-55, -10), colorbar_title = "σ, dB")
  display(heatmap!(xlabel = "x, km", ylabel = "y, km"))
end

function draw_z(Z::Matrix{Float64}, cell_geometry::CellGeometry)
  im_width = 1000
  im_height = 1000
  X_ = linspace(cell_geometry.xlim[1], cell_geometry.xlim[2], im_width)
  Y_ = linspace(cell_geometry.ylim[1], cell_geometry.ylim[2], im_height)

  Z = imfilter(Z, Kernel.gaussian(1)) .* (Z .!= 0)
  Z = 10log10.(Z)
  Z[isinf.(Z)] = 0
  Z = imresize(Z, (im_width, im_height))

  colormap = load(joinpath(@__DIR__, "..\\..\\res\\nsw-grad.png"))
  heatmap(X_ * 1e-3, Y_ * 1e-3, Z', c = cgrad(vec(colormap)), clim = (0, 80), colorbar_title = "Z, dBZ")
  display(heatmap!(xlabel = "x, km", ylabel = "y, km"))
end

end
