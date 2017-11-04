module AlphaBeta

using FileIO, JLD2, Wxr, Descriptors, NLopt, Plots

export KtoZ, grant_ab

immutable KtoZ
  alpha::Float64
  beta::Float64
  alpha_db::Float64
  beta_db::Float64
end

function calc_alpha_beta(; debug::Bool=false)
  # Getting K-Z relation
  r_par = Descriptors.grant()[1]

  w_min = 0.005
  w_max = 5.0

  K_closure(w::Float64) = K(w, r_par.wave_len)
  Z_closure(w::Float64) = Z(w, r_par.wave_len)

  Q(x::Vector{Float64}, grad::Vector{Float64}) =
    QuadGK.quadgk(w::Float64 -> (K_closure(w) - x[1] * Z_closure(w)^x[2])^2, w_min, w_max)[1]

  Q_db(x::Vector{Float64}, grad::Vector{Float64}) =
    QuadGK.quadgk(w::Float64 -> (K_closure(w) - x[1] * (log10(Z_closure(w)))^x[2])^2, w_min, w_max)[1]

  if debug
    alp = linspace(0., 0.0002, 10)
    bet = linspace(0.75, 0.85, 10)

    tt = broadcast(vcat, alp, bet')
    vals = (v -> Q(v, Vector{Float64}([]))).(tt)
    display(surface(collect(bet), collect(alp), vals))

    alp = linspace(0.0, 4.0e-7, 10)
    bet = linspace(9.0, 10.0, 10)

    tt = broadcast(vcat, alp, bet')
    vals = (v -> Q_db(v, Vector{Float64}([]))).(tt)
    display(surface(collect(bet), collect(alp * 1e7), vals))
  end

  opt = Opt(:GN_DIRECT, 2)
  xtol_rel!(opt, 1e-4)
  ftol_rel!(opt, 1e-4)
  lower_bounds!(opt, [0.00005, 0.75])
  upper_bounds!(opt, [0.0002, 0.85])
  min_objective!(opt, Q)
  (f_min, x_min, ret) = optimize(opt, [0.00016, 0.828])
  print(f_min, x_min, ret)

  (kz_alpha, kz_beta) = (x_min...)

  opt = Opt(:GN_DIRECT, 2)
  xtol_rel!(opt, 1e-4)
  ftol_rel!(opt, 1e-4)
  lower_bounds!(opt, [0.0, 9.0])
  upper_bounds!(opt, [4e-7, 10.0])
  min_objective!(opt, Q_db)
  (f_min, x_min, ret) = optimize(opt, [3.1e-7, 9.78])
  print(f_min, x_min, ret)

  kz_alpha_db = x_min[1] * 10^(-x_min[2])
  kz_beta_db = x_min[2]

  KtoZ(kz_alpha, kz_beta, kz_alpha_db, kz_beta_db)
end

function grant_ab(force_refresh::Bool=false)
  file_path = joinpath(@__DIR__, "..\\res\\alpha_beta.jld2")
  if !isfile(file_path) || force_refresh
    alpha_beta = calc_alpha_beta()
    save(file_path, "alpha_beta", alpha_beta)
  else
    alpha_beta = load(file_path, "alpha_beta")
  end

  alpha_beta
end

end
