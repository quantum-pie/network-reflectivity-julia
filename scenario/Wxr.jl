module Wxr

using SpecialFunctions, QuadGK, Constants
export K, Z

const ddmin = 0
const ddmax = 0.6e-2

function mie_abcd{S <: Number, R <: Real}(m::S, x::R)
  #=
  % Computes a matrix of Mie coefficients, a_n, b_n, c_n, d_n,
  % of orders n=1 to nmax, complex refractive index m=m'+im",
  % and size parameter x=k0*a, where k0= wave number
  % in the ambient medium, a=sphere radius;
  % p. 100, 477 in Bohren and Huffman (1983) BEWI:TDD122
  % C. M?tzler, June 2002
  =#

  nmax = round(Int64, 2 + x + 4x^(1/3))
  n = collect(1:nmax)
  nu = (n + 0.5)
  z = m * x
  m2 = m^2
  sqx = sqrt(0.5 * pi / x)
  sqz = sqrt(0.5 * pi / z)
  bx = besselj.(nu, x) * sqx
  bz = besselj.(nu, z) * sqz
  yx = bessely.(nu, x) * sqx
  hx = bx + im*yx
  b1x = vcat(sin(x)/x, bx[1:nmax-1])
  b1z = vcat(sin(z)/z, bz[1:nmax-1])
  y1x = vcat(-cos(x)/x, yx[1:nmax-1])
  h1x = b1x + im*y1x
  ax = x .* b1x - n .* bx
  az = z .* b1z - n .* bz
  ahx = x .* h1x - n .* hx
  an = (m2 .* bz .* ax - bx .* az)  ./ (m2 .* bz .* ahx - hx .* az)
  bn = (bz .* ax - bx .* az)        ./ (bz .* ahx - hx .* az)
  cn = (bx .* ahx - hx .* ax)       ./ (bz .* ahx - hx .* az)
  dn = m .* (bx .* ahx - hx .* ax)  ./ (m2 .* bz .* ahx - hx .* az)
  Dict(:a => an, :b => bn, :c => cn, :d => dn)
end

function mie{S <: Number, R <: Real}(m::S, x::R)::Vector{Float64}
  #=
  % Computation of Mie Efficiencies for given
  % complex refractive-index ratio m=m'+im"
  % and size parameter x=k0*a, where k0= wave number in ambient
  % medium, a=sphere radius, using complex Mie Coefficients
  % an and bn for n=1 to nmax,
  % s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
  % Result: m', m", x, efficiencies for extinction (qext),
  % scattering (qsca), absorption (qabs), backscattering (qb),
  % asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).
  % Uses the function "mie_abcd" for an and bn, for n=1 to nmax.
  % C. M?tzler, May 2002.
  =#

  if x == 0 # To avoid a singularity at x=0
    return [real(m), imag(m), 0, 0, 0, 0, 0, 0, 1.5]
  elseif x > 0 # This is the normal situation
    nmax = round(Int64, 2 + x + 4x^(1/3))
    n1 = nmax - 1
    n = collect(1:nmax)
    cn = 2n + 1
    c1n = n .* (n + 2) ./ (n + 1)
    c2n = cn ./ n ./ (n + 1)
    x2 = x^2
    f = mie_abcd(m, x)
    anp = real.(f[:a])
    anpp = imag.(f[:a])
    bnp = real.(f[:b])
    bnpp = imag.(f[:b])
    g11 = vcat(anp[2:nmax], 0)
    g12 = vcat(anpp[2:nmax], 0)
    g13 = vcat(bnp[2:nmax], 0)
    g14 = vcat(bnpp[2:nmax], 0)
    dn=cn .* (anp + bnp)
    q = sum(dn)
    qext = 2q / x2
    en = cn .* (anp .* anp + anpp .* anpp + bnp .* bnp + bnpp .* bnpp)
    q = sum(en)
    qsca = 2q / x2
    qabs = qext - qsca
    fn=(f[:a] - f[:b]) .* cn
    gn = (-1).^n
    f[:c] = fn .* gn
    q = sum(f[:c])
    qb = q * conj(q) / x2
    asy1 = c1n .* (anp .* g11 + anpp .* g12 + bnp .* g13 + bnpp .* g14)
    asy2 = c2n .* (anp .* bnp + anpp .* bnpp)
    asy = 4 / x2 * sum(asy1 + asy2) / qsca
    qratio = qb / qsca
    return [real(m), imag(m), x, qext, qsca, qabs, qb, asy, qratio]
  end
end

function drop_size_distr{S <: Real, R <: Real}(d::S, w::R)
  # Drop size distribution model, returns in cm^-4
  # d - drop diameter in cm
  # w - liquid water content in g/m^3

  6.92e-2 * w^0.038 * exp(-21.6 * w^(-0.24) * d)
end

function K{S <: Real, R <: Real}(w::S, lambda::R)
  #=
  % Attenuation coefficient in dB/km
  % W - liquid water content in g/m^3
  % lambda - radar wavelength in m
  % ddmin - minimum drop diameter in m
  % ddmax - maximum drop diameter in m
  =#

  wn = 2pi / lambda

  # function to integrate
  function fint{P <: Real}(d::P)
    mie_result = mie(Constants.m, wn * d * 1e-2 / 2)

    # extinction cross-section in cm^2
    sigb = mie_result[4] * pi * (d/2)^2
    drop_size_distr(d, w) * sigb
  end

  0.4343e6 * QuadGK.quadgk(fint, ddmin * 1e2, ddmax * 1e2)[1]
end

function Z{S <: Real, R <: Real}(w::S, lambda::R)
  #=
  % Radar reflectivity in mm^6 * m^-3
  % W - liquid water content in g/m^3
  % lambda - radar wavelength in m
  % ddmin - minimum drop diameter in m
  % ddmax - maximum drop diameter in m
  =#

  wn = 2pi / lambda

  # function to integrate
  function fint{P <: Real}(d::P)
    mie_result = mie(Constants.m, wn * d * 1e-2 / 2);

    # backscattering cross-section in cm^2
    sigb = mie_result[7] * pi * (d/2)^2;
    drop_size_distr(d, w) * sigb;
  end

  1e12 * (lambda * 1e2)^4 / (pi^5 * Constants.Km_sq_abs) * QuadGK.quadgk(fint, ddmin * 1e2, ddmax * 1e2)[1]
end

end
