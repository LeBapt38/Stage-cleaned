## Construction of the amplitude and method to project on plane waves

#%% Dependencies

using SpecialFunctions
using HypergeometricFunctions

# Mathematica shorthands/precisions.  BigFloat's precision is binary; use the
# requested decimal digits converted to a conservative number of bits whenever a
# high-precision constant/grid is constructed.
# Precision reduced for testing; the original values are GP=256, WP=512, MP=32.
const GP = 60
const WP = 120
const MP = 32
const GP_BITS = ceil(Int, GP * log2(10))
const WP_BITS = ceil(Int, WP * log2(10))
const IEPSILON = setprecision(GP_BITS) do
    complex(big"0", big(10.0)^(-WP))
end

"""
    AlphaIndex(i, j, n, m)

Numerical replacement for the Mathematica symbol `Alpha[i,j,n,m]`, which was
renamed there to symbols like `Alpha1x2x0x1`.
"""
struct AlphaIndex
    i::Int
    j::Int
    n::Int
    m::Int
end

struct AmplitudeTerm{T}
    alpha::AlphaIndex
    value::T
end

#%% Useful functions

_ascomplex(x) = x isa Complex ? x : complex(x, zero(x))
_msqrt(x) = sqrt(_ascomplex(x))      # principal branch, avoiding Julia's real negative DomainError
_mlog(x) = log(_ascomplex(x))
_mpow(x, a) = _ascomplex(x)^a

"Mathematica `Chop[#, 10^-GP]&` followed by high precision, applied elementwise."
function _chop_gp(x)
    tol = big(10.0)^(-GP)
    if x isa Complex
        re = abs(real(x)) < tol ? zero(real(x)) : real(x)
        im = abs(imag(x)) < tol ? zero(imag(x)) : imag(x)
        return complex(re, im)
    else
        return abs(x) < tol ? zero(x) : x
    end
end

function _sort_key(x)
    z = _ascomplex(x)
    return (real(z), imag(z))
end

function set_prec_grid(values)
    setprecision(WP_BITS) do
        out = [_chop_gp(v) for v in vec(values)]
        out = unique(out)
        sort!(out, by=_sort_key)
        return out
    end
end

"Legendre-Q-like kernel used later by the Froissart-Gribov projection."
function polQ(J, d, z)
    zc = _ascomplex(z)
    return (4 / pi) *
        (sqrt(pi) * gamma(J + 1) * gamma((d - 2) / 2)) /
        (2^(J + 1) * gamma(J + (d - 1) / 2) * zc^(J + 1)) *
        HypergeometricFunctions._₂F₁((J + 1) / 2, (J + 2) / 2, J + (d - 1) / 2, inv(zc^2))
end

"Numerical derivative corresponding to `D[PolQ[J,d,zz], zz] /. zz -> z`."
function polQ_derivative(J, d, z; h=nothing)
    zc = _ascomplex(z)
    hh = h === nothing ? sqrt(eps(float(real(one(zc))))) * (one(zc) + abs(zc)) : h
    return (polQ(J, d, zc + hh) - polQ(J, d, zc - hh)) / (2hh)
end

z0t(s) = (s + 4) / (s - 4)
z0u(s) = (s + 4) / (4 - s)
tt(s, z) = ((4 - s) / 2) * (1 - z)
uu(s, z) = ((4 - s) / 2) * (1 + z)

half_compactify(s0, phi) = (2s0) / (1 + cos(phi))
compact_jacobian(s0, phi) = (2s0 * sin(phi)) / (1 + cos(phi))^2
compact_integrand_value(f, s0, phi) = compact_jacobian(s0, phi) * f(half_compactify(s0, phi))

#%% Construction of the convergent piece

rho(s, sigma) = (_msqrt(sigma - 4) - _msqrt(4 - s)) / (_msqrt(sigma - 4) + _msqrt(4 - s))
disc_rho(s, sigma) = (2 * _msqrt(s - 4) * _msqrt(sigma - 4)) / (s + sigma - 8)
inv_rho_s(r, sigma) = 4 - (sigma - 4) * ((1 - r) / (1 + r))^2

tx(s, x) = ((s - 4) / 2) * (x - 1)
ux(s, x) = 4 - s - tx(s, x)

const CSP = big(4) / big(3)
const MOTHERPOINT = big(20) / big(3)

function phi_grid(n::Integer)
    setprecision(WP_BITS) do
        set_prec_grid([(big(pi) / 2) * (1 + cos(big(pi) * i / (n + 1))) for i in n:-1:1])
    end
end

rho_grid(n::Integer) = set_prec_grid([exp(complex(big"0", phi)) for phi in phi_grid(n)])
s_grid(n::Integer) = set_prec_grid([inv_rho_s(r, MOTHERPOINT) for r in rho_grid(n)]) .+ IEPSILON

function t_grid(n::Integer)
    nhalf = round(Int, n / 2)
    phis = vcat(phi_grid(nhalf), phi_grid(2000)[end - nhalf + 1:end])
    return set_prec_grid((4 / big(pi)) .* phis)
end

const _sigma_grid_cache = Dict{Int, Vector{BigFloat}}()
function sigma_grid(n::Integer)
    n < 1 && throw(ArgumentError("sigma_grid is defined for n >= 1"))
    return get!(_sigma_grid_cache, n) do
        if n == 1
            return set_prec_grid(real.(s_grid(1)))
        end

        previous = sigma_grid(n - 1)
        for i in 1:100
            for candidate in real.(s_grid(i))
                if !(candidate in previous)
                    return set_prec_grid(vcat(previous, candidate))
                end
            end
        end
        error("sigma_grid error ...")
    end
end

p_grid(n::Integer) = ones(Int, length(sigma_grid(n)))
ell_grid(n::Integer) = collect(0:2:n)

phi_norm(d, s) = _mpow(s - 4, (d - 3) / 4) * _mpow(s, -1 / 4)
Nd(d) = (16pi)^((2 - d) / 2) / gamma((d - 2) / 2)
nl(d, l::Integer) = l == 0 ?
    2^(-3 + 2d) * pi^(0.5 * (-3 + d)) * gamma(0.5 * (-1 + d)) :
    ((4pi)^(d / 2) * (d + 2l - 3) * gamma(d + l - 3)) /
        (pi * gamma((d - 2) / 2) * gamma(l + 1))
# nd(d) = big(2.0)^(5 - d) * 2^(-3 + 2d) * pi^(0.5 * (-3 + d)) * gamma(0.5 * (-1 + d))
function nd(d::Integer)
    setprecision(WP_BITS) do
        big(2.0)^(5 - d) *
        big(2.0)^(-3 + 2d) *
        big(pi)^((big(d) - 3) / 2) *
        gamma((big(d) - 1) / 2)
    end
end

const _ORDERED_CHANNEL_PAIRS = ((1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2))

function _channel_value(channel::Int, s, t, u)
    channel == 1 && return s
    channel == 2 && return t
    channel == 3 && return u
    throw(ArgumentError("unknown channel $channel"))
end

"Crossing-symmetric rho ansatz: 1/6 sum over Permutations[{s,t,u},{2}]."
function rho_ansatz(s, t, u, n::Integer, m::Integer, sigma_i, sigma_j)
    total = zero(rho(s, sigma_i)^n * rho(t, sigma_j)^m)
    for (a, b) in _ORDERED_CHANNEL_PAIRS
        total += rho(_channel_value(a, s, t, u), sigma_i)^n * rho(_channel_value(b, s, t, u), sigma_j)^m
    end
    return total / 6
end

# These discontinuity rules intentionally reproduce the Mathematica replacement
# implementation and are valid only for n,m in {0,1}.  The replacement is made on
# the already-specialised symbolic ansatz
#     rho[T] -> disc_rho[T], rho[S] -> 0, rho[U] -> 0,
#     rho[T]rho[S] -> disc_rho[T]rho[S],
#     rho[T]rho[U] -> disc_rho[T]rho[U], rho[U]rho[S] -> 0.
# function _disc_rho_ansatz(channel::Int, s, t, u, n::Integer, m::Integer, sigma_i, sigma_j)
#     (n in (0, 1) && m in (0, 1)) ||
#         throw(ArgumentError("disc_rho_ansatz matches the Mathematica code only for n,m in {0,1}"))

#     x = _channel_value(channel, s, t, u)
#     total = zero(disc_rho(x, sigma_i))
#     if n == 0 && m == 0
#         return total
#     end

#     for (a, b) in _ORDERED_CHANNEL_PAIRS
#         if n == 1 && m == 0
#             a == channel && (total += disc_rho(_channel_value(a, s, t, u), sigma_i))
#         elseif n == 0 && m == 1
#             b == channel && (total += disc_rho(_channel_value(b, s, t, u), sigma_j))
#         elseif a == channel && b != channel
#             total += disc_rho(_channel_value(a, s, t, u), sigma_i) * rho(_channel_value(b, s, t, u), sigma_j)
#         elseif b == channel && a != channel
#             total += rho(_channel_value(a, s, t, u), sigma_i) * disc_rho(_channel_value(b, s, t, u), sigma_j)
#         end
#     end
#     return total / 6
# end

function discT_rho_ansatz(s, t, u, n::Integer, m::Integer, sigma_i, sigma_j)
    (n in (0, 1) && m in (0, 1)) ||
        throw(ArgumentError("disc_rho_ansatz matches the Mathematica code only for n,m in {0,1}"))

    total = zero(disc_rho(_channel_value(2, s, t, u), sigma_i))
    if n == 0 && m == 0
        return total
    end

    for (a, b) in _ORDERED_CHANNEL_PAIRS
        if n == 1 && m == 0 && a==2
            total += disc_rho(_channel_value(2, s, t, u), sigma_i)
        elseif n == 0 && m == 1 && b==2
            total += disc_rho(_channel_value(2, s, t, u), sigma_j)
        elseif n == 1 && m == 1 && a==2
            total += disc_rho(_channel_value(2, s, t, u), sigma_i) * rho(_channel_value(b, s, t, u), sigma_j)
        elseif n == 1 && m == 1 && b==2
            total += rho(_channel_value(a, s, t, u), sigma_i) * disc_rho(_channel_value(2, s, t, u), sigma_j)
        end
    end
    return total / 6
end

# discT_rho_ansatz(s, t, u, n::Integer, m::Integer, sigma_i, sigma_j) = _disc_rho_ansatz(2, s, t, u, n, m, sigma_i, sigma_j)
# discU_rho_ansatz(s, t, u, n::Integer, m::Integer, sigma_i, sigma_j) = _disc_rho_ansatz(3, s, t, u, n, m, sigma_i, sigma_j)

function cond(i::Integer, j::Integer, n::Integer, m::Integer, sigma_i, sigma_j, pg_i::Integer, pg_j::Integer)
    return Int(abs(sigma_i) <= abs(sigma_j) &&
        n <= m &&
        (m != 0 || j <= 1) &&
        (n != 0 || i <= 1) &&
        n <= min(pg_i, pg_j) &&
        m <= max(pg_i, pg_j))
end

function Mraw(d, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer})
    terms = AmplitudeTerm[]
    for i in eachindex(sigmas), j in eachindex(sigmas), n in 0:pgs[i], m in 0:pgs[j]
        cnd = cond(i, j, n, m, sigmas[i], sigmas[j], pgs[i], pgs[j])
        push!(terms, AmplitudeTerm(AlphaIndex(i, j, n, m), nd(d) * rho_ansatz(s, t, u, n, m, sigmas[i], sigmas[j]) * cnd))
    end
    return terms
end

function discTMraw(d, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer})
    terms = AmplitudeTerm[]
    for i in eachindex(sigmas), j in eachindex(sigmas), n in 0:pgs[i], m in 0:pgs[j]
        cnd = cond(i, j, n, m, sigmas[i], sigmas[j], pgs[i], pgs[j])
        push!(terms, AmplitudeTerm(AlphaIndex(i, j, n, m), nd(d) * discT_rho_ansatz(s, t, u, n, m, sigmas[i], sigmas[j]) * cnd))
    end
    return terms
end

# function discUMraw(d, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer})
#     terms = AmplitudeTerm[]
#     for i in eachindex(sigmas), j in eachindex(sigmas), n in 0:pgs[i], m in 0:pgs[j]
#         cnd = cond(i, j, n, m, sigmas[i], sigmas[j], pgs[i], pgs[j])
#         push!(terms, AmplitudeTerm(AlphaIndex(i, j, n, m), nd(d) * discU_rho_ansatz(s, t, u, n, m, sigmas[i], sigmas[j]) * cnd))
#     end
#     return terms
# end

#%% Construction of the divergent piece

_isodd_integer(d) = isodd(Integer(d))

function Delta(d::Integer, s)
    exponent = -((d - 3) / 2)
    pref = _mpow((4 - s) / 4, exponent)
    if _isodd_integer(d)
        return pref * (pi * (-1)^((d - 3) ÷ 2)) / _mlog((4 - s) / 4) * rho(s, 8)^2
    else
        return pref * sin(pi * (d - 3) / 2)
    end
end

Delta(d::Integer, s, t, u) = Delta(d, s) + Delta(d, t) + Delta(d, u)
DeltaM(d::Integer, s, t, u) = [AmplitudeTerm(AlphaIndex(0, 0, 0, 0), nd(d) * Delta(d, s, t, u))]

pstep(d::Integer) = ceil((d - 3) / 2) - (d - 3) / 2

function p_powers(d::Integer)
    if !_isodd_integer(d)
        step = pstep(d)
        vals = Float64[]
        x = -((d - 3) / 2) + step
        while x <= -1e-9 + 10eps(Float64)
            push!(vals, x)
            x += step
        end
        return vals
    else
        return collect(1:((d - 3) ÷ 2))
    end
end

function delta_Delta(d::Integer, s, p)
    if _isodd_integer(d)
        return _mpow((4 - s) / 4, -((d - 3) / 2) + p) *
            (pi * (-1)^((d - 3) ÷ 2)) / _mlog((4 - s) / 4) * rho(s, 8)^2
    else
        return _mpow(4 - s, p)
    end
end

delta_Delta(d::Integer, s, t, u, p) = delta_Delta(d, s, p) + delta_Delta(d, t, p) + delta_Delta(d, u, p)

function deltaM(d::Integer, s, t, u)
    ps = p_powers(d)
    return [AmplitudeTerm(AlphaIndex(0, 0, 0, i), nd(d) * delta_Delta(d, s, t, u, ps[i])) for i in eachindex(ps)]
end

#%% Putting it together and expression in term of s, z

# This follows the active Mathematica definition exactly: the convergent Mraw
# term is commented out in `M`, while the discontinuities use only the convergent
# rho-ansatz discontinuities.
M(d::Integer, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer}) = vcat(Mraw(d,s,t,u, sigmas, pgs),DeltaM(d, s, t, u), deltaM(d, s, t, u))
discTM(d::Integer, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer}) = discTMraw(d, s, t, u, sigmas, pgs)
# discUM(d::Integer, s, t, u, sigmas::AbstractVector, pgs::AbstractVector{<:Integer}) = discUMraw(d, s, t, u, sigmas, pgs)

myM(d::Integer, Nmax::Integer, s, t, u) = M(d, s, t, u, sigma_grid(Nmax), p_grid(Nmax))
myDiscTM(d::Integer, Nmax::Integer, s, t, u) = discTM(d, s, t, u, sigma_grid(Nmax), p_grid(Nmax))
# myDiscUM(d::Integer, Nmax::Integer, s, t, u) = discUM(d, s, t, u, sigma_grid(Nmax), p_grid(Nmax))

_coeff_get(coeffs, idx::AlphaIndex, default) = get(coeffs, idx, get(coeffs, (idx.i, idx.j, idx.n, idx.m), default))

"""
    evaluate_terms(terms, coeffs; strict=false)

Numerical counterpart of Mathematica's coefficient replacement followed by a dot
product. `coeffs` can be keyed either by `AlphaIndex(i,j,n,m)` or by the tuple
`(i,j,n,m)`. Missing coefficients are treated as zero by default, which is
convenient for terms killed by `cond`; set `strict=true` to require every
coefficient explicitly.
"""
function evaluate_terms(terms::AbstractVector{<:AmplitudeTerm}, coeffs; strict::Bool=false)
    isempty(terms) && return 0
    total = zero(first(terms).value)
    for term in terms
        coeff = _coeff_get(coeffs, term.alpha, nothing)
        if coeff === nothing
            strict && throw(KeyError(term.alpha))
            coeff = zero(total)
        end
        total += coeff * term.value
    end
    return total
end

# The Mathematica notebook hard-codes Nmax=20 in these definitions; keep that as
# the default while allowing it to be overridden for tests.
function exprSTU(d::Integer, coeffs, s, t, u; Nmax::Integer=20, strict::Bool=false)
    return evaluate_terms(myM(d, Nmax, s, t, u), coeffs; strict=strict) / (32*big(pi))
end

function discTexprSTU(d::Integer, coeffs, s, t, u; Nmax::Integer=20, strict::Bool=false)
    return evaluate_terms(myDiscTM(d, Nmax, s, t, u), coeffs; strict=strict) / (32*big(pi))
end

# function discUexprSTU(d::Integer, coeffs, s, t, u; Nmax::Integer=20, strict::Bool=false)
#     return evaluate_terms(myDiscUM(d, Nmax, s, t, u), coeffs; strict=strict) / (32pi)
# end

expr(d::Integer, coeffs, s, z; Nmax::Integer=20, strict::Bool=false) =
    exprSTU(d, coeffs, s, tt(s, z), uu(s, z); Nmax=Nmax, strict=strict)

discTExpr(d::Integer, coeffs, s, z; Nmax::Integer=20, strict::Bool=false) =
    discTexprSTU(d, coeffs, s, tt(s, z), uu(s, z); Nmax=Nmax, strict=strict)

discUExpr(d::Integer, coeffs, s, z; Nmax::Integer=20, strict::Bool=false) =
    discUexprSTU(d, coeffs, s, tt(s, z), uu(s, z); Nmax=Nmax, strict=strict)

