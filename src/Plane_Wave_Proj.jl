## Construction of the amplitude and method to project on plane waves

#%% Dependencies

using SpecialFunctions
using HypergeometricFunctions
using QuadGK

# Mathematica shorthands/precisions.  BigFloat's precision is binary; use the
# requested decimal digits converted to a conservative number of bits whenever a
# high-precision constant/grid is constructed.
# Precision reduced for testing; the original values are GP=256, WP=512, MP=32.
const GP = 8
const WP = 32
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
polQ(J, d, z) = polQ(PolQKernel(J, d), z)

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
    model = cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict, numeric_type=:big)
    return exprSTU(model, s, t, u)
end

function discTexprSTU(d::Integer, coeffs, s, t, u; Nmax::Integer=20, strict::Bool=false)
    model = cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict, numeric_type=:big)
    return discTexprSTU(model, s, t, u)
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

#%% Froissart-Gribov projection

const P0 = parse(Int, get(ENV, "PW_PROJ_P", "4"))

"""
    my_n_integrate_faster(f, a, b, p=P0; maxevals=nothing)

QuadGK-backed replacement for the previous handwritten Gauss-Kronrod driver.
The precision is intentionally tied to the requested accuracy, so exploratory
runs can stay cheap (`PW_PROJ_P=4` by default) while validation can request the
old setting with `p=8` or `PW_PROJ_P=8`.
"""
function my_n_integrate_faster(f, a, b, p::Integer=P0; maxevals=nothing)
    if p <= 4
        aa = Float64(a)
        bb = Float64(b)
        tol = 10.0^(-p)
        effective_maxevals = maxevals === nothing ? 20_000 : maxevals
        value, _ = QuadGK.quadgk(f, aa, bb; rtol=tol, atol=tol,
            maxevals=effective_maxevals, order=3)
        return value
    end

    # High-precision evaluations should all run at the file-level working
    # precision WP.  Do not silently lower precision to MP or 4p inside the
    # quadrature path.
    setprecision(WP_BITS) do
        aa = BigFloat(a)
        bb = BigFloat(b)
        tol = big(10.0)^(-p)
        effective_maxevals = maxevals === nothing ? 10^7 : maxevals
        value, _ = QuadGK.quadgk(f, aa, bb; rtol=tol, atol=tol,
            maxevals=effective_maxevals)
        return value
    end
end

# Precomputed coefficient-weighted models used by the projection integrands.
# They keep the hypergeometric kernels and all branch-sensitive building blocks
# unchanged, but avoid rebuilding `AmplitudeTerm` arrays at every quadrature node.

struct WeightedRhoTerm{T,S}
    weight::T
    i::Int
    j::Int
    n::Int
    m::Int
    sigma_i::S
    sigma_j::S
end

struct WeightedThresholdTerm{T,P}
    weight::T
    p::P
end

struct ProjectionModel{S,Q,R,L,T,W,C}
    d::Int
    Nmax::Int
    sigmas::S
    sigma_roots::Q
    terms::R
    leading_threshold_weight::L
    subleading_threshold_terms::T
    threshold_weights::W
    coeffs::C
end

mutable struct ProjectionWorkspace{T}
    rs::Vector{T}
    rt::Vector{T}
    ru::Vector{T}
    dt::Vector{T}
end

struct PolQKernel{P,A,B,C,E}
    prefactor::P
    a::A
    b::B
    c::C
    exponent::E
end

struct PgenKernel{A,B,C}
    a::A
    b::B
    c::C
end

const _projection_model_cache = IdDict{Any, Dict{Tuple{Int, Int, Bool, Symbol}, Any}}()
const _projection_model_cache_lock = ReentrantLock()

function _coefficient_or_missing(coeffs, idx::AlphaIndex, strict::Bool)
    coeff = _coeff_get(coeffs, idx, nothing)
    if coeff === nothing
        strict && throw(KeyError(idx))
        return nothing
    end
    return coeff
end

function make_projection_model(d::Integer, coeffs; Nmax::Integer=20, strict::Bool=false,
                               numeric_type::Symbol=:big)
    return setprecision(numeric_type == :big ? WP_BITS : precision(BigFloat)) do
        sigmas = sigma_grid(Nmax)
        pgs = p_grid(Nmax)
        norm = nd(d) / (32 * big(pi))
        if numeric_type == :float64
            sigmas = Float64.(sigmas)
            norm = Float64(norm)
        elseif numeric_type != :big
            throw(ArgumentError("numeric_type must be :big or :float64, got $numeric_type"))
        end
        sigma_roots = [_msqrt(sigma - 4) for sigma in sigmas]

        weight_type = typeof(complex(norm, zero(norm)))
        sigma_type = eltype(sigmas)
        terms = WeightedRhoTerm{weight_type, sigma_type}[]

    for i in eachindex(sigmas), j in eachindex(sigmas), n in 0:pgs[i], m in 0:pgs[j]
        cnd = cond(i, j, n, m, sigmas[i], sigmas[j], pgs[i], pgs[j])
        cnd == 0 && continue

        idx = AlphaIndex(i, j, n, m)
        coeff = _coefficient_or_missing(coeffs, idx, strict)
        coeff === nothing && continue
        iszero(coeff) && continue

        weight = weight_type(_ascomplex(coeff) * norm * cnd)
        push!(terms, WeightedRhoTerm(weight, i, j, n, m, sigmas[i], sigmas[j]))
    end

    threshold_weights = Dict{Int, weight_type}()
    leading_coeff = _coefficient_or_missing(coeffs, AlphaIndex(0, 0, 0, 0), strict)
    leading_threshold_weight = leading_coeff === nothing ? nothing : weight_type(_ascomplex(leading_coeff) * norm)
    leading_threshold_weight !== nothing && (threshold_weights[0] = leading_threshold_weight)

    ps = p_powers(d)
    subleading_threshold_terms = WeightedThresholdTerm{weight_type, eltype(ps)}[]
    for k in eachindex(ps)
        idx = AlphaIndex(0, 0, 0, k)
        coeff = _coefficient_or_missing(coeffs, idx, strict)
        coeff === nothing && continue
        iszero(coeff) && continue
        weight = weight_type(_ascomplex(coeff) * norm)
        threshold_weights[k] = weight
        push!(subleading_threshold_terms, WeightedThresholdTerm(weight, ps[k]))
    end

        return ProjectionModel(Int(d), Int(Nmax), sigmas, sigma_roots, terms,
            leading_threshold_weight, subleading_threshold_terms, threshold_weights, coeffs)
    end
end

function cached_projection_model(d::Integer, coeffs; Nmax::Integer=20,
                                 strict::Bool=false, numeric_type::Symbol=:big)
    key = (Int(d), Int(Nmax), strict, numeric_type)
    lock(_projection_model_cache_lock) do
        models = get!(_projection_model_cache, coeffs) do
            Dict{Tuple{Int, Int, Bool, Symbol}, Any}()
        end
        return get!(models, key) do
            make_projection_model(d, coeffs; Nmax=Nmax, strict=strict, numeric_type=numeric_type)
        end
    end
end

function _cached_rho_ansatz(term::WeightedRhoTerm, rs, rt, ru)
    i = term.i
    j = term.j
    if term.n == 0 && term.m == 0
        return one(rs[i])
    elseif term.n == 1 && term.m == 0
        return (rs[i] + rt[i] + ru[i]) / 3
    elseif term.n == 0 && term.m == 1
        return (rs[j] + rt[j] + ru[j]) / 3
    elseif term.n == 1 && term.m == 1
        return (rs[i] * rt[j] + rs[i] * ru[j] + rt[i] * rs[j] +
            rt[i] * ru[j] + ru[i] * rs[j] + ru[i] * rt[j]) / 6
    end
    return nothing
end

function _cached_discT_rho_ansatz(term::WeightedRhoTerm, rs, ru, dt)
    i = term.i
    j = term.j
    if term.n == 0 && term.m == 0
        return zero(dt[i])
    elseif term.n == 1 && term.m == 0
        return dt[i] / 3
    elseif term.n == 0 && term.m == 1
        return dt[j] / 3
    elseif term.n == 1 && term.m == 1
        return (dt[i] * (rs[j] + ru[j]) + (rs[i] + ru[i]) * dt[j]) / 6
    end
    return nothing
end

_rho_with_root(s, sigma_root) = begin
    root4ms = _msqrt(4 - s)
    (sigma_root - root4ms) / (sigma_root + root4ms)
end

_disc_rho_with_root(s, sigma, sigma_root) =
    (2 * _msqrt(s - 4) * sigma_root) / (s + sigma - 8)

function ProjectionWorkspace(model::ProjectionModel, s)
    sample = _ascomplex(_rho_with_root(s, first(model.sigma_roots)))
    T = typeof(sample)
    n = length(model.sigmas)
    return ProjectionWorkspace(Vector{T}(undef, n), Vector{T}(undef, n),
        Vector{T}(undef, n), Vector{T}(undef, n))
end

function _fill_rho_s!(ws::ProjectionWorkspace, model::ProjectionModel, s)
    root4ms = _msqrt(4 - s)
    @inbounds for k in eachindex(model.sigmas)
        root = model.sigma_roots[k]
        ws.rs[k] = (root - root4ms) / (root + root4ms)
    end
    return ws
end

function _fill_rho_tu!(ws::ProjectionWorkspace, model::ProjectionModel, t, u)
    root4mt = _msqrt(4 - t)
    root4mu = _msqrt(4 - u)
    @inbounds for k in eachindex(model.sigmas)
        root = model.sigma_roots[k]
        ws.rt[k] = (root - root4mt) / (root + root4mt)
        ws.ru[k] = (root - root4mu) / (root + root4mu)
    end
    return ws
end

function _fill_discT_tu!(ws::ProjectionWorkspace, model::ProjectionModel, t, u)
    roottm4 = _msqrt(t - 4)
    root4mu = _msqrt(4 - u)
    @inbounds for k in eachindex(model.sigmas)
        sigma = model.sigmas[k]
        root = model.sigma_roots[k]
        ws.ru[k] = (root - root4mu) / (root + root4mu)
        ws.dt[k] = (2 * roottm4 * root) / (t + sigma - 8)
    end
    return ws
end

function exprSTU_prepared(model::ProjectionModel, s, t, u, ws::ProjectionWorkspace)
    total = zero(_ascomplex(s + t + u))
    if !isempty(model.terms)
        _fill_rho_tu!(ws, model, t, u)
        @inbounds for term in model.terms
            cached = _cached_rho_ansatz(term, ws.rs, ws.rt, ws.ru)
            total += term.weight * (cached === nothing ?
                rho_ansatz(s, t, u, term.n, term.m, term.sigma_i, term.sigma_j) : cached)
        end
    end

    if model.leading_threshold_weight !== nothing && !iszero(model.leading_threshold_weight)
        total += model.leading_threshold_weight * Delta(model.d, s, t, u)
    end
    @inbounds for term in model.subleading_threshold_terms
        total += term.weight * delta_Delta(model.d, s, t, u, term.p)
    end
    return total
end

function _threshold_s_contribution(model::ProjectionModel, s)
    total = zero(_ascomplex(s))
    if model.leading_threshold_weight !== nothing && !iszero(model.leading_threshold_weight)
        total += model.leading_threshold_weight * Delta(model.d, s)
    end
    @inbounds for term in model.subleading_threshold_terms
        total += term.weight * delta_Delta(model.d, s, term.p)
    end
    return total
end

function exprSTU_prepared_scache(model::ProjectionModel, s, t, u,
                                 ws::ProjectionWorkspace, threshold_s)
    total = zero(_ascomplex(s + t + u))
    if !isempty(model.terms)
        _fill_rho_tu!(ws, model, t, u)
        @inbounds for term in model.terms
            cached = _cached_rho_ansatz(term, ws.rs, ws.rt, ws.ru)
            total += term.weight * (cached === nothing ?
                rho_ansatz(s, t, u, term.n, term.m, term.sigma_i, term.sigma_j) : cached)
        end
    end

    total += threshold_s
    if model.leading_threshold_weight !== nothing && !iszero(model.leading_threshold_weight)
        total += model.leading_threshold_weight * (Delta(model.d, t) + Delta(model.d, u))
    end
    @inbounds for term in model.subleading_threshold_terms
        total += term.weight * (delta_Delta(model.d, t, term.p) + delta_Delta(model.d, u, term.p))
    end
    return total
end

function exprSTU(model::ProjectionModel, s, t, u, ws::ProjectionWorkspace)
    _fill_rho_s!(ws, model, s)
    return exprSTU_prepared(model, s, t, u, ws)
end

function exprSTU(model::ProjectionModel, s, t, u)
    ws = ProjectionWorkspace(model, s)
    return exprSTU(model, s, t, u, ws)
end

function discTexprSTU_prepared(model::ProjectionModel, s, t, u, ws::ProjectionWorkspace)
    total = zero(_ascomplex(s + t + u))
    isempty(model.terms) && return total

    _fill_discT_tu!(ws, model, t, u)
    @inbounds for term in model.terms
        cached = _cached_discT_rho_ansatz(term, ws.rs, ws.ru, ws.dt)
        total += term.weight * (cached === nothing ?
            discT_rho_ansatz(s, t, u, term.n, term.m, term.sigma_i, term.sigma_j) : cached)
    end
    return total
end

function discTexprSTU(model::ProjectionModel, s, t, u, ws::ProjectionWorkspace)
    _fill_rho_s!(ws, model, s)
    return discTexprSTU_prepared(model, s, t, u, ws)
end

function discTexprSTU(model::ProjectionModel, s, t, u)
    ws = ProjectionWorkspace(model, s)
    return discTexprSTU(model, s, t, u, ws)
end

function _expr_at_z_prepared(model::ProjectionModel, s, z, ws::ProjectionWorkspace)
    return exprSTU_prepared(model, s, tt(s, z), uu(s, z), ws)
end

function _expr_at_z_prepared(model::ProjectionModel, s, z, ws::ProjectionWorkspace, threshold_s)
    return exprSTU_prepared_scache(model, s, tt(s, z), uu(s, z), ws, threshold_s)
end

function _expr_at_z(model::ProjectionModel, s, z, ws::ProjectionWorkspace)
    _fill_rho_s!(ws, model, s)
    return _expr_at_z_prepared(model, s, z, ws)
end

function _discT_at_z_prepared(model::ProjectionModel, s, z, ws::ProjectionWorkspace)
    return discTexprSTU_prepared(model, s, tt(s, z), uu(s, z), ws)
end

expr(model::ProjectionModel, s, z) = exprSTU(model, s, tt(s, z), uu(s, z))
discTExpr(model::ProjectionModel, s, z) = discTexprSTU(model, s, tt(s, z), uu(s, z))

function _real_type_like(x)
    T = typeof(real(_ascomplex(x)))
    return T <: Integer ? Float64 : T
end

function _real_one_like(x)
    return one(_real_type_like(x))
end

function _real_const_like(x, y)
    return convert(_real_type_like(x), y)
end

function _eval_precision(f, model::ProjectionModel)
    if eltype(model.sigmas) <: BigFloat
        return setprecision(WP_BITS) do
            f()
        end
    end
    return f()
end

function PolQKernel(J, d::Integer)
    seed = _ascomplex(J + d)
    R = _real_type_like(seed)
    one_r = one(R)
    two_r = one_r + one_r
    four_r = two_r + two_r
    pi_r = convert(R, pi)
    half = inv(two_r)
    exponent = J + one_r
    a = (J + one_r) * half
    b = (J + two_r) * half
    c = J + (convert(R, d) - one_r) * half
    prefactor = (four_r / pi_r) *
        (sqrt(pi_r) * gamma(J + one_r) * gamma((convert(R, d) - two_r) * half)) /
        (two_r^exponent * gamma(c))
    return PolQKernel(prefactor, a, b, c, exponent)
end

function polQ(kernel::PolQKernel, z)
    zc = _ascomplex(z)
    return kernel.prefactor * zc^(-kernel.exponent) *
        HypergeometricFunctions._₂F₁(kernel.a, kernel.b, kernel.c, inv(zc^2))
end

function PgenKernel(J, d::Integer)
    R = _real_type_like(_ascomplex(J + d))
    return PgenKernel(-J, J + convert(R, d - 3), convert(R, (d - 2) / 2))
end

Pgen(kernel::PgenKernel, z) = HypergeometricFunctions._₂F₁(kernel.a, kernel.b, kernel.c, (1 - z) / 2)

function _polQ_series_radius(z0)
    zc = _ascomplex(z0)
    one_r = _real_one_like(zc)
    distances = [abs(zc), abs(zc - one_r), abs(zc + one_r)]
    finite_distances = filter(x -> isfinite(float(x)) && x > 0, distances)
    nearest = isempty(finite_distances) ? one_r : minimum(finite_distances)
    return min(_real_const_like(zc, 1e-4) * (one_r + abs(zc)), nearest / 10)
end

"Cauchy-coefficient implementation of Mathematica `SeriesCoefficient[PolQ[J,d,z], {z,z0,k}]`."
function polQ_series_coefficient(kernel::PolQKernel, z0, k::Integer; samples::Integer=64, radius=nothing)
    k < 0 && throw(ArgumentError("Taylor coefficient order must be non-negative"))
    k == 0 && return polQ(kernel, z0)

    r = radius === nothing ? _polQ_series_radius(z0) : radius
    total = zero(polQ(kernel, z0 + r))
    for q in 0:(samples - 1)
        theta = 2 * oftype(r, pi) * q / samples
        w = r * exp(complex(zero(r), theta))
        total += polQ(kernel, z0 + w) / w^k
    end
    return total / samples
end

polQ_series_coefficient(J, d::Integer, z0, k::Integer; samples::Integer=64, radius=nothing) =
    polQ_series_coefficient(PolQKernel(J, d), z0, k; samples=samples, radius=radius)

polQ_derivative_order(kernel::PolQKernel, z0, order::Integer) =
    order == 0 ? polQ(kernel, z0) : factorial(order) * polQ_series_coefficient(kernel, z0, order)

polQ_derivative_order(J, d::Integer, z0, order::Integer) =
    polQ_derivative_order(PolQKernel(J, d), z0, order)

function exprJmain1(qkernel::PolQKernel, model::ProjectionModel, s; p::Integer=P0)
    return _eval_precision(model) do
        srat = s
        s0 = z0t(srat)
        two_s0 = 2 * s0
        ws = ProjectionWorkspace(model, srat)
        _fill_rho_s!(ws, model, srat)
        pi_bound = _real_const_like(s0, pi)
        integrand(phi) = begin
            c = cos(phi)
            denom = 1 + c
            z = two_s0 / denom
            jac = two_s0 * sin(phi) / denom^2
            jac * polQ(qkernel, z) * _discT_at_z_prepared(model, srat, z, ws)
        end
        my_n_integrate_faster(integrand, 0, pi_bound, p)
    end
end

exprJmain1(J, model::ProjectionModel, s; p::Integer=P0) =
    exprJmain1(PolQKernel(J, model.d), model, s; p=p)

function exprJmain1(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
                    strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big))
    return exprJmain1(J, cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict,
        numeric_type=numeric_type), s; p=p)
end

function _divT_threshold(n::Integer, ss)
    alpha = n + _real_const_like(ss, 0.5)
    return (-1)^n * _mpow(ss - 4, -alpha)
end

function _disc_threshold_T(n::Integer, ss, z)
    return _divT_threshold(n, tt(ss, z))
end

function _keyhole_c0(n::Integer, ss, z0, dir)
    alpha = n + _real_const_like(ss, 0.5)
    lambda = _real_const_like(ss, 10.0)^(-max(30, min(WP ÷ 2, 80)))
    dz = lambda * dir
    return dz^alpha * _disc_threshold_T(n, ss, z0 + dz)
end

function exprJkeyhole(qkernel::PolQKernel, d::Integer, s, n::Integer; p::Integer=P0, ord::Integer=6)
    alpha = n + _real_const_like(s, 0.5)
    z0 = z0t(s)

    # Main integral displaced exactly as in the Mathematica keyhole code:
    # z0t[s] + Exp[I Arg[-1 + z0t[s]]] 10^-5.
    z0_one = _real_one_like(z0)
    main_dir = exp(complex(zero(z0_one), angle(-z0_one + z0)))
    shifted_s0 = z0 + main_dir * _real_const_like(z0, 1e-5)
    two_shifted_s0 = 2 * shifted_s0
    pi_bound = _real_const_like(z0, pi)
    main_integrand(phi) = begin
        c = cos(phi)
        denom = 1 + c
        z = two_shifted_s0 / denom
        jac = two_shifted_s0 * sin(phi) / denom^2
        jac * polQ(qkernel, z) * _disc_threshold_T(n, s, z)
    end
    main = my_n_integrate_faster(main_integrand, 0, pi_bound, p)

    # Local keyhole correction around z0t[s].
    dir = exp(complex(zero(z0_one), angle(z0 - z0_one)))
    eps = _real_const_like(z0, 1e-5)
    c0 = _keyhole_c0(n, s, z0, dir)
    keyhole = zero(main)
    for k in 0:ord
        qk = polQ_series_coefficient(qkernel, z0, k)
        keyhole += qk * dir^(k + 1 - alpha) * eps^(k + 1 - alpha) / (k + 1 - alpha)
    end

    return main + c0 * keyhole
end

exprJkeyhole(J, d::Integer, s, n::Integer; p::Integer=P0, ord::Integer=6) =
    exprJkeyhole(PolQKernel(J, d), d, s, n; p=p, ord=ord)

function exprJthresh(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
                     strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big))
    return exprJthresh(J, cached_projection_model(d, coeffs; Nmax=Nmax,
        strict=strict, numeric_type=numeric_type), s; p=p)
end

function exprJthresh(qkernel::PolQKernel, model::ProjectionModel, s; p::Integer=P0)
    return _eval_precision(model) do
        z0 = z0t(s)
        residue = zero(polQ(qkernel, z0))
        nmax_residue = floor(Int, (model.d - 3) / 2)
        pi_like = _real_const_like(s, pi)
        sminus4_half = (s - 4) / 2
        for n in nmax_residue:-1:1
            m_idx = Int(iseven(model.d)) - 2 * (n - nmax_residue)
            weight = get(model.threshold_weights, m_idx, nothing)
            weight === nothing && continue
            factor = (m_idx == 0) ? sin(pi_like * (model.d - 3) / 2) : 1
            residue += factor * weight * (pi_like * (-1)^(n + 1)) / factorial(n - 1) *
                polQ_derivative_order(qkernel, z0, n - 1) * _mpow(sminus4_half, -n)
        end

        keyhole = zero(residue)
        nmax_keyhole = floor(Int, (model.d - 4) / 2)
        for n in nmax_keyhole:-1:0
            m_idx = Int(isodd(model.d)) - 2 * (n - nmax_keyhole)
            weight = get(model.threshold_weights, m_idx, nothing)
            weight === nothing && continue
            factor = (m_idx == 0) ? sin(pi_like * (model.d - 3) / 2) : 1
            keyhole += factor * weight * exprJkeyhole(qkernel, model.d, s, n; p=p)
        end

        residue + keyhole
    end
end

exprJthresh(J, model::ProjectionModel, s; p::Integer=P0) =
    exprJthresh(PolQKernel(J, model.d), model, s; p=p)

"Froissart-Gribov projection used for generic complex spin J."
function exprJ_FG(J, model::ProjectionModel, s; p::Integer=P0)
    return _eval_precision(model) do
        qkernel = PolQKernel(J, model.d)
        exprJmain1(qkernel, model, s; p=p) + exprJthresh(qkernel, model, s; p=p)
    end
end

function exprJ_FG(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
                  strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big))
    return exprJ_FG(J, cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict,
        numeric_type=numeric_type), s; p=p)
end

function IPt(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
             strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big))
    return IPt(J, cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict,
        numeric_type=numeric_type), s; p=p)
end

#%% Usual Pj projection and optimized projection

Pgen(J, d::Integer, z) = Pgen(PgenKernel(J, d), z)

function IPt(J, model::ProjectionModel, s; p::Integer=P0)
    return _eval_precision(model) do
        pgen_kernel = PgenKernel(J, model.d)
        ws = ProjectionWorkspace(model, s)
        _fill_rho_s!(ws, model, s)
        threshold_s = _threshold_s_contribution(model, s)
        pi_bound = _real_const_like(s, pi)
        integrand(theta) = begin
            c = cos(theta)
            sin(theta) * _mpow(1 - c^2, (model.d - 4) / 2) * Pgen(pgen_kernel, c) *
                _expr_at_z_prepared(model, s, c, ws, threshold_s)
        end
        my_n_integrate_faster(integrand, 0, pi_bound, p)
    end
end

function _is_even_positive_integer(J)
    if J isa Integer
        return J > 0 && iseven(J)
    elseif J isa Real
        return J > 0 && isinteger(J) && iseven(Integer(J))
    else
        # Keep the complex-J plane on the Froissart-Gribov branch even at
        # points with zero imaginary part; pass a real integer to request IPt.
        return false
    end
end

function _as_big_complex(x)
    z = _ascomplex(x)
    return complex(BigFloat(real(z)), BigFloat(imag(z)))
end

"Final projection: use `IPt` on positive even real-integer J, otherwise the FG representation."
function exprJ(J, model::ProjectionModel, s; p::Integer=P0)
    use_ipt = _is_even_positive_integer(J)
    if eltype(model.sigmas) <: BigFloat
        return setprecision(WP_BITS) do
            JJ = _as_big_complex(J)
            ss = _as_big_complex(s)
            use_ipt ? IPt(JJ, model, ss; p=p) : exprJ_FG(JJ, model, ss; p=p)
        end
    end
    return use_ipt ? IPt(J, model, s; p=p) : exprJ_FG(J, model, s; p=p)
end

function exprJ(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
               strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big))
    model = cached_projection_model(d, coeffs; Nmax=Nmax, strict=strict, numeric_type=numeric_type)
    return exprJ(J, model, s; p=p)
end

plane_wave_projection(J, d::Integer, coeffs, s; Nmax::Integer=20, p::Integer=P0,
                      strict::Bool=false, numeric_type::Symbol=(p <= 4 ? :float64 : :big)) =
    exprJ(J, d, coeffs, s; Nmax=Nmax, p=p, strict=strict, numeric_type=numeric_type)