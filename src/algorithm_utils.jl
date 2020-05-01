function update_iterate!(iter::Iterate{T}, Δ::Iterate{T}, α::T) where T
    iter.x .+= α .* Δ.x
    iter.v .+= α .* Δ.v
    iter.λ .+= α .* Δ.λ
    iter.s .+= α .* Δ.s
    iter.w .+= α .* Δ.w

    iter.τ += α * Δ.τ
    iter.κ += α * Δ.κ

    nothing
end

function solve_newton_system!(
    Δ::Iterate{T}, # overwritten during solve
    ls::AbstractLinearSolver{T},
    θvw::Vector{T},
    b::Vector{T},
    c::Vector{T},
    ubi::Vector{Int},
    uval::Vector{T},
    δx::Vector{T},
    δy::Vector{T},
    δz::Vector{T},
    δ0::T,
    iter::Iterate{T},
    ξp::Vector{T},
    ξu::Vector{T},
    ξd::Vector{T},
    ξg::T,
    ξxs::Vector{T},
    ξvw::Vector{T},
    ξτκ::T
) where T
    _ξd   = copy(ξd)
    _ξu   = copy(ξu)
    _ξd .-= (ξxs ./ iter.x)
    _ξu .-= (ξvw ./ iter.w)

    solve_augsys!(Δ.x, Δ.λ, Δ.w, ls, θvw, ubi, ξp, _ξd, _ξu)

    Δ.τ = (ξg + (ξτκ / iter.τ) + dot(c, Δ.x) - dot(b, Δ.λ) + dot(ubv, Δ.w)) / δ0
    Δ.κ = (ξτκ - iter.κ  * Δ.τ) / iter.τ

    Δ.x .+= (Δ.τ .* δx)
    Δ.λ .+= (Δ.τ .* δy)
    Δ.w .+= (Δ.τ .* δw)

    # δs = inv(X) * (ξxs - Sδx)
    Δ.s  .= ξxs
    Δ.s .-= (iter.s .* Δ.x)
    Δ.s ./= iter.x

    # δw = inv(W) * (ξvw - v .* Δx)
    Δ.v  .= ξvw
    Δ.v .-= (iter.v .* Δ.w)
    Δ.v ./= iter.w

    nothing
end

function solve_augsys!(
    δx::Vector{T}, # overwritten during solve
    δy::Vector{T}, # overwritten during solve
    δz::Vector{T}, # overwritten during solve
    ls::AbstractLinearSolver{T},
    θvw::Vector{T},
    ubi::Vector{Int},
    ξp::Vector{T},
    ξd::Vector{T},
    ξu::Vector{T}
) where T
    δz .= zero(T)

    _ξd = copy(ξd)
    @views _ξd[ubi] .-= (ξu .* θvw)

    solve_augmented_system!(δx, δy, ls, ξp, _ξd)

    δz .-= ξu
    @views δz .+= dx[ubi]
    δz .*= θvw

    nothing
end

function update_residuals!(
    res::Residuals{T},
    iter::Iterate{T},
    A::AbstractMatrix{T},
    b::Vector{T},
    c::Vector{T},
    ubi::Vector{T},
    ubv::Vector{T}
) where T
    # calculate `rp` and its norm
    mul!(res.rp, A, iter.x)    # rp = A * x
    rmul!(res.rp, -oneunit(T)) # rp = - rp = - A * x
    axpy!(iter.τ, b, res.rp)   # rp = b * τ - rp
    res.rpn = norm(res.rp)

    # calculate `ru` and its norm
    rmul!(res.ru, zero(T))                       # zero out ru
    axpy!(-oneunit(T), iter.v, res.ru)           # ru = -w + ru = -w
    @views axpy!(-oneunit(T), pt.x[ubi], res.ru) # ru = ru - x
    axpy!(iter.τ, ubv, res.ru)                   # ru = τ * u + ru
    res.run = norm(res.ru)

    # calculate `rd` and its norm
    mul!(res.rd, transpose(A), iter.λ)           # rd = A * λ
    rmul!(res.rd, -oneunit(T))                   # rd = - rd = - A * λ
    axpy!(iter.τ, c, res.rd)                     # rd = τ * c + rd
    axpy!(-oneunit(T), iter.s, res.rd)           # rd = rd - s
    @views axpy!(oneunit(T), pt.w, res.rd[uind]) # rd = w + rd

    # calculate `rg` and its norm
    # rg = c'x - (b'λ + u'w) + κ
    res.rg = iter.κ + dot(c, iter.x) - dot(b, iter.λ) + dot(ubv, iter.w)
    res.rgn = norm(rg)

    nothing
end

function check_status!(
    solv::Solver{T},
    iter::Iterate{T},
    res::Residuals{T},
    tols::Tolerances{T},
    A::AbstractMatrix{T},
    b::Vector{T},
    c::Vector{T},
    ubi::Vector{T},
    ubv::Vector{T}
) where T
    # _p = max(|rp| / (τ * (1 + |b|)), |ru| / (τ * (1 + |u|)))
    _p = max(
        res.rpn / (iter.τ * (oneunit(T) + norm(b, Inf)),
        res.run / (iter.τ * (oneunit(T) + norm(ubv, Inf))))
    )

    # _d = |rd| / (τ * (1 + |c|))
    _d = res.rdn / (iter.τ * (oneunit(T) + norm(c, Inf)))

    primal_bound = dot(c, iter.x)                     # c'x
    dual_bound   = dot(b, iter.λ) - dot(ubv, iter.w)  # b'λ - u'w

    # _g = |c'x - b'λ| / (τ + |b'λ|)
    _g = abs(primal_bound - dual_bound) / (iter.τ + abs(dual_bound))

    solv.status_primal = _p <= tols.εp ? IterateFeasible : IterateUndefined
    solv.status_dual   = _d <= tols.εd ? IterateFeasible : IterateUndefined

    # check optimality
    if _p <= tols.εp && _d <= tols.εd && _g <= εg
        solv.status        = SolverOptimal
        solv.status_primal = IterateOptimal
        solv.status_dual   = IterateOptimal

        return nothing
    end

    # check infeasibility
    _ax_n = max(
        norm(A * iter.x, Inf),
        norm(iter.x[ubi] + iter.v, Inf)
    )
    _cb_n = norm(c, Inf) / max(1, norm(b, Inf))
    if _ax_n * _cb_n < - tols.εi * dot(c, iter.x)
        solv.status_primal = IterateInfeasible
        solv.status        = SolverDualInfeasible

        return nothing
    end

    δ = transpose(A) * iter.λ + iter.s
    δ[ubi] .-= iter.w
    if norm(δ, Inf) * norm(b, Inf) / max(1, norm(c, Inf)) < tols.εi * dual_bound
        solv.status_dual = IterateInfeasible
        solv.status      = SolverPrimalInfeasible

        return nothing
    end

    nothing
end

function _max_alpha(v::Vector{T}, dv::Vector{T}) where T
    n = size(x, 1)
    size(dv, 1) == n || throw(DimensionMismatch(
        "expected vector of length $n"
    ))

    α = T(Inf)
    for i in Base.OneTo(n)
        if dv[i] <  zero(T)
            _pot_α = -v[i] / dv[i]
            if _pot_α < α
                α = _pot_α
            end
        end
    end

    α
end

function max_alpha(iter::Iterate{T}, Δ::Iterate{T}) where T
    α_τ = Δ.τ < zero(T) ? (-iter.τ / Δ.τ) : oneunit(T)
    α_κ = Δ.κ < zero(T) ? (-iter.κ / Δ.κ) : oneunit(T)

    α = min(
        oneunit(T),
        _max_alpha(iter.x, Δ.x),
        _max_alpha(iter.v, Δ.v),
        _max_alpha(iter.s, Δ.s),
        _max_alpha(iter.w, Δ.w),
        α_τ, α_κ
    )
end
