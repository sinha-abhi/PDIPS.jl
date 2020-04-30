"""
    reformulate(lp::Problem) where T

Convert from IplpProblem to IplpStandardProblem.
"""
function reformulate(lp::Problem) where T
    if isa(lp, StandardProblem)
        return lp
    end

    # count free and upper bounded variables
    free = 0
    ub = 0
    nzv = 0 # non-zero values
    for (i, (l, h)) in enumerate(zip(lp.lo, lp.hi))
        if l == -Inf && h == Inf
            # free
            free += 1
            nzv += length(lp.cols[i].ind)
        elseif isfinite(l) && isfinite(h)
            # l ⩽ x ⩽ h
            ub += 1
        elseif l == -Inf && isfinite(h)
            # x ⩽ h
            # will be dealt with later
        elseif isfinite(l) && h == Inf
            # l ⩽ x
            # will be dealt with later
        else
            error("unexpected bounds ($l, $h)")
        end

        nzv += length(lp.cols[i].ind)
    end

    nv = lp.nv

    I = Vector{Int}(undef, nzv) # row indices
    J = Vector{Int}(undef, nzv) # column indices
    V = Vector{T}(undef, nzv)
    ind_ub = Vector{Int}(undef, ub)
    val_ub = Vector{T}(undef, ub)

    b = Vector{T}(undef, lp.nc)
    c = Vector{T}(undef, nv + free)

    free = 0
    up = 0
    nzv = 0
    for (j, (l, h)) in enumerate(zip(lp.lo, lp.hi))
        column = lp.cols[i] # current column

        if l == -Inf && h == Inf
            # free variable
            # split into positive and negative parts: x = x⁺ - x⁻
            c[j + free]  = lp.c[i]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end

            c[j + free + 1] = -lp.c[i]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free + 1
                V[nzv] = -v
            end

            free += 1
        elseif isfinite(l) && isfinite(h)
            # l ⩽ x ⩽ h
            c[j + free] = lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                b[i] -= (v * l)

                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end

            up += 1
            ind_ub[ub] = j + free
            val_ub[ub] = h - l
        elseif l == -Inf && isfinite(h)
            # x ⩽ h
            c[j + free] = -lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                b[i] -= (-v * u)

                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = -v
            end
        elseif isfinite(l) && h == Inf
            # l ⩽ x
            c[j + free] = lp.c[j]
            for (i, v) in zip(column.ind, column.nzval)
                nzv += 1
                I[nzv] = i
                J[nzv] = j + free
                V[nzv] = v
            end
        else
            # this error _should_ have already been caught
            error("unexpected bounds ($l, $h)")
        end
    end

    # construct A out of I, J, V
    # FIXME: currently assuming that `A` will be sparse
    nv += free
    A = sparse(I, J, V, lp.nc, nv)

    IplpStandardProblem{T}(lp.nc, nv, ub, ind_ub, val_ub, A, b, c)
end

"""
    starting_pt!(iter::Iterate{T}) where T

Convert `iter` to starting point. It is faster to first allocate memory and then
fill the desired values.
The starting point is:
    (x, λ, s, τ, κ) = (e, 0, e, 1, 1).
"""
function starting_pt!(iter::Iterate{T}) where T
    # primal
    iter.x .= oneunit(T)
    iter.v .= oneunit(T)

    # dual
    iter.λ .= zero(T)
    iter.s .= oneunit(T)
    iter.w .= oneunit(T)

    # homogenous constants
    iter.τ = oneunit(T)
    iter.κ = oneunit(T)

    nothing
end

function duality!(i::Iterate{T}) where T
    iter.μ = (
        (i.τ * i.κ
        + dot(i.x, i.s)
        + dot(i.v, i.w))
        / (i.nv + i.nu + 1)
    )

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
    solv::IplpSolver{T},
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

    solv.status_primal = (_p <= tols.εp) ? IterateFeasible : IterateUndefined
    solv.status_dual   = (_d <= tols.εd) ? IterateFeasible : IterateUndefined

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

function get_solution(solv::IplpSolver{T}, org::Problem) where T
    sln = IplpSolution{T}()
    sln.flag = solv.status

    # TODO: recover solution to original problem from std form
    # sln.x =

    slp = solv.lp
    sln.A_std = slp.A
    sln.b_std = slp.b
    sln.c_std = slp.c

    f_iter = solv.iter
    sln.x_std = f_iter.x
    sln.λ_std = f_iter.λ
    sln.s_std = f_iter.s

    return sln
end
