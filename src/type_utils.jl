"""
    reformulate(lp::Problem) where T

Convert from IplpProblem to IplpStandardProblem.
"""
function reformulate(lp::AbstractProblem)
    if isa(lp, AbstractStandardProblem)
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

    StandardProblem{T}(lp.nc, nv, ub, ind_ub, val_ub, A, b, c)
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

function get_solution(solv::Solver{T}, org::AbstractProblem) where T
    sln = Solution{T}()
    sln.status = solv.status

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