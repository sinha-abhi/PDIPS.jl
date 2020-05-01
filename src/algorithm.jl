# TODO: finish implementation of steps
function ipm_step!(solv::Solver{T}) where T
    iter = solv.iter
    lp   = solv.lp
    res  = solv.res

    nc  = lp.nc # number of constraints
    nv  = lp.nv # number of variables
    nu  = lp.nu # number of upper bounded variables
    A   = lp.A
    b   = lp.b
    c   = lp.c
    ubi = lp.upb_ind
    ubv = lp.upb_val

    # scaling for problem data
    θxs = iter.s ./ iter.x
    θvw = iter.w ./ iter.v
    @views θxs[ubi] .+= θvw

    # regularization
    r_min = sqrt(eps(T)) # approx 1e-8
    solv.regP .= max.(r_min, solv.regP ./ 10)
    solv.regD .= max.(r_min, solv.regD ./ 10)
    solv.regG  = max(r_min, solv.regG / 10)

    # factorization
    # make three attempts, after increasing regularization
    # see Tulip.jl for details on factorization
    attempt = 0
    while attempt < 4
        try
            update_linear_solver!(solv.ls, θxs, solv.regP, solv.regD)
            break
        catch e
            if !(isa(e, PosDefException) || isa(e, ZeroPivotException))
                # unexpected error
                rethrow(e)
            end

            attempt += 1
            solv.regP .*= 100
            solv.regD .*= 100
            solv.regG  *= 100
        end
    end

    # factorization may not have been complete
    if attempt >= 3
        throw(PosDefException(0))
    end

    # predictor search direction
    Δ  = Iterate{T}(nc, nv, nu)
    # Δc = Iterate{T}(nc, nv, nu)

    δx = zeros(T, nv)
    δy = zeros(T, nc)
    δz = zeros(T, nu)
    solve_augsys!(δx, δy, δz, solv.ls, θvw, ubi, b, c, ubv)
    δ0 = solv.regG + iter.κ / iter.τ - dot(c, δx) + dot(b, δy) - dot(ubv, δz)

    # affine-scaling search direction
    solve_newton_system!(
        Δ, solv.ls, θvw, b, c, ubi, ubv, δx, δy, δz, δ0, iter,
        res.rp, res.ru, res.rd, res.rg,
        - iter.x .* iter.s, # ξxs
        - iter.v .* iter.w, # ξvw
        - iter.τ  * iter.κ  # ξτκ
    )

    # step length for affine-scaling direction
    # γ = (1 - α)^2 * min(β, 1 - α)
    # η = 1 - γ
    α = max_alpha(iter, Δ)
    γ = (oneunit(T) - α)^2 * min(T(1e-1), oneunit(T) - α)
    η = 1 - γ

    # apply mehrotra corrector
    solve_newton_system!(
        Δ, solv.ls, θvw, b, c, ubi, ubv, δx, δy, δz, δ0, iter,
        η .* res.rp, η .* res.ru, η .* res.rd, η .* res.rg,
        - iter.x .* iter.s .+ γ * iter.μ .- Δ.x .* Δ.s,
        - iter.v .* iter.w .+ γ * iter.μ .- Δ.v .* Δ.w,
        - iter.τ  * iter.κ  + γ * iter.μ  - Δ.τ  * Δ.κ
    )

    # calculate new step length, and dampen alpha
    α = max_alpha(iter, Δ)
    α *= T(9995 / 10_000)

    # update iterate
    update_iterate!(iter, Δ, α)

    nothing
end

function solve!(solv::Solver{T}, maxiter::Int) where T
    slp = solv.lp

    while true
        # update residuals
        update_residuals!(
            solv.res, solv.iter,
            slp.A, slp.b, slp.c,
            slp.upb_ind, slp.upb_val
        )

        # update duality measure
        duality!(solv.iter)

        # check solver status -- exit if needed
        check_status!(
            solv, solv.iter, solv.res, solv.tols,
            slp.A, slp.b, slp.c, slp.upb_ind, slp.upb_val
        )
        if (
            solv.status == SolverOptimal          ||
            solv.status == SolverPrimalInfeasible ||
            solv.status == SolverDualInfeasible
        )
            break
        elseif solv.niter >= maxiter
            solv.status = SolverExceededIterations
            break
        end

        # take next ipm_step
        try
            ipm_step!(solv)
        catch e
            if isa(e, OutOfMemoryError)
                solv.status = SolverExceededMemory
            elseif isa(e, PosDefException) || isa(e, SingularException)
                solv.status = SolverNumericalInstability
            else
                solv.status = SolverUndefined
                rethrow(e)
            end

            break
        end

        solv.niter += 1
    end

    nothing
end
