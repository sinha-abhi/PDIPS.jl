# TODO: implement algorithm

function ipm_step!()

end

function solve!(solv::IplpSolver{T}, maxiter::Int) where T
    while true
        # update residuals
        # update duality measure

        # check solver status -- exit if needed
        # take next ipm_step

        solv.iter += 1
    end
end

"""
    solve(lp::Problem, maxiter::Int = 100, tol::T) where T

Solve the linear program. Returns an `IplpSolution`.
"""
function solve(lp::Problem, maxiter::Int = 100, tol::T = 1e-8) where T <: Real
    # convert problem to standard form
    slp = reformulate(lp)

    # starting point
    iter = Iterate{T}(slp.nc, slp.nv, slp.nu)
    duality!(iter)

    tols = Tolerances{T}(tol)
    solv = IplpSolver{T}(slp, iter, tol)

    # call internal HSD algorithm
    solve!(solv, maxiter)

    # extract solution from solver
    get_solution(solv, lp)
end
