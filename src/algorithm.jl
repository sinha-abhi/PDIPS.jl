# TODO: implement steps
function ipm_step!(solv::IplpSolver{T}) where T

end

function solve!(solv::IplpSolver{T}, maxiter::Int) where T
    slp = solv.slp

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
        catch err
            # TODO: catch errors
        end

        solv.niter += 1
    end

    nothing
end
