using Statistics
using Optim: optimize, Brent, Options
using DataFrames

const EPS = 1e-9

"""
    objective(dormancy; phi, dt, growth, simfn=simulateGrowth,
              total_days=500.0, eval_window=300.0)

- `growth` is known max growth (optParameters[1])
- `dormancy` is optParameters[2]
- `phi` forwards to model
- `dt` is  `dilution_interval`
"""
function objective(dormancy; phi::Float64, dt::Float64, growth::Float64,
                   simfn=simulateGrowth, total_days::Float64=500.0, eval_window::Float64=300.0)
    active, _, _, _, _, time =
        simfn([growth, dormancy];
              phi=phi,
              withDormancy=true,
              dilution_interval=dt,
              total_time=total_days,
              run_partial_last=true)

    t_end = time[end]
    keep = time .>= (t_end - eval_window)
    tail = @view active[keep]

    vals = filter(isfinite, collect(tail))
    if isempty(vals)
        return Inf
    end
    m = mean(vals)
    return (isfinite(m) && m > 0) ? (1.0 / (m + EPS)) : Inf
end

"""
adding in an intermediary between Brent and the fallback grid, where it tries to bracket around the minimum
"""

function bracketMinimum(f, lo, hi; n::Int=25)
    xs = collect(range(lo, hi; length=n))
    ys = map(f, xs)
    k  = findmin(ys)[2]
    a  = (k == 1)            ? xs[1]       : xs[k-1]
    b  = (k == length(xs))   ? xs[end]     : xs[k+1]
    δ  = 1e-9*(hi - lo)                      #nudge to avoid equal endpoints
    a, b = max(lo, a + δ), min(hi, b - δ)
    return a, b, xs, ys, k
end

function bracketBrent(f, lo, hi; n::Int=25, reltol=1e-3, abstol=1e-6, maxiters=200)
    a, b, xs, ys, k = bracketMinimum(f, lo, hi; n=n)
    res = Optim.optimize(f, a, b, Optim.Brent();
                         rel_tol=reltol, abs_tol=abstol, iterations=maxiters)
    return res, xs, ys, k
end

"""
    optimizeDormancy(phi, dt; growth, bounds=(0.0,1.0),
                     total_days=500.0, eval_window=300.0,
                     reltol=1e-3, abstol=1e-6, maxiters=200)
"""

function optimizeDormancy(phi::Float64, dt::Float64;
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0, 3.0),
    total_days::Float64=500.0,
    eval_window::Float64=300.0,
    reltol::Real=1e-3,
    abstol::Real=1e-6,
    maxiters::Int=200,
)
    lo, hi = float(bounds[1]), float(bounds[2])

    f(d) = objective(d; phi=phi, dt=dt, growth=growth,
                        total_days=total_days, eval_window=eval_window)

    result = optimize(f, lo, hi, Brent();
        rel_tol=reltol, abs_tol=abstol, iterations=maxiters,
        show_trace=false, store_trace=false)
  
    dopt = Optim.minimizer(result)
    fmin = Optim.minimum(result)
  
    if !Optim.converged(result) || !isfinite(fmin)
        #bracket brent
        res2, xs, ys, k = bracketBrent(f, lo, hi; n=31,
                                              reltol=reltol, abstol=abstol, maxiters=maxiters)
        if Optim.converged(res2) && isfinite(Optim.minimum(res2))
            return (Optim.minimizer(res2), Optim.minimum(res2),
                    (; method="Brent(bracketed)", converged=true,
                       iterations=Optim.iterations(res2)))
        end
    
        #final safety: adaptive grid (coarse → zoom around best cell)
        a = (k == 1) ? xs[1] : xs[k-1]
        b = (k == length(xs)) ? xs[end] : xs[k+1]
        fine = collect(range(a, b; length=25))
        v2   = map(f, fine)
        j    = findmin(v2)[2]
        return (fine[j], v2[j],
                (; method="adaptive_grid", converged=isfinite(v2[j]),
                   iterations=length(xs) + length(fine)))
    end
    
    return (dopt, fmin, (; method="Brent", converged=true, iterations=Optim.iterations(result)))
end

"""
    heatMap(phiRange::Tuple, dtRange::Tuple; nphi=21, ndt=21, growth,
            bounds=(0.0,1.0), total_days=500.0, eval_window=300.0,
            threaded=true, reltol=1e-3, abstol=1e-6, maxiters=200)
"""
function heatMap(phiRange::Tuple, dtRange::Tuple;
    nphi::Int=21, ndt::Int=21,
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0,5.0),
    total_days::Float64=500.0, eval_window::Float64=300.0,
    threaded::Bool=true,
    reltol::Real=1e-3, abstol::Real=1e-6, maxiters::Int=200,
)
    phis = exp10.(range(log10(float(phiRange[1])), log10(float(phiRange[2])); length=nphi))
    dts  = collect(range(float(dtRange[1]),  float(dtRange[2]);  length=ndt))

    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[], loss=Float64[])

    lock = ReentrantLock()

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:nphi
            for j in 1:ndt
                dopt, fmin, _ = optimizeDormancy(phis[i], dts[j];
                    growth=growth, bounds=bounds,
                    total_days=total_days, eval_window=eval_window,
                    reltol=reltol, abstol=abstol, maxiters=maxiters)
                lock(lock) do
                    push!(results, (phis[i], dts[j], dopt, fmin))
                end
            end
        end
    else
        for i in 1:nphi, j in 1:ndt
            dopt, fmin, _ = optimizeDormancy(phis[i], dts[j];
                growth=growth, bounds=bounds,
                total_days=total_days, eval_window=eval_window,
                reltol=reltol, abstol=abstol, maxiters=maxiters)
            push!(results, (phis[i], dts[j], dopt, fmin))
        end
    end

    return results
end


using DataFrames

"""
    heatMapTracked(phiRange, dtRange; ...) -> (df, summary)
"""
function heatMapTracked(phiRange::Tuple, dtRange::Tuple;
    nphi::Int=21, ndt::Int=21,
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0,5.0),
    total_days::Float64=500.0, eval_window::Float64=300.0,
    threaded::Bool=true, reltol::Real=1e-3, abstol::Real=1e-6, maxiters::Int=200,
)
    phis = exp10.(range(log10(float(phiRange[1])), log10(float(phiRange[2])); length=nphi))
    dts  = collect(range(float(dtRange[1]),  float(dtRange[2]);  length=ndt))

    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[], loss=Float64[],
    method=String[], iters=Int[])

    mylock = ReentrantLock()

    work(i,j) = begin
        dopt, fmin, st = optimizeDormancy(phis[i], dts[j];
            growth=growth, bounds=bounds,
            total_days=total_days, eval_window=eval_window,
            reltol=reltol, abstol=abstol, maxiters=maxiters)
        lock(mylock) do
            push!(results, (phis[i], dts[j], dopt, fmin, String(st.method), Int(st.iterations)))
        end
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:nphi
            for j in 1:ndt
                work(i,j)
            end
        end
    else
        for i in 1:nphi, j in 1:ndt
            work(i,j)
        end
    end

    summary = combine(groupby(results, :method), nrow => :count)
    return results, summary
end
