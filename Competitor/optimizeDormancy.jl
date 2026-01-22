using Statistics
import Optim
using DataFrames

const EPS = 1e-9 #to avoid div by zero

function objective(dormancy; phi::Float64, dt::Float64, growth::Float64,
                   simfn=simulateGrowth, total_days::Float64=1000.0, eval_window::Float64=300.0,
                   dilute_dormant::Bool=true, dilution_effect::Float64=dilEffect)
    active, _, _, _, _, _, _, time =
        simfn([growth, dormancy];
              phi=phi,
              withDormancy=true,
              dilution_interval=dt,
              total_time=total_days,
              run_partial_last=true,
              dilute_dormant=dilute_dormant,
              dilEffect=dilution_effect)

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

function _tail_mean_active_menace(phi::Float64, dt::Float64, q::Float64;
    growth::Float64,
    total_days::Float64=1000.0,
    eval_window::Float64=300.0,
    dilute_dormant::Bool=true,
    dilution_effect::Float64=dilEffect,
    simfn=simulateGrowth,
)
    A, _, _, _, _, M, _, t =
        simfn([growth, q];
              phi=phi,
              withDormancy=true,
              dilution_interval=dt,
              total_time=total_days,
              run_partial_last=true,
              dilute_dormant=dilute_dormant,
              dilEffect=dilution_effect)

    t_end = t[end]
    keep  = t .>= (t_end - eval_window)
    A_keep = filter(isfinite, collect(@view A[keep]))
    M_keep = filter(isfinite, collect(@view M[keep]))
    mA = isempty(A_keep) ? NaN : mean(A_keep)
    mM = isempty(M_keep) ? NaN : mean(M_keep)
    return mA, mM
end

function bracketMinimum(f, lo, hi; n::Int=25)
    xs = collect(range(lo, hi; length=n))
    ys = map(f, xs)
    k  = findmin(ys)[2]
    a  = (k == 1)            ? xs[1]       : xs[k-1]
    b  = (k == length(xs))   ? xs[end]     : xs[k+1]
    δ  = 1e-9*(hi - lo)
    a, b = max(lo, a + δ), min(hi, b - δ)
    return a, b, xs, ys, k
end

function bracketBrent(f, lo, hi; n::Int=25, reltol=1e-3, abstol=1e-6, maxiters=200)
    a, b, xs, ys, k = bracketMinimum(f, lo, hi; n=n)
    res = Optim.optimize(f, a, b, Optim.Brent();
                         rel_tol=reltol, abs_tol=abstol, iterations=maxiters)
    return res, xs, ys, k
end

function optimizeDormancy(phi::Float64, dt::Float64;
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0, 3.0),
    total_days::Float64=1000.0,
    eval_window::Float64=300.0,
    reltol::Real=1e-3,
    abstol::Real=1e-6,
    maxiters::Int=200,
    dilute_dormant::Bool=true,
    dilution_effect::Float64=dilEffect)

    lo, hi = float(bounds[1]), float(bounds[2])
    hi = min(hi, growth)  #dormancy can’t exceed growth

    f(d) = objective(d; phi=phi, dt=dt, growth=growth,
                     total_days=total_days, eval_window=eval_window,
                     dilute_dormant=dilute_dormant,
                     dilution_effect=dilution_effect)

    f_lo = f(lo)
    f_hi = f(hi)

    resB, xs, ys, k = bracketBrent(f, lo, hi; n=31,
                                   reltol=reltol, abstol=abstol, maxiters=maxiters)

    use_brent = Optim.converged(resB) && isfinite(Optim.minimum(resB))
    if use_brent
        qB   = Optim.minimizer(resB)
        fB   = Optim.minimum(resB)

        tol = max(reltol * abs(fB), abstol)
        qbest, fbest, chosen = qB, fB, "Brent(bracketed)"

        return (qbest, fbest, (; method=chosen, converged=true,
                               iterations=Optim.iterations(resB)))
    end

    ncoarse = 61
    coarse  = collect(range(lo, hi; length=ncoarse))  # includes endpoints
    vC      = map(f, coarse)
    for i in eachindex(vC); if !isfinite(vC[i]); vC[i] = 1e12; end; end
    jC      = argmin(vC)

    a = (jC == 1)            ? coarse[1]      : coarse[jC-1]
    b = (jC == length(coarse)) ? coarse[end]   : coarse[jC+1]

    nfine = 21
    fine  = collect(range(a, b; length=nfine))
    vF    = map(f, fine)
    for i in eachindex(vF); if !isfinite(vF[i]); vF[i] = 1e12; end; end
    jF    = argmin(vF)
    qG    = fine[jF]
    fG    = vF[jF]

    tolG = max(reltol * abs(fG), abstol)
    qbest, fbest, chosen = qG, fG, "grid_adaptive"


    return (qbest, fbest, (; method=chosen, converged=isfinite(fbest),
                           iterations=ncoarse + nfine))
end

#classic heatmap
function heatMap(phiRange::Tuple, dtRange::Tuple;
    nphi::Int=21, ndt::Int=21,
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0,5.0),
    total_days::Float64=1000.0, eval_window::Float64=500.0,
    threaded::Bool=true,
    reltol::Real=1e-3, abstol::Real=1e-6, maxiters::Int=200,
    dilute_dormant::Bool=false, dilution_effect::Float64=dilEffect
)
    phis = exp10.(range(log10(float(phiRange[1])), log10(float(phiRange[2])); length=nphi))
    dts  = collect(range(float(dtRange[1]),  float(dtRange[2]);  length=ndt))

# results with two extra columns
    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[], loss=Float64[],
                        avg_active=Float64[], avg_menace=Float64[])

    lk = ReentrantLock()

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:nphi
            for j in 1:ndt
                dopt, fmin, _ = optimizeDormancy(phis[i], dts[j];
                    growth=growth, bounds=bounds,
                    total_days=total_days, eval_window=eval_window,
                    reltol=reltol, abstol=abstol, maxiters=maxiters,
                    dilute_dormant=dilute_dormant,
                    dilution_effect=dilution_effect)
                mA, mM = _tail_mean_active_menace(phis[i], dts[j], dopt;
                    growth=growth, total_days=total_days, eval_window=eval_window,
                    dilute_dormant=dilute_dormant, dilution_effect=dilution_effect)
                Base.lock(lk) do
                    push!(results, (phis[i], dts[j], dopt, fmin, mA, mM))
                end
            end
        end
    else
        for i in 1:nphi, j in 1:ndt
            dopt, fmin, _ = optimizeDormancy(phis[i], dts[j];
                growth=growth, bounds=bounds,
                total_days=total_days, eval_window=eval_window,
                reltol=reltol, abstol=abstol, maxiters=maxiters,
                dilute_dormant=dilute_dormant,
                dilution_effect=dilution_effect)
            mA, mM = _tail_mean_active_menace(phis[i], dts[j], dopt;
                growth=growth, total_days=total_days, eval_window=eval_window,
                dilute_dormant=dilute_dormant, dilution_effect=dilution_effect)
            push!(results, (phis[i], dts[j], dopt, fmin, mA, mM))
        end
    end

    return results
end


using DataFrames

#track type of solver used
function heatMapTracked(phiRange::Tuple, dtRange::Tuple;
    nphi::Int=21, ndt::Int=21,
    growth::Float64,
    bounds::Tuple{<:Real,<:Real}=(0.0,5.0),
    total_days::Float64=1000.0, eval_window::Float64=500.0,
    threaded::Bool=true, reltol::Real=1e-3, abstol::Real=1e-6, maxiters::Int=200,
    dilute_dormant::Bool=false, dilution_effect::Float64=dilEffect
)
    phis = exp10.(range(log10(float(phiRange[1])), log10(float(phiRange[2])); length=nphi))
    dts  = collect(range(float(dtRange[1]),  float(dtRange[2]);  length=ndt))

    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[], loss=Float64[],
        method=String[], iters=Int[], avg_active=Float64[], avg_menace=Float64[])

    mylock = ReentrantLock()

    work(i,j) = begin
        dopt, fmin, st = optimizeDormancy(phis[i], dts[j];
            growth=growth, bounds=bounds,
            total_days=total_days, eval_window=eval_window,
            reltol=reltol, abstol=abstol, maxiters=maxiters,
            dilute_dormant=dilute_dormant,
            dilution_effect=dilution_effect)
        mA, mM = _tail_mean_active_menace(phis[i], dts[j], dopt;
            growth=growth, total_days=total_days, eval_window=eval_window,
            dilute_dormant=dilute_dormant, dilution_effect=dilution_effect)
        lock(mylock) do
            push!(results, (phis[i], dts[j], dopt, fmin, String(st.method), Int(st.iterations), mA, mM))
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

#for dt / dileffect variation
function heatMapDilution(dtRange::Tuple;
    neff::Int=11, ndt::Int=12,
    phi::Float64,
    growth::Float64,
    effRange::Tuple{<:Real,<:Real}=(0.20, 1.00),
    bounds::Tuple{<:Real,<:Real}=(0.0, 3.0),
    total_days::Float64=1000.0,
    eval_window::Float64=500.0,
    threaded::Bool=true,
    reltol::Real=1e-3, abstol::Real=1e-6, maxiters::Int=200,
    dilute_dormant::Bool=false)

    effs = collect(range(float(effRange[1]), float(effRange[2]); length=neff))
    dts  = collect(range(float(dtRange[1]),  float(dtRange[2]);  length=ndt))

    bins = [Vector{NTuple{4,Float64}}() for _ in 1:Threads.nthreads()]

    work(e, d, tid) = begin
        dopt, fmin, _ = optimizeDormancy(phi, d;
            growth=growth, bounds=bounds,
            total_days=total_days, eval_window=eval_window,
            reltol=reltol, abstol=abstol, maxiters=maxiters,
            dilute_dormant=dilute_dormant,
            dilution_effect=e)
        push!(bins[tid], (e, d, dopt, fmin))
    end
    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:neff
            tid = Threads.threadid()
            for j in 1:ndt
                work(effs[i], dts[j], tid)
            end
        end
    else
        for i in 1:neff, j in 1:ndt
            work(effs[i], dts[j], 1)
        end
    end

    flat = reduce(vcat, bins)
    mAs = Float64[]; mMs = Float64[]
    for (e, d, dopt, _) in flat
        mA, mM = _tail_mean_active_menace(phi, d, dopt;
            growth=growth, total_days=total_days, eval_window=eval_window,
            dilute_dormant=dilute_dormant, dilution_effect=e)
        push!(mAs, mA); push!(mMs, mM)
    end
    return DataFrame(dilEffect = getindex.(flat, 1),
                     dt        = getindex.(flat, 2),
                     q         = getindex.(flat, 3),
                     loss      = getindex.(flat, 4),
                     avg_active = mAs,
                     avg_menace = mMs)
end
