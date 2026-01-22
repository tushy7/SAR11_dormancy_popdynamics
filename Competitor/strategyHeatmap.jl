include("simFunctions.jl")       # simulateGrowth + modelParams
include("optimizeDormancy.jl")   # optimization functions
using Statistics
using Plots
using Optim

# Your varibles here
phi_min, phi_max = 1e-10, 1e-6
nphi = 12  # number of phi samples (increase for higher resolution)

dt_min, dt_max = 5.0, 30.0
ndt = 12  # number of dt samples (increase for higher resolution)

dilution_effect = 0.1
umax = 0.69  # max growth rate

resus_stoch = 0.10  # "stochastic dormancy"
resus_spore = 0.01  # "spore-forming dormancy"

total_days = 1000.0
eval_window = 500.0  # last X days to average over
dilute_dormant = false

# End user defined variables

phi_vals = exp10.(range(log10(phi_min), log10(phi_max); length=nphi))
dt_vals = range(dt_min, dt_max; length=ndt)

mp_stoch = copy(modelParams)
mp_stoch[:resuscitationRate] = resus_stoch

mp_spore = copy(modelParams)
mp_spore[:resuscitationRate] = resus_spore


const EPS = 1e-9

function tail_mean(x, t; window::Float64)
    t_end = t[end]
    keep = t .>= (t_end - window)
    vals = filter(isfinite, collect(@view x[keep]))
    return isempty(vals) ? NaN : mean(vals)
end

function objective_with_params(dormancy; phi::Float64, dt::Float64, growth::Float64,
                               mp::Dict,
                               total_days::Float64,
                               eval_window::Float64,
                               dilute_dormant::Bool,
                               dilution_effect::Float64,
                               isSporeFormer::Bool=false)

    A, _, _, _, _, _, time =
        simulateGrowth([growth, dormancy];
            phi=phi,
            withDormancy=true,
            dilution_interval=dt,
            total_time=total_days,
            run_partial_last=true,
            dilute_dormant=dilute_dormant,
            dilEffect=dilution_effect,
            customModelParams=mp,
            isSporeFormer=isSporeFormer)

    t_end = time[end]
    keep = time .>= (t_end - eval_window)
    tail = filter(isfinite, collect(@view A[keep]))

    if isempty(tail)
        return Inf
    end

    m = mean(tail)
    return (isfinite(m) && m > 0) ? 1.0 / (m + EPS) : Inf
end

function optimizeDormancy_params(phi::Float64, dt::Float64;
    growth::Float64,
    mp::Dict,
    bounds::Tuple{<:Real,<:Real}=(0.0, 3.0),
    total_days::Float64=total_days,
    eval_window::Float64=eval_window,
    dilute_dormant::Bool=dilute_dormant,
    dilution_effect::Float64=dilution_effect,
    isSporeFormer::Bool=false,
    reltol::Real=1e-3,
    abstol::Real=1e-6,
    maxiters::Int=200)

    lo, hi = float(bounds[1]), float(bounds[2])
    hi = min(hi, growth)

    f(d) = objective_with_params(d;
        phi=phi, dt=dt, growth=growth, mp=mp,
        total_days=total_days, eval_window=eval_window,
        dilute_dormant=dilute_dormant, dilution_effect=dilution_effect,
        isSporeFormer=isSporeFormer)

    res = Optim.optimize(f, lo, hi, Optim.Brent();
                         rel_tol=reltol, abs_tol=abstol,
                         iterations=maxiters)

    return Optim.minimizer(res), Optim.minimum(res)
end

function tail_active_menace_params(phi::Float64, dt::Float64, q::Float64;
    growth::Float64,
    mp::Dict,
    total_days::Float64=total_days,
    eval_window::Float64=eval_window,
    dilute_dormant::Bool=dilute_dormant,
    dilution_effect::Float64=dilution_effect,
    isSporeFormer::Bool=false)

    A, Q, _, _, I, M, _, t =
        simulateGrowth([growth, q];
            phi=phi,
            withDormancy=true,
            dilution_interval=dt,
            total_time=total_days,
            run_partial_last=true,
            dilute_dormant=dilute_dormant,
            dilEffect=dilution_effect,
            customModelParams=mp,
            isSporeFormer=isSporeFormer)

    mA = tail_mean(A, t; window=eval_window)
    mQ = tail_mean(Q, t; window=eval_window)
    mI = tail_mean(I, t; window=eval_window)
    mM = tail_mean(M, t; window=eval_window)

    mT = mA + mQ + mI
    return mT, mM
end

println("Computing strategy heatmap...")
println("Phi range: $phi_min to $phi_max ($nphi points)")
println("Dt range: $dt_min to $dt_max ($ndt points)")
println("Total simulations: $(nphi * ndt * 3) (this may take a while)")

winner_matrix = zeros(Int, ndt, nphi)  # rows = dt, cols = phi
cell_counts_no_dorm = zeros(ndt, nphi)
cell_counts_stoch = zeros(ndt, nphi)
cell_counts_spore = zeros(ndt, nphi)

for (i, dt) in enumerate(dt_vals)
    for (j, φ) in enumerate(phi_vals)
        
        total_sims = nphi * ndt
        current_sim = (i-1) * nphi + j
        if current_sim % 10 == 0
            println("Progress: $current_sim / $total_sims")
        end
        
        q_stoch_temp, _ = optimizeDormancy_params(φ, dt;
            growth=umax, mp=mp_stoch, isSporeFormer=false)
        _, mM = tail_active_menace_params(φ, dt, q_stoch_temp;
            growth=umax, mp=mp_stoch, isSporeFormer=false)
        cell_counts_no_dorm[i, j] = mM
        
        q_stoch, _ = optimizeDormancy_params(φ, dt;
            growth=umax, mp=mp_stoch, isSporeFormer=false)
        mT_stoch, _ = tail_active_menace_params(φ, dt, q_stoch;
            growth=umax, mp=mp_stoch, isSporeFormer=false)
        cell_counts_stoch[i, j] = mT_stoch
        
        q_spore, _ = optimizeDormancy_params(φ, dt;
            growth=umax, mp=mp_spore, isSporeFormer=true)
        mT_spore, _ = tail_active_menace_params(φ, dt, q_spore;
            growth=umax, mp=mp_spore, isSporeFormer=true)
        cell_counts_spore[i, j] = mT_spore
        
        max_cells = max(mM, mT_stoch, mT_spore)
        if mM == max_cells
            winner_matrix[i, j] = 0  # no dormancy
        elseif mT_stoch == max_cells
            winner_matrix[i, j] = 1  # stochastic
        else
            winner_matrix[i, j] = 2  # spore-forming
        end
    end
end

println("Simulation complete!")

using Colors
colors = [colorant"#1b9e77", colorant"#d95f02", colorant"#7570b3"]
cmap = cgrad(colors, 3, categorical=true)

plt = heatmap(
    phi_vals, dt_vals, winner_matrix;
    xlabel="Viral stress (phi)",
    ylabel="Dilution interval (days)",
    title="Winning dormancy strategy",
    xscale=:log10,
    color=cmap,
    colorbar=false,  # Turn off the colorbar
    size=(900, 600),
    clims=(0, 2)  # Force color limits to cover all three strategies
)

legend_x = log10(phi_max) * 0.7  # Adjust position as needed
legend_y_start = dt_max * 0.85
legend_spacing = (dt_max - dt_min) * 0.08

annotate!(plt, legend_x, legend_y_start, 
    text("No dormancy", 10, :left, colorant"#1b9e77"))
annotate!(plt, legend_x, legend_y_start - legend_spacing, 
    text("Stochastic", 10, :left, colorant"#d95f02"))
annotate!(plt, legend_x, legend_y_start - 2*legend_spacing, 
    text("Spore-forming", 10, :left, colorant"#7570b3"))

display(plt)

savefig(plt, "strategy_heatmap.png")
println("Heatmap saved as strategy_heatmap.png")

total_points = nphi * ndt
n_no_dorm = count(==(0), winner_matrix)
n_stoch = count(==(1), winner_matrix)
n_spore = count(==(2), winner_matrix)

println("\n--- Summary Statistics ---")
println("Total parameter combinations: $total_points")
println("No dormancy wins: $n_no_dorm ($(round(100*n_no_dorm/total_points, digits=1))%)")
println("Stochastic dormancy wins: $n_stoch ($(round(100*n_stoch/total_points, digits=1))%)")
println("Spore-forming dormancy wins: $n_spore ($(round(100*n_spore/total_points, digits=1))%)")