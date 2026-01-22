include("simFunctions.jl")       # simulateGrowth + modelParams
include("optimizeDormancy.jl")   # bracketBrent etc.
using Statistics
using Plots

# User input here. Toggle either :phi or :dt
x_axis_variable = :dt  # Change this to :dt to vary dilution interval instead

chosen_dt = 5.0                          # fixed dilution interval (days)
phi_min, phi_max = 1e-8, 1e-6           # stress range (phi)
nphi = 100                                # number of phi samples

chosen_phi = 1e-6                         # fixed viral stress
dt_min, dt_max = 5.0, 30.0               # dilution interval range (days)
ndt = 100                                 # number of dt samples

dilution_effect = 0.1                    # dilEffect
umax = 0.69                              # max growth rate (opt parameter 1)

resus_stoch = 0.10                        # "stochastic dormancy"
resus_spore = 0.01                        # "spore-forming dormancy"

total_days = 1000.0
eval_window = 500.0                      # last X days to average over
dilute_dormant = false

# End user input variables

if x_axis_variable == :phi
    x_vals = exp10.(range(log10(phi_min), log10(phi_max); length=nphi))
    x_label = "Stress (phi)"
    x_scale = :log10
elseif x_axis_variable == :dt
    x_vals = range(dt_min, dt_max; length=ndt)
    x_label = "Dilution interval (days)"
    x_scale = :identity
else
    error("x_axis_variable must be :phi or :dt")
end

mp_stoch = copy(modelParams)
mp_stoch[:resuscitationRate] = resus_stoch

mp_spore = copy(modelParams)
mp_spore[:resuscitationRate] = resus_spore


const EPS = 1e-9  

function tail_mean(x, t; window::Float64)
    t_end = t[end]
    keep  = t .>= (t_end - window)
    vals  = filter(isfinite, collect(@view x[keep]))
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
            phi              = phi,
            withDormancy     = true,
            dilution_interval= dt,
            total_time       = total_days,
            run_partial_last = true,
            dilute_dormant   = dilute_dormant,
            dilEffect        = dilution_effect,
            customModelParams= mp,
            isSporeFormer    = isSporeFormer)

    t_end = time[end]
    keep  = time .>= (t_end - eval_window)
    tail  = filter(isfinite, collect(@view A[keep]))

    if isempty(tail)
        return Inf
    end

    m = mean(tail)
    return (isfinite(m) && m > 0) ? 1.0 / (m + EPS) : Inf
end

function optimizeDormancy_params(phi::Float64, dt::Float64;
    growth::Float64,
    mp::Dict,
    bounds::Tuple{<:Real,<:Real} = (0.0, 3.0),
    total_days::Float64 = total_days,
    eval_window::Float64 = eval_window,
    dilute_dormant::Bool = dilute_dormant,
    dilution_effect::Float64 = dilution_effect,
    isSporeFormer::Bool = false,
    reltol::Real = 1e-3,
    abstol::Real = 1e-6,
    maxiters::Int = 200)

    lo, hi = float(bounds[1]), float(bounds[2])
    hi = min(hi, growth)  # dormancy can't exceed growth

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
    total_days::Float64 = total_days,
    eval_window::Float64 = eval_window,
    dilute_dormant::Bool = dilute_dormant,
    dilution_effect::Float64 = dilution_effect,
    isSporeFormer::Bool = false)

    A, Q, _, _, I, M, _, t =
        simulateGrowth([growth, q];
            phi              = phi,
            withDormancy     = true,
            dilution_interval= dt,
            total_time       = total_days,
            run_partial_last = true,
            dilute_dormant   = dilute_dormant,
            dilEffect        = dilution_effect,
            customModelParams= mp,
            isSporeFormer    = isSporeFormer)

    mA = tail_mean(A, t; window=eval_window)
    mQ = tail_mean(Q, t; window=eval_window)
    mI = tail_mean(I, t; window=eval_window)
    mM = tail_mean(M, t; window=eval_window)

    mT = mA + mQ + mI
    return mT, mM
end

no_dormancy = similar(x_vals)  # menace = no-dormancy potential
stochastic  = similar(x_vals)  # optimal stochastic dormancy
spore_form  = similar(x_vals)  # optimal spore-formers

for (i, x) in pairs(x_vals)
    if x_axis_variable == :phi
        φ = x
        Δt = chosen_dt
    else  # :dt
        φ = chosen_phi
        Δt = x
    end
    
    if i % 10 == 0
        println("Progress: $i / $(length(x_vals))")
    end
    
    q_stoch, _ = optimizeDormancy_params(φ, Δt;
        growth=umax, mp=mp_stoch)

    mA_stoch, mM = tail_active_menace_params(φ, Δt, q_stoch;
        growth=umax, mp=mp_stoch)

    stochastic[i]  = mA_stoch
    no_dormancy[i] = mM

    q_spore, _ = optimizeDormancy_params(φ, Δt;
        growth=umax, mp=mp_spore, isSporeFormer=true)

    mA_spore, _ = tail_active_menace_params(φ, Δt, q_spore;
        growth=umax, mp=mp_spore, isSporeFormer=true)

    spore_form[i] = mA_spore
end

# Making the plots

if x_axis_variable == :phi
    plot_title = "Dormancy strategies vs viral stress (dt = $chosen_dt days)"
else
    plot_title = "Dormancy strategies vs dilution interval (phi = $chosen_phi)"
end

plt = plot(
    x_vals, no_dormancy;
    lw      = 2,
    label   = "No dormancy (menace)",
    xlabel  = x_label,
    ylabel  = "Total cells",
    xscale  = x_scale,
    title   = plot_title,
    legend  = :topright,
    grid    = false,
)

plot!(plt, x_vals, stochastic; lw=2, label="Stochastic dormancy")
plot!(plt, x_vals, spore_form; lw=2, ls=:dash, label="Spore-forming dormancy")

display(plt)

if x_axis_variable == :phi
    savefig(plt, "strategy_curves_vs_phi.png")
    println("Plot saved as strategy_curves_vs_phi.png")
else
    savefig(plt, "strategy_curves_vs_dt.png")
    println("Plot saved as strategy_curves_vs_dt.png")
end