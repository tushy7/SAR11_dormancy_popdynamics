include("simFunctions.jl")
include("optimizeDormancy.jl")

phiRange = (1.0e-10, 1.0e-6)
dtRange  = (0.25, 21.0)
umax   = 0.69

results = heatMap(phiRange, dtRange;
                  nphi=2, ndt=2, growth=umax,
                  dilute_dormant=false,
                  dilution_effect=dilEffect)

using CSV
CSV.write("heatmap_results.csv", results)





df, summary = heatMapTracked((1e-10, 1e-6), (0.25, 21.0);
                             nphi=9, ndt=9, growth=umax, dilute_dormant = true, dilution_effect=dilEffect)

@show summary

using Plots
# pick one that used the adaptive grid (or any specific cell)
row = first(filter(:method => ==("Brent"), df))
plt, info = visualize_optimization(row.phi, row.dt; growth=growth, bounds=(0.0, umax))
display(plt)
@show info

# counts by method
combine(groupby(df, :method), nrow => :count)

# fraction by dt or phi
combine(groupby(df, :dt), :method => (m->count(==("adaptive_grid"), m)) => :grid,
                           nrow => :N) |>
    x -> transform(x, [:grid, :N] => ByRow(/) => :frac_grid)

# plot a “method map” (0 = Brent, 1 = grid)
using Statistics
using DataFrames
pivot = unstack(transform(df, :method => ByRow(m->m=="adaptive_grid" ? 1 : 0) => :gridflag),
                :dt, :phi, :gridflag)  # table to heatmap

#dense sweep & inspect both metrics
phi, dt = 1e-7, 5.0
qs   = range(0.0, 2.5; length=101)
loss = [objective(q; phi=phi, dt=dt, growth=0.69) for q in qs]

#compute the tail mean directly for context
function tail_mean(q)
    act,_,_,_,_, t = simulateGrowth([0.69, q]; phi=phi, dilution_interval=dt,
                                    total_time=500.0, run_partial_last=true)
    t_end = t[end]; keep = t .>= (t_end - 300.0)
    vals = filter(isfinite, collect(act[keep]))
    isempty(vals) ? NaN : mean(vals)
end
mact = [tail_mean(q) for q in qs]

#print a few sentinel points to sanity check
for (q, L, m) in zip(qs[1:10:end], loss[1:10:end], mact[1:10:end])
    @info "q=$(round(q, digits=2))" loss=L mean_tail=m
end

scanDormancyFromCSV("heatmap_results.csv", 3, 4)




include("simFunctions.jl")
include("optimizeDormancy.jl")

phi     = 0.0
umax    = 0.69
dtRange = (0.5, 21.0)

df_qstar = heatMapDilution(dtRange;
    neff=8, ndt=8,
    phi=phi, growth=umax,
    effRange=(0.0, 0.40),
    dilute_dormant=false)   # false = dilute ACTIVE only

using CSV
CSV.write("qstar_dilution_map_$(phi).csv", df_qstar)



dtRange = (0.25, 21.0)
ndt     = 8

# φ×Δt map
df_phi_dt = heatMap((1e-10, 1e-7), dtRange;
    nphi=2, ndt=ndt, growth=umax,
    dilute_dormant=false,           # or true — but the same in both calls
    dilution_effect=0.25,
    total_days=1000.0, eval_window=300.0)

# Dilution map (single ϕ row)
df_dil = heatMapDilution(dtRange;
    neff=1, ndt=ndt,
    phi=1e-9, growth=umax,
    effRange=(0.25, 0.25),          # just the 0.25 row
    dilute_dormant=false,           # same as above
    total_days=1000.0, eval_window=300.0)

using DataFrames, Statistics

left  = sort(filter(:phi => ==(1e-7), df_phi_dt), :dt)[:, [:dt, :q]]
right = sort(filter(:dilEffect => ==(0.25), df_dil), :dt)[:, [:dt, :q]]
cmp   = innerjoin(left, right; on=:dt, makeunique=true)
rename!(cmp, [:dt, :q_phi_dt, :q_dil025])
cmp.diff = cmp.q_phi_dt .- cmp.q_dil025
println(cmp)
@show maximum(abs, cmp.diff)
    