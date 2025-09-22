include("growthSimulator.jl")
include("oldObjecFunc.jl")

phiRange = (1.0e-10, 1.0e-6)
dtRange  = (0.25, 21.0)
growth   = 0.69

results = heatMap(phiRange, dtRange; nphi=10, ndt=8, growth=growth)

using CSV
CSV.write("heatmap_results.csv", results)





df, summary = heatMapTracked((1e-10, 1e-6), (0.25, 21.0);
                             nphi=9, ndt=9, growth=0.69)

@show summary

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
