using Random, BlackBoxOptim, Plots, CSV, DataFrames

# Constraint check (returns true if point is in a forbidden area)
function is_forbidden(x, y)
    # Circle at center
    if (x - 5)^2 + (y - 5)^2 < 2^2
        return true
    end
    # Lower-left quarter circle
    if x ≤ 1 && y ≤ 1 && x^2 + y^2 < 1^2
        return true
    end
    # Upper-left quarter circle
    if x ≤ 1 && y ≥ 9 && x^2 + (y - 10)^2 < 1^2
        return true
    end
    return false
end

# Objective function: negative of the minimum pairwise distance
function min_dist_obj(vec::Vector{Float64})
    N = length(vec) ÷ 2
    coords = [vec[2i-1:2i] for i in 1:N]

    # If any point is forbidden, penalize
    if any(p -> is_forbidden(p[1], p[2]), coords)
        return 1e6
    end

    # Compute all pairwise distances
    dists = [sqrt(sum((coords[i] .- coords[j]).^2)) for i in 1:N for j in i+1:N]
    return -minimum(dists)  # maximize min dist = minimize -min dist
end

results = DataFrame(n=Int[], mindist=Float64[])

for n in 2:10
    println("Optimizing layout for n = $n...")
    res = bboptimize(min_dist_obj;
        SearchRange = (0.0, 10.0),
        NumDimensions = 2n,
        MaxSteps = 1000000,
        Method = :adaptive_de_rand_1_bin_radiuslimited,
        TraceMode = :silent,
        AbsoluteTolerance = 1e-5
    )

    best_layout = best_candidate(res)
    best_mindist = -best_fitness(res)
    push!(results, (n, best_mindist))

    if n in 2:10
        coords = [best_layout[2i-1:2i] for i in 1:n]
        scatter(first.(coords), last.(coords),
            xlim=(0,10), ylim=(0,10),
            xlabel="x", ylabel="y", label="Turbines",
            title="Turbine Layout for n = $n", aspect_ratio=:equal)

        # Full circle at (5,5)
        θ = range(0, 2π, length=100)
        plot!(5 .+ 2*cos.(θ), 5 .+ 2*sin.(θ), lc=:red, label="Forbidden Zone")

        # Quarter circle at (0,0) — bottom-left
        θ_bl = range(0, π/2, length=100)
        plot!(1*cos.(θ_bl), 1*sin.(θ_bl), lc=:red, label="")

        # Quarter circle at (0,10) — top-left
        θ_tl = range(3π/2, 2π, length=100)
        plot!(1*cos.(θ_tl), 10 .+ 1*sin.(θ_tl), lc=:red, label="")

        savefig("layout_optim_delivs/layout_n$n.png")
    end
end

# Save results
CSV.write("layout_optim_delivs/wind_farm_results.csv", results; writeheader=false)

# Plot 1: Line plot of min distance vs n
plot(results.n, results.mindist,
    xlabel="Number of Turbines", ylabel="Min Pairwise Distance",
    title="Minimum Distance vs Number of Turbines", lw=2, marker=:circle)
savefig("layout_optim_delivs/mindist_vs_n.png")

# Plot 2: Efficiency plot
eff = 1 ./(1 .+ 1 ./results.mindist)
scatter(results.mindist, eff, xlabel="Minimum Distance", ylabel="Efficiency", title="Efficiency vs Min Distance",
        legend=false, marker=:circle, ms=10)

# Annotate each point with its n value
for (i, n) in enumerate(results.n)
    annotate!(results.mindist[i], eff[i], text(string(n), 0, 0, 8))
end
savefig("layout_optim_delivs/efficiency_scatter.png")