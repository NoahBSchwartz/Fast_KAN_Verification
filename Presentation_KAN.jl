using LinearAlgebra
using Plots

struct BSpline
    grid_knots::Vector{Float64}
    ctrl_points::Vector{Float64}
    degree::Int
end

function BSpline(lb::Float64, ub::Float64, degree::Int; n::Int = 20, rng::AbstractRNG=Random.GLOBAL_RNG)
    k = degree
    cpoints::Vector{Float64} = randn(rng, Float64, n+k)  # Use provided RNG
    grid_pts::Vector{Float64} = LinRange(lb, ub, n+1)
    
    # Add k points before and after the bounds
    for i in 1:k
        pushfirst!(grid_pts, lb - (i * (ub-lb)/(2n)))
        push!(grid_pts, ub + (i * (ub-lb)/(2n)))
    end
    
    BSpline(grid_pts, cpoints, k)
end

function basis(bs::BSpline, x::Float64)
    n = size(bs.grid_knots)[1] - 1
    k = bs.degree
    # Initialize basis functions array for all degrees
    b = Vector{Vector{Float64}}(undef, k+1)
    # zeroth basis function (piecewise constant)
    b[1] = (x .>= bs.grid_knots[1:end-1]) .& (x .< bs.grid_knots[2:end])
    
    # compute basis functions for degrees 1 to k
    for degree in 1:k
        b[degree+1] = zeros(Float64, n-degree)
        for i in 1:(n-degree)
            numerator_left = x - bs.grid_knots[i]
            numerator_right = bs.grid_knots[i+degree+1] - x
            denominator_left = bs.grid_knots[i+degree] - bs.grid_knots[i]
            denominator_right = bs.grid_knots[i+degree+1] - bs.grid_knots[i+1]
            b[degree+1][i] = (numerator_left ./ denominator_left) .* b[degree][i] + (numerator_right ./ denominator_right) .* b[degree][i+1]
        end
    end
    
    return b
end

function eval_spline(bs::BSpline, x)
    b = basis(bs, x)[bs.degree+1]
    b ⋅ bs.ctrl_points
end

bs = BSpline(-10.0, 10.0, 3)
x = LinRange(-11.0, 11.0, 1000)
y = [eval_spline(bs, x_elt) for x_elt in x]
plot(x, y, label="spline")
scatter!(bs.grid_knots[bs.degree+1:end-bs.degree-1], bs.ctrl_points, label="control points")

function basis_derivative(bs::BSpline, x::Float64)
    n = size(bs.grid_knots)[1] - 1
    k = bs.degree
    b = basis(bs, x)
    # Initialize derivative arrays for all degrees
    d = Vector{Vector{Float64}}(undef, k+1)
    # derivative of degree 0 is zero everywhere
    d[1] = zeros(Float64, n)

    for degree in 1:k
        d[degree+1] = zeros(Float64, n-degree)
        for i in 1:(n-degree)
            term1 = b[degree][i] / (bs.grid_knots[i+degree] - bs.grid_knots[i])
            term2 = b[degree][i+1] / (bs.grid_knots[i+degree+1] - bs.grid_knots[i+1])
            d[degree+1][i] = degree * (term1 - term2)
        end
    end
    return d
end

function eval_derivative(bs::BSpline, x::Float64)
    d = basis_derivative(bs, x)[bs.degree+1]
    return d ⋅ bs.ctrl_points
end

function plot_spline_with_derivative(bs::BSpline, x)
    y = [eval_spline(bs, x_elt) for x_elt in x]
    dy = [eval_derivative(bs, x_elt) for x_elt in x]
    p1 = plot(x, y, label="spline", title="B-spline")
    scatter!(p1, bs.grid_knots[bs.degree+1:end-bs.degree-1], bs.ctrl_points, label="control points")
    p2 = plot(x, dy, label="derivative", title="B-spline derivative")
    plot(p1, p2, layout=(2,1))
end

function find_critical_points(bs::BSpline, tolerance::Float64=1e-6)
    x = LinRange(-20.0, 20.0, 10000)
    y_deriv = [eval_derivative(bs, x_i) for x_i in x]
    sign_changes = Float64[]
    for i in 1:(length(y_deriv)-1)
        if sign(y_deriv[i+1]) != sign(y_deriv[i])
            push!(sign_changes, x[i+1])
        end
    end
    return sign_changes
end

critical_points = find_critical_points(bs)

x = LinRange(-20.0, 20.0, 1000)
y = [eval_spline(bs, x_elt) for x_elt in x]
crit_points = find_critical_points(bs)
crit_values = [eval_spline(bs, x) for x in crit_points]
p = plot(x, y, label="spline", legend=:outertopright)
scatter!(crit_points, crit_values, color=:red, 
        markersize=6, label="critical points")
y_deriv = [eval_derivative(bs, x_elt) for x_elt in x]
plot!(x, y_deriv, label="derivative", linestyle=:dash)

function interval_analysis_spline(bs::BSpline, x_min::Float64, x_max::Float64)
    x = LinRange(-20.0, 20.0, 1000)
    y = [eval_spline(bs, x_elt) for x_elt in x]
    crit_points = find_critical_points(bs)
    crit_values = Float64[]
    for x in crit_points
          if x > x_min && x < x_max
                push!(crit_values, eval_spline(bs, x))
          end
    end
    max_y, min_y = maximum(crit_values), minimum(crit_values)
    return max_y, min_y
end

x_min, x_max = -5.0, 5.0
max_y, min_y = interval_analysis_spline(bs, x_min, x_max)
p = plot(x, y, label="spline", legend=:outertopright)
plot!([x_min, x_max], [max_y, max_y], color=:green, 
    label="maximum", linewidth=2, linestyle=:dash)
plot!([x_min, x_max], [min_y, min_y], color=:purple, 
    label="minimum", linewidth=2, linestyle=:dash)
return p

struct KAN
    inner_splines::Vector{BSpline}
    outer_spline::BSpline
end

function KAN(n_inner::Int=3; seed::Int=42)
    rng = Random.MersenneTwister(seed)
    inner_splines = [BSpline(-10.0, 10.0, 3, rng=rng) for _ in 1:n_inner]
    outer_spline = BSpline(-10.0, 10.0, 3, rng=rng)
    KAN(inner_splines, outer_spline)
end

function (kan::KAN)(x::Float64)
    # Evaluate KAN at point x (Sum inner values and evaluate through outer spline)
    inner_values = [eval_spline(spline, x) for spline in kan.inner_splines]
    inner_sum = sum(inner_values)
    return eval_spline(kan.outer_spline, inner_sum)
end

x_range = LinRange(-10, 10, 1000)
kan = KAN(3) # Make a KAN with 3 inner splines
kan_values = [kan(x) for x in x_range]
p = plot(layout=(5,1), size=(800,800))
for (i, spline) in enumerate(kan.inner_splines)
    inner_values = [eval_spline(spline, x) for x in x_range]
    plot!(p[i], x_range, inner_values, 
          label="Inner Spline $i", 
          title="Inner Spline $i",
          ylabel="y")
end
outer_values = [eval_spline(kan.outer_spline, x) for x in x_range]
plot!(p[4], x_range, outer_values, 
      label="Outer Spline", 
      title="Outer Spline",
      ylabel="y")
plot!(p[5], x_range, kan_values, 
      label="KAN Output", 
      title="Complete KAN",
      xlabel="x",
      ylabel="y")
display(p)

function kan_interval_analysis_vanilla(x_range, kan)
    p1 = plot(layout=(5,1), size=(800,800))
    inner_intervals = Vector{Tuple{Float64, Float64}}()
    for (i, spline) in enumerate(kan.inner_splines)
          width = rand() * 4 + 0.5
          placer = 7 - rand() * 14
          x_min, x_max = placer - 0.5 * width, placer + 0.5 * width
          max_y, min_y = interval_analysis_spline(spline, x_min, x_max)
          push!(inner_intervals, (min_y, max_y))
          inner_values = [eval_spline(spline, x) for x in x_range]
          plot!(p1[i], x_range, inner_values,
                label="Inner Spline $i",
                title="Inner Spline $i",
                ylabel="y")
          plot!(p1[i], [x_min, x_max], [max_y, max_y],
                color=:green,
                label="Maximum",
                linewidth=2,
                linestyle=:dash)
          plot!(p1[i], [x_min, x_max], [min_y, min_y],
                color=:purple,
                label="Minimum",
                linewidth=2,
                linestyle=:dash)
          vline!(p1[i], [x_min, x_max], 
                color=:gray, 
                linestyle=:dot, 
                label=false)
    end

    sum_min = sum(interval[1] for interval in inner_intervals)
    sum_max = sum(interval[2] for interval in inner_intervals)
    outer_max, outer_min = interval_analysis_spline(kan.outer_spline, sum_min, sum_max)
    inner_values = [eval_spline(kan.outer_spline, x) for x in x_range]
    plot!(p1[4], x_range, inner_values, 
          label="Outer Spline", 
          title="Outer Spline",
          xlabel="x",
          ylabel="y")
    plot!(p1[4], [sum_min, sum_max], [outer_max, outer_max],
            color=:green,
            label="Maximum",
            linewidth=2,
            linestyle=:dash)
    plot!(p1[4], [sum_min, sum_max], [outer_min, outer_min],
            color=:purple,
            label="Minimum",
            linewidth=2,
            linestyle=:dash)
    vline!(p1[4], [sum_min, sum_max], 
          color=:gray, 
          linestyle=:dot, 
          label=false)

      kan_vals = [kan(x) for x in x_range]
      plot!(p1[5], x_range, kan_vals, 
            label="KAN Output", 
            linewidth=2)
      hline!(p1[5], [outer_max, outer_min], 
            color=:green, 
            linestyle=:dash, 
            label="Guaranteed Output Interval")
      title!(p1[5], "Complete KAN Output")
      xlabel!(p1[5], "x")
      ylabel!(p1[5], "Output")
    return p1, (outer_min, outer_max)
end

x_range = LinRange(-10, 10, 1000)
p, (outer_min, outer_max) = kan_interval_analysis_vanilla(x_range, kan)
println("Output is guaranteed to be in [$(outer_min), $(outer_max)]")
display(p)

function fit_line_through_points(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    return slope, intercept
end

function compute_max_error(x_points::Vector{Float64}, y_points::Vector{Float64}, slope::Float64, intercept::Float64)
    y_predicted = slope .* x_points .+ intercept
    return maximum(abs.(y_predicted - y_points))
end

struct DynamicSplineApproximation
    segments::Vector{Tuple{Float64, Float64, Float64, Float64}}  # (start_x, end_x, slope, intercept)
    max_error::Float64
end

function compute_segment_error(x_points::Vector{Float64}, y_points::Vector{Float64}, start_idx::Int, end_idx::Int)
    if end_idx <= start_idx
        return Inf
    end
    slope, intercept = fit_line_through_points(x_points[start_idx], y_points[start_idx], x_points[end_idx], y_points[end_idx])
    return compute_max_error(x_points[start_idx:end_idx], y_points[start_idx:end_idx], slope, intercept)
end

function find_bspline_segments_given_max_segments(bs::BSpline, max_segments::Int, min_x::Float64, max_x::Float64)
    n_samples = 200
    x_points = collect(LinRange(min_x, max_x, n_samples))
    y_points = [eval_spline(bs, x) for x in x_points]
    
    # Initialize dynamic programming tables
    # dp[i,j] = min error using i segments up to point j
    dp = fill(Inf, (max_segments, n_samples))
    # back[i,j] = endpoint of previous segment in the optimal solution
    back = zeros(Int, (max_segments, n_samples))
    
    # Base case: first segment
    for j in 2:n_samples
        error = compute_segment_error(x_points, y_points, 1, j)
        dp[1,j] = error
    end
    
    # Other Segments
    for i in 2:max_segments
        for j in i:n_samples
            # Try all possible previous segment endpoints and write the maximum to the slot
            for k in (i-1):(j-1)
                prev_error = dp[i-1,k]
                curr_error = compute_segment_error(x_points, y_points, k, j)
                total_error = max(prev_error, curr_error)
                
                if total_error < dp[i,j]
                    dp[i,j] = total_error
                    back[i,j] = k
                end
            end
        end
    end
    
    segments = Tuple{Float64, Float64, Float64, Float64}[]
    curr_seg = n_samples
    
    # Backtrack through the segments
    for i in max_segments:-1:1
        prev_seg = i > 1 ? back[i,curr_seg] : 1
        
        # Get segment endpoints
        x1, y1 = x_points[prev_seg], y_points[prev_seg]
        x2, y2 = x_points[curr_seg], y_points[curr_seg]
        
        # Compute line parameters
        slope, intercept = fit_line_through_points(x1, y1, x2, y2)
        
        # Add segment to result
        pushfirst!(segments, (x1, x2, slope, intercept))
        
        curr_seg = prev_seg
    end
    
    return segments, dp[max_segments,n_samples]  # Return segments and achieved max error
end

function fit_bspline_segments_given_max_segments(kan::KAN, max_segments::Int)
    min_x = minimum(kan.inner_splines[1].grid_knots)
    max_x = maximum(kan.inner_splines[1].grid_knots)

    inner_approximations = DynamicSplineApproximation[]
    max_achieved_error = 0.0
    for spline in kan.inner_splines
        segments, error = find_bspline_segments_given_max_segments(spline, max_segments, min_x, max_x)
        max_achieved_error = max(max_achieved_error, error)
        push!(inner_approximations, DynamicSplineApproximation(segments, max_error))
    end
    
    min_x = minimum(kan.outer_spline.grid_knots)
    max_x = maximum(kan.outer_spline.grid_knots)
    outer_segments, outer_error = find_bspline_segments_given_max_segments(kan.outer_spline, max_segments, min_x, max_x)
    max_achieved_error = max(max_achieved_error, outer_error)
    outer_approx = DynamicSplineApproximation(outer_segments, max_error)
    
    return inner_approximations, outer_approx, max_achieved_error
end

kan = KAN(3)
max_segments = 20
inner_approx, outer_approx, max_error = fit_bspline_segments_given_max_segments(kan, max_segments)
x_sample = LinRange(-10, 10, 1000)
p = plot_fitted_segments(kan, inner_approx, outer_approx, x_sample)
display(p)

function verify_kan_with_dp_2(kan::KAN, target_error::Float64, max_total_segments::Int)
    min_x = minimum(kan.inner_splines[1].grid_knots)
    max_x = maximum(kan.inner_splines[1].grid_knots)
    
    # Step 1: Analyze each spline to determine error vs. segment count relationship
    all_splines = vcat(kan.inner_splines, [kan.outer_spline])
    n_splines = length(all_splines)
    
    # Precompute errors for each spline with different segment counts
    error_tables = Vector{Tuple{Vector{Int}, Vector{Float64}}}(undef, n_splines)
    
    # For each spline, compute errors for 1 up to max_segments_per_spline segments
    max_segments_per_spline = 20  # Limit per spline for computation efficiency
    
    for (i, spline) in enumerate(all_splines)
        segment_counts = Int[]
        errors = Float64[]
        
        # Calculate error for different segment counts
        for seg_count in 1:max_segments_per_spline
            segments, error = find_bspline_segments_given_max_segments(spline, seg_count, min_x, max_x)
            push!(segment_counts, seg_count)
            push!(errors, error)
            
            # Early stopping if we've reached low enough error
            if error < target_error / 2  # Conservative threshold
                break
            end
        end
        
        error_tables[i] = (segment_counts, errors)
    end
    
    # Step 2: Use dynamic programming to find optimal segment allocation
    # We'll use a DP table where:
    # - States are (spline_index, segments_used)
    # - Values are minimum achievable max error
    
    # Initialize DP table with infinity
    dp = Dict{Tuple{Int,Int}, Float64}()
    # Track decisions for reconstruction
    decisions = Dict{Tuple{Int,Int}, Int}()
    
    # Base case: allocate 0 to all splines (infinite error)
    for s in 0:max_total_segments
        dp[(0, s)] = Inf
    end
    dp[(0, 0)] = 0.0  # Start with 0 error
    
    # Fill the DP table
    for i in 1:n_splines
        segment_counts, errors = error_tables[i]
        
        for s in 0:max_total_segments
            dp[(i, s)] = Inf  # Default to infinity
            
            # Try different allocations for current spline
            for (j, seg_count) in enumerate(segment_counts)
                if seg_count <= s
                    # Use seg_count segments for spline i
                    prev_error = dp[(i-1, s-seg_count)]
                    new_max_error = max(prev_error, errors[j])
                    
                    if new_max_error < dp[(i, s)]
                        dp[(i, s)] = new_max_error
                        decisions[(i, s)] = seg_count
                    end
                end
            end
        end
    end
    
    # Step 3: Find minimum segments required to achieve target error
    min_segments_required = max_total_segments
    for s in 1:max_total_segments
        if dp[(n_splines, s)] <= target_error
            min_segments_required = s
            break
        end
    end
    
    if dp[(n_splines, min_segments_required)] > target_error
        @warn "Could not achieve target error of $target_error. Minimum error: $(dp[(n_splines, min_segments_required)])"
    end
    
    # Step 4: Reconstruct the optimal allocation
    allocation = zeros(Int, n_splines)
    remaining_segments = min_segments_required
    
    for i in n_splines:-1:1
        allocation[i] = decisions[(i, remaining_segments)]
        remaining_segments -= allocation[i]
    end
    
    # Step 5: Create the approximations with optimal allocations
    inner_approximations = DynamicSplineApproximation[]
    for (i, spline) in enumerate(kan.inner_splines)
        segments, error = find_bspline_segments_given_max_segments(spline, allocation[i], min_x, max_x)
        push!(inner_approximations, DynamicSplineApproximation(segments, error))
    end
    
    outer_segments, outer_error = find_bspline_segments_given_max_segments(
        kan.outer_spline, allocation[end], min_x, max_x)
    outer_approx = DynamicSplineApproximation(outer_segments, outer_error)
    
    achieved_error = dp[(n_splines, min_segments_required)]
    
    # Return details about the optimization
    println("Optimized segment allocation to achieve error ≤ $target_error:")
    println("  Total segments required: $min_segments_required")
    println("  Achieved maximum error: $achieved_error")
    println("  Segment allocation:")
    for i in 1:length(kan.inner_splines)
        println("    Inner spline $i: $(allocation[i]) segments")
    end
    println("    Outer spline: $(allocation[end]) segments")
    
    return inner_approximations, outer_approx, achieved_error, allocation
end

function plot_optimized_kan(kan::KAN, inner_approx::Vector{DynamicSplineApproximation}, 
    outer_approx::DynamicSplineApproximation, achieved_error::Float64, allocation::Vector{Int64}, target_error::Float64)
    # Create visualization
    x_sample = LinRange(-10, 10, 1000)
    p = plot(layout=(5,1), size=(800,1000), legend=:topright)
    
    # Plot inner splines
    for (i, (spline, approx)) in enumerate(zip(kan.inner_splines, inner_approx))
        spline_vals = [eval_spline(spline, x) for x in x_sample]
        plot!(p[i], x_sample, spline_vals, label="Inner Spline $i", linewidth=2)
        
        # Plot segments
        for (j, (start_x, end_x, slope, intercept)) in enumerate(approx.segments)
            x_segment = range(start_x, end_x, length=50)
            y_segment = [slope * x + intercept for x in x_segment]
            plot!(p[i], x_segment, y_segment,
                linestyle=:dash,
                label=(j==1 ? "Approx ($(length(approx.segments)) segments)" : false),
                color=:orange)
            
            # Plot error bounds
            plot!(p[i], x_segment, y_segment .+ approx.max_error,
                linestyle=:dot, color=:red, alpha=0.3,
                label=(j==1 ? "Error Bound: $(round(approx.max_error, digits=3))" : false))
            plot!(p[i], x_segment, y_segment .- approx.max_error,
                linestyle=:dot, color=:red, alpha=0.3, label=false)
        end
        
        title!(p[i], "Inner Spline $i ($(length(approx.segments)) segments)")
        ylabel!(p[i], "Output")
    end
    
    # Plot outer spline
    outer_vals = [eval_spline(kan.outer_spline, x) for x in x_sample]
    plot!(p[4], x_sample, outer_vals, label="Outer Spline", linewidth=2)
    
    for (j, (start_x, end_x, slope, intercept)) in enumerate(outer_approx.segments)
        x_segment = range(start_x, end_x, length=50)
        y_segment = [slope * x + intercept for x in x_segment]
        plot!(p[4], x_segment, y_segment,
            linestyle=:dash,
            label=(j==1 ? "Approx ($(length(outer_approx.segments)) segments)" : false),
            color=:orange)
        
        # Plot error bounds
        plot!(p[4], x_segment, y_segment .+ outer_approx.max_error,
            linestyle=:dot, color=:red, alpha=0.3,
            label=(j==1 ? "Error Bound: $(round(outer_approx.max_error, digits=3))" : false))
        plot!(p[4], x_segment, y_segment .- outer_approx.max_error,
            linestyle=:dot, color=:red, alpha=0.3, label=false)
    end
    
    title!(p[4], "Outer Spline ($(length(outer_approx.segments)) segments)")
    ylabel!(p[4], "Output")
    
    # Plot the whole KAN output
    kan_vals = [kan(x) for x in x_sample]
    plot!(p[5], x_sample, kan_vals,
        label="KAN Output", linewidth=2)
    
    title!(p[5], "KAN Output (Target Error: $target_error, Achieved: $(round(achieved_error, digits=4)))")
    xlabel!(p[5], "x")
    ylabel!(p[5], "Output")
    
    # Add text annotation with allocation summary
    annotation_text = "Segment allocation:\n"
    for i in 1:length(kan.inner_splines)
        annotation_text *= "Inner $i: $(allocation[i])  "
    end
    annotation_text *= "\nOuter: $(allocation[end])  Total: $(sum(allocation))"
    annotate!(p[5], -5, minimum(kan_vals) + 0.1, text(annotation_text, 10, :left, :top))
    
    return p
end

kan = KAN(3)
target_error = 0.01
inner_approx, outer_approx, achieved_error, allocation = verify_kan_with_dp_2(kan, target_error, 90)
p = plot_optimized_kan(kan, inner_approx, outer_approx, achieved_error, allocation, target_error)
display(p)