using LinearAlgebra
using Random
using JuMP

struct BSpline
    grid_knots::Vector{Float64}
    ctrl_points::Vector{Float64}
    degree::Int
end

function BSpline(lb::Float64, ub::Float64, degree::Int; n::Int=20, rng::AbstractRNG=Random.GLOBAL_RNG)
    k = degree
    cpoints::Vector{Float64} = randn(rng, Float64, n + k)
    grid_pts::Vector{Float64} = LinRange(lb, ub, n + 1)
    for i in 1:k
        pushfirst!(grid_pts, lb - (i * (ub - lb) / (2n)))
        push!(grid_pts, ub + (i * (ub - lb) / (2n)))
    end

    BSpline(grid_pts, cpoints, k)
end

function basis(bs::BSpline, x::Float64)
    n = size(bs.grid_knots)[1] - 1
    k = bs.degree
    b = Vector{Vector{Float64}}(undef, k + 1)
    b[1] = (x .>= bs.grid_knots[1:end-1]) .& (x .< bs.grid_knots[2:end])
    for degree in 1:k
        b[degree+1] = zeros(Float64, n - degree)
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
    inner_values = [eval_spline(spline, x) for spline in kan.inner_splines]
    inner_sum = sum(inner_values)
    return eval_spline(kan.outer_spline, inner_sum)
end

function fit_line_through_points(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    return slope, intercept
end

struct DynamicSplineApproximation
    segments::Vector{Tuple{Float64,Float64,Float64,Float64}}
    max_error::Float64
end

function find_bspline_segments_given_max_segments(bs::BSpline, max_segments::Int, min_x::Float64, max_x::Float64)
    n_samples = min(max_segments * 10, 500)
    x_points = collect(LinRange(min_x, max_x, n_samples))
    y_points = [eval_spline(bs, x) for x in x_points]
    n = length(x_points)
    errors = fill(Inf, (n, n))
    for start_idx in 1:n-1
        max_len = min(n-start_idx+1, round(Int, max_x - min_x))
        end_idx = min(start_idx + max_len-1, n)
        end_indices = collect(start_idx:end_idx)
        x_start = x_points[start_idx]
        y_start = y_points[start_idx]
        x_ends = x_points[end_indices]
        y_ends = y_points[end_indices]
        valid_segments = abs.(x_ends .- x_start) .> 1e-10
        slopes = zeros(length(end_indices))
        intercepts = zeros(length(end_indices))
        slopes[valid_segments] = (y_ends[valid_segments] .- y_start) ./ (x_ends[valid_segments] .- x_start)
        intercepts[valid_segments] = y_start .- slopes[valid_segments] .* x_start
        for (end_idx_position, end_idx) in enumerate(end_indices)
            if valid_segments[end_idx_position]
                segment_x = x_points[start_idx:end_idx]
                segment_y = y_points[start_idx:end_idx]
                predicted_y = slopes[end_idx_position] .* segment_x .+ intercepts[end_idx_position]
                errors[start_idx,end_idx] = maximum(abs.(predicted_y .- segment_y))
            end
        end
    end
    dp = fill(Inf, (max_segments, n_samples))
    back = zeros(Int, (max_segments, n_samples))
    for j in 1:n_samples
        dp[1,j] = errors[1,j]
    end
    for i in 2:max_segments
        for j in i:n_samples
                prev_errors = dp[i-1,1:j-1]
                curr_errors = [errors[k,j] for k in 1:j-1]
                total_errors = max.(prev_errors, curr_errors)
                min_error, min_idx = findmin(total_errors)
                dp[i,j] = min_error
                back[i,j] = min_idx
        end
    end
    segments = Tuple{Float64, Float64, Float64, Float64}[]
    curr_seg = n_samples
    for i in max_segments:-1:1
        prev_seg = i > 1 ? back[i,curr_seg] : 1
        x1, y1 = x_points[prev_seg], y_points[prev_seg]
        x2, y2 = x_points[curr_seg], y_points[curr_seg]
        slope, intercept = fit_line_through_points(x1, y1, x2, y2)
        pushfirst!(segments, (x1, x2, slope, intercept))
        curr_seg = prev_seg
    end
    return segments, dp[max_segments,n_samples]
end

function calculate_lipschitz_constant(bs::BSpline, min_x::Float64, max_x::Float64)
    # For a more accurate Lipschitz constant, we sample the derivative
    n_samples = 1000
    x_points = collect(LinRange(min_x, max_x, n_samples))
    
    # Approximate the derivative using finite differences
    max_derivative = 0.0
    delta = (max_x - min_x) / (n_samples * 10)
    
    for x in x_points
        y1 = eval_spline(bs, x)
        y2 = eval_spline(bs, x + delta)
        derivative = abs((y2 - y1) / delta)
        max_derivative = max(max_derivative, derivative)
    end
    
    return max_derivative
end


function fit_kan_given_total_error_lipschitz(kan::KAN, target_error::Float64, max_total_segments::Int)
    min_x = minimum(kan.inner_splines[1].grid_knots)
    max_x = maximum(kan.inner_splines[1].grid_knots)
    
    # Step 1: Calculate Lipschitz constant for the outer spline
    outer_lipschitz = calculate_lipschitz_constant(kan.outer_spline, min_x, max_x)
    
    # Step 2: Separate inner splines and outer spline
    inner_splines = kan.inner_splines
    outer_spline = kan.outer_spline
    n_inner_splines = length(inner_splines)
    
    # Step 3: Compute error tables for all splines
    inner_error_tables = Vector{Tuple{Vector{Int}, Vector{Float64}}}(undef, n_inner_splines)
    
    # Calculate error tables for inner splines
    for (i, spline) in enumerate(inner_splines)
        segment_counts = Int[]
        errors = Float64[]
        for seg_count in 1:max_total_segments
            segments, error = find_bspline_segments_given_max_segments(spline, seg_count, min_x, max_x)
            push!(segment_counts, seg_count)
            push!(errors, error)
            if error < target_error / (2 * outer_lipschitz * n_inner_splines)
                # We can stop if the error is small enough that even with propagation
                # it would be less than half the target error
                break
            end
        end
        inner_error_tables[i] = (segment_counts, errors)
    end
    
    # Calculate error table for outer spline
    outer_segment_counts = Int[]
    outer_errors = Float64[]
    for seg_count in 1:max_total_segments
        segments, error = find_bspline_segments_given_max_segments(outer_spline, seg_count, min_x, max_x)
        push!(outer_segment_counts, seg_count)
        push!(outer_errors, error)
        if error < target_error / 2
            # We can stop if the error is less than half the target error
            break
        end
    end
    
    # Step 4: Modified dynamic programming approach to account for error propagation
    
    # The key difference from the original function:
    # Instead of tracking max error across all splines, we'll track:
    # 1. The sum of inner spline errors (which will be multiplied by Lipschitz constant)
    # 2. The outer spline error separately
    
    # Define DP state: dp[i][s] = minimum sum of errors for first i inner splines using s segments
    dp = fill(Inf, (n_inner_splines + 1, max_total_segments + 1))
    decisions = fill(0, (n_inner_splines + 1, max_total_segments + 1))
    
    # Base case: no splines, no error
    dp[1, 1] = 0.0
    
    # Fill DP table for inner splines only
    for i in 1:n_inner_splines
        segment_counts, errors = inner_error_tables[i]
        for s in 0:max_total_segments
            dp[i+1, s+1] = Inf
            for (j, seg_count) in enumerate(segment_counts)
                if seg_count <= s
                    # Previous sum of errors
                    prev_sum_error = dp[i, s-seg_count+1]
                    
                    # Add current spline's error to the sum
                    new_sum_error = prev_sum_error + errors[j]
                    
                    if new_sum_error < dp[i+1, s+1]
                        dp[i+1, s+1] = new_sum_error
                        decisions[i+1, s+1] = seg_count
                    end
                end
            end
        end
    end
    
    # Step 5: Find optimal allocation considering outer spline and error propagation
    min_segments_required = max_total_segments
    min_total_error = Inf
    best_outer_segments = 1
    
    # Try different allocations for outer spline
    for outer_seg_idx in 1:length(outer_segment_counts)
        outer_seg = outer_segment_counts[outer_seg_idx]
        outer_err = outer_errors[outer_seg_idx]
        
        remaining_for_inner = max_total_segments - outer_seg
        
        if remaining_for_inner >= n_inner_splines  # Need at least one segment per inner spline
            # For each possible inner segment allocation
            for inner_seg in n_inner_splines:remaining_for_inner
                inner_sum_error = dp[n_inner_splines+1, inner_seg+1]
                
                if inner_sum_error < Inf
                    # Calculate total error with propagation
                    propagated_error = outer_lipschitz * inner_sum_error
                    total_error = outer_err + propagated_error
                    
                    # Check if this allocation is better
                    if total_error <= target_error && (inner_seg + outer_seg) < min_segments_required
                        min_segments_required = inner_seg + outer_seg
                        min_total_error = total_error
                        best_outer_segments = outer_seg
                    end
                end
            end
        end
    end
    
    # Step 6: Reconstruct the allocation from DP decisions
    if min_total_error > target_error
        @warn "Could not achieve target error of $target_error. Minimum error: $min_total_error"
    end
    
    inner_allocation = zeros(Int, n_inner_splines)
    remaining_segments = min_segments_required - best_outer_segments
    
    # Backtrack to find inner allocations
    for i in n_inner_splines:-1:1
        inner_allocation[i] = decisions[i+1, remaining_segments+1]
        remaining_segments -= inner_allocation[i]
    end
    
    # Step 7: Build the approximations with the optimal allocation
    inner_approximations = DynamicSplineApproximation[]
    inner_errors = Float64[]
    
    for (i, spline) in enumerate(inner_splines)
        segments, error = find_bspline_segments_given_max_segments(
            spline, inner_allocation[i], min_x, max_x)
        push!(inner_approximations, DynamicSplineApproximation(segments, error))
        push!(inner_errors, error)
    end
    
    outer_segments, outer_error = find_bspline_segments_given_max_segments(
        outer_spline, best_outer_segments, min_x, max_x)
    outer_approx = DynamicSplineApproximation(outer_segments, outer_error)
    
    # Calculate final achieved error with propagation
    propagated_inner_error = outer_lipschitz * sum(inner_errors)
    achieved_error = outer_error + propagated_inner_error
    
    # For debugging and analysis
    println("Lipschitz constant of outer spline: ", outer_lipschitz)
    println("Inner spline allocations: ", inner_allocation)
    println("Inner spline errors: ", inner_errors)
    println("Sum of inner errors: ", sum(inner_errors))
    println("Propagated inner error: ", propagated_inner_error)
    println("Outer spline allocation: ", best_outer_segments)
    println("Outer spline error: ", outer_error)
    println("Total achieved error: ", achieved_error)
    println("Total segments used: ", sum(inner_allocation) + best_outer_segments)
    
    # Output the dynamic programming table
    println("\n--- Dynamic Programming Table ---")
    println("dp[i,s] = minimum sum of errors for first i inner splines using s segments")
    
    # Determine a reasonable limit for display to avoid overwhelming output
    max_rows_to_show = min(n_inner_splines + 1, 10)
    max_cols_to_show = min(min_segments_required + 10, 30)
    
    # Print column headers (segments)
    print("i\\s |")
    for s in 0:max_cols_to_show-1
        @printf("%8d |", s)
    end
    println()
    println("-" ^ (9 * max_cols_to_show + 5))
    
    # Print each row (splines)
    for i in 1:max_rows_to_show
        @printf("%3d |", i-1)
        for s in 0:max_cols_to_show-1
            if dp[i, s+1] == Inf
                print("     Inf |")
            else
                @printf("%8.4f |", dp[i, s+1])
            end
        end
        println()
    end
    
    # Print decisions table if needed
    println("\n--- Decisions Table ---")
    println("decisions[i,s] = number of segments allocated to spline i when using s total segments")
    
    # Print column headers (segments)
    print("i\\s |")
    for s in 0:max_cols_to_show-1
        @printf("%8d |", s)
    end
    println()
    println("-" ^ (9 * max_cols_to_show + 5))
    
    # Print each row (splines)
    for i in 1:max_rows_to_show
        @printf("%3d |", i-1)
        for s in 0:max_cols_to_show-1
            @printf("%8d |", decisions[i, s+1])
        end
        println()
    end
    
    return inner_approximations, outer_approx, achieved_error, vcat(inner_allocation, best_outer_segments)
end

function plot_optimized_kan(kan::KAN, inner_approx::Vector{DynamicSplineApproximation}, 
    outer_approx::DynamicSplineApproximation, achieved_error::Float64, allocation::Vector{Int64}, target_error::Float64)
    x_sample = LinRange(-10, 10, 1000)
    p = plot(layout=(5,1), size=(800,1000), legend=:topright)
    for (i, (spline, approx)) in enumerate(zip(kan.inner_splines, inner_approx))
        spline_vals = [eval_spline(spline, x) for x in x_sample]
        plot!(p[i], x_sample, spline_vals, label="Inner Spline $i", linewidth=2)
        for (j, (start_x, end_x, slope, intercept)) in enumerate(approx.segments)
            x_segment = range(start_x, end_x, length=50)
            y_segment = [slope * x + intercept for x in x_segment]
            plot!(p[i], x_segment, y_segment,
                linestyle=:dash,
                label=(j==1 ? "Approx ($(length(approx.segments)) segments)" : false),
                color=:orange)
            plot!(p[i], x_segment, y_segment .+ approx.max_error,
                linestyle=:dot, color=:red, alpha=0.3,
                label=(j==1 ? "Error Bound: $(round(approx.max_error, digits=3))" : false))
            plot!(p[i], x_segment, y_segment .- approx.max_error,
                linestyle=:dot, color=:red, alpha=0.3, label=false)
        end
        title!(p[i], "Inner Spline $i ($(length(approx.segments)) segments)")
        ylabel!(p[i], "Output")
    end
    outer_vals = [eval_spline(kan.outer_spline, x) for x in x_sample]
    plot!(p[4], x_sample, outer_vals, label="Outer Spline", linewidth=2)
    
    for (j, (start_x, end_x, slope, intercept)) in enumerate(outer_approx.segments)
        x_segment = range(start_x, end_x, length=50)
        y_segment = [slope * x + intercept for x in x_segment]
        plot!(p[4], x_segment, y_segment,
            linestyle=:dash,
            label=(j==1 ? "Approx ($(length(outer_approx.segments)) segments)" : false),
            color=:orange)
        plot!(p[4], x_segment, y_segment .+ outer_approx.max_error,
            linestyle=:dot, color=:red, alpha=0.3,
            label=(j==1 ? "Error Bound: $(round(outer_approx.max_error, digits=3))" : false))
        plot!(p[4], x_segment, y_segment .- outer_approx.max_error,
            linestyle=:dot, color=:red, alpha=0.3, label=false)
    end
    title!(p[4], "Outer Spline ($(length(outer_approx.segments)) segments)")
    ylabel!(p[4], "Output")
    kan_vals = [kan(x) for x in x_sample]
    plot!(p[5], x_sample, kan_vals,
        label="KAN Output", linewidth=2)
    title!(p[5], "KAN Output (Target Error: $target_error, Achieved: $(round(achieved_error, digits=4)))")
    xlabel!(p[5], "x")
    ylabel!(p[5], "Output")
    return p
end

kan = KAN(3, seed=42)
target_error = 0.5
inner_approx, outer_approx, achieved_error, allocation = fit_kan_given_total_error_lipschitz(kan, target_error, 300)
p = plot_optimized_kan(kan, inner_approx, outer_approx, achieved_error, allocation, target_error)
display(p)