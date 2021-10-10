@shorthands bivariate
@shorthands bivariatelegend


"""
test 1
"""
@recipe function plot(layer::T) where {T <: SimpleSDMLayer}
    seriestype --> :heatmap
    if get(plotattributes, :seriestype, :heatmap) in [:heatmap, :contour]
        eltype(layer) <: AbstractFloat || throw(ArgumentError("This plot is only supported for layers with number values ($(eltype(layer)))"))
        aspect_ratio --> 1
        xlims --> extrema(longitudes(layer))
        ylims --> extrema(latitudes(layer))
        lg = copy(layer.grid)
        replace!(lg, nothing => NaN)
        lg = convert(Matrix{Float64}, lg)
        longitudes(layer), latitudes(layer), lg
    elseif get(plotattributes, :seriestype, :histogram) in [:histogram, :density]
        collect(layer)
    end
end

"""
test 2
"""
@recipe function plot(l1::FT, l2::ST; classes::Int=3, p0=colorant"#e8e8e8", p1=colorant"#64acbe", p2=colorant"#c85a5a") where {FT <: SimpleSDMLayer, ST <: SimpleSDMLayer}
    eltype(l1) <: Number || throw(ArgumentError("Plotting is only supported for layers with number values ($(eltype(l1)))"))
    eltype(l2) <: Number || throw(ArgumentError("Plotting is only supported for layers with number values ($(eltype(l2)))"))
    seriestype --> :scatter
    if get(plotattributes, :seriestype, :scatter) in [:scatter, :histogram2d]
        SimpleSDMLayers._layers_are_compatible(l1, l2)
        valid_i = findall(.!isnothing.(l1.grid) .& .!isnothing.(l2.grid))
        l1.grid[valid_i], l2.grid[valid_i]
    elseif get(plotattributes, :seriestype, :bivariate) in [:bivariate]
        SimpleSDMLayers._layers_are_compatible(l1, l2)
        c1 = LinRange(p0, p1, classes)
        c2 = LinRange(p0, p2, classes)
        breakpoints = LinRange(0.0, 1.0, classes+1)
        classified = similar(l1, Int)
        cols = Vector{typeof(p0)}(undef, classes^2)
        for i in 2:length(breakpoints)
            m1 = broadcast(v -> breakpoints[i - 1] <= v <= breakpoints[i], l1)
            for j in 2:length(breakpoints)
                current_class = (i-2)*classes + (j-2) + 1
                @info i j current_class
                m2 = broadcast(v -> breakpoints[j - 1] <= v <= breakpoints[j], l2)
                m = reduce(*, [m1, m2])
                replace!(m, false => nothing)
                classified[keys(m)] = fill(current_class, sum(m))
                cols[current_class] = ColorBlendModes.BlendMultiply(c1[i - 1], c2[j - 1])
            end
        end
        @series begin
            seriescolor := reverse(cols)
            seriestype := :heatmap
            subplot := 1
            legend := false
            convert(Float16, classified)
        end
    elseif get(plotattributes, :seriestype, :bivariatelegend) in [:bivariatelegend]
        c1 = LinRange(p0, p1, classes)
        c2 = LinRange(p0, p2, classes)
        grid --> :none
        ticks --> :none
        legend --> false
        frametype --> :none
        xlims --> (1-0.5, classes+0.5)
        ylims --> (1-0.5, classes+0.5)
        aspect_ratio --> 1
        cols = Vector{typeof(p0)}(undef, classes^2)
        current_class = 1
        m = zeros(Float64, classes, classes)
        for i in 1:classes
            for j in 1:classes
                m[i,j] = current_class
                cols[current_class] = ColorBlendModes.BlendMultiply(c1[i], c2[j])
                current_class += 1
            end
        end
        @series begin
            seriescolor := cols
            seriestype := :heatmap
            m
        end
    end
end
