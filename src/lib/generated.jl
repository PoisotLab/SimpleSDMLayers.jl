ops = Dict(
    :Base => [:sum, :maximum, :minimum],
    :Statistics => [:mean, :median, :std]
    )

for op in ops
    mod = op.first
    for fun in op.second
        @info fun
        eval(quote
            using $mod
            import $mod: $fun
        end)
        for ty in (:SimpleSDMResponse, :SimpleSDMPredictor)
            @info ty
            eval(quote
            """
                $($mod).$($fun)(l::$($ty){T}) where {T <: Number}

            Applies `$($fun)` (from `$($mod)`) to an object of type `$($ty)`. This function has been
            automatically generated. Note that this function is only applied to the
            non-`nothing` elements of the layer, and has no method to work on the `dims`
            keyword; the grid itself can be extracted with `convert(Matrix, l)`.
            """
            function $mod.$fun(l::$ty{T}) where {T <: Number}
                return $mod.$fun(filter(!isnothing, l.grid))
            end
            end)
        end
    end
end
