function validate_grid_size_and_length(sz, len, x, y, z)
    if isnothing(len) && (isnothing(x) || isnothing(y) || isnothing(z))
        throw(ArgumentError("Must supply length or x, y, z keyword arguments."))
    end

    if !isnothing(len) && !isnothing(x) && !isnothing(y) && !isnothing(z)
        throw(ArgumentError("Cannot specify both length and x, y, z keyword arguments."))
    end

    length(sz) == 3        || throw(ArgumentError("length($sz) must be 3."))
    all(isa.(sz, Integer)) || throw(ArgumentError("size=$sz should contain integers."))
    all(sz .>= 1)          || throw(ArgumentError("size=$sz must be nonzero and positive!"))

    if !isnothing(len)
        length(len) == 3       || throw(ArgumentError("length($len) must be 3."))
        all(isa.(len, Number)) || throw(ArgumentError("length=$len should contain numbers."))
        all(len .>= 0)         || throw(ArgumentError("length=$len must be nonzero and positive!"))
    end

    if isnothing(len)
        function coord2xyz(c)
            c == 1 && return "x"
            c == 2 && return "y"
            c == 3 && return "z"
        end

        for (i, c) in enumerate((x, y, z))
            name = coord2xyz(i)
            length(c) == 2       || throw(ArgumentError("$name length($c) must be 2."))
            all(isa.(c, Number)) || throw(ArgumentError("$name=$c should contain numbers."))
            c[2] >= c[1]         || throw(ArgumentError("$name=$c should be an increasing interval."))
        end
    end
end
