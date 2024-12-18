function subalgred(strings::AbstractMatrix{Int8})::Matrix{Int8}
    abstrings = strings[:, 1]
    multstrings = strings[:, 1:1]
    i = 1
    while !all(string in eachcol(multstrings) for string in eachcol(strings))
        while true
            string = @view strings[:, i]
            i += 1
            if !(string in eachcol(multstrings))
                append!(abstrings, string)
                break
            end
        end
        for multstring in eachcol(multstrings)
            multstrings = [multstrings pauliprod(view(strings, :, i-1), multstring)[1]]
        end
        multstrings = [multstrings view(strings, :, i-1)]
    end
    return reshape(abstrings, size(multstrings, 1), :)
end

function symsubspaces(kstrings::AbstractMatrix{Int8}, abstrings::AbstractMatrix{Int8})::Vector{Matrix{Int8}}
    strings = kstrings
    symstrings = Vector{Matrix{Int8}}(undef, size(abstrings, 2))
    for i in axes(abstrings, 2)
        temp = Vector{Int8}(undef, 0)
        j = 1
        while j <= size(strings, 2)
            if !(pauliprod(view(abstrings, :, i), view(strings, :, j))[3])
                append!(temp, view(strings, :, j))
                strings = strings[:, 1:end .!= j]
            else
                j += 1
            end
        end
        symstrings[i] = reshape(temp, size(abstrings, 1), :)
    end
    return symstrings
end

function cleangenerators!(symgenerators::AbstractVector{Matrix{Int8}}, abstrings::AbstractMatrix{Int8})::Matrix{Int8}
    
    temp = Matrix{Int8}(undef, size(abstrings, 1), 0)
    while temp in symgenerators
        ind = findfirst(==(temp), symgenerators)
        popat!(symgenerators, ind)
        abstrings = abstrings[:, 1:end .!= ind]
    end
    return abstrings
end