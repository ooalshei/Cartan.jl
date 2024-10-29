function subalgred(strings::Matrix{Int8})::Matrix{Int8}
    abstrings = strings[:, 1:1]
    multstrings = abstrings
    i = 1
    while !all(string in eachcol(multstrings) for string in eachcol(strings))
        while true
            string = strings[:, i]
            if !(string in eachcol(multstrings))
                abstrings = [abstrings string]
                i += 1
                break
            else
                i += 1
            end
        end
        for j in axes(multstrings, 2)
            multstrings = [multstrings pauliprod(strings[:, i-1], multstrings[:, j])[1]]
        end
        multstrings = [multstrings strings[:, i-1]]
    end
    return abstrings
end

function symsubspaces(kstrings::Matrix{Int8}, abstrings::Matrix{Int8})::Vector{Matrix{Int8}}
    strings = copy(kstrings)
    symstrings = Vector{Matrix{Int8}}(undef, size(abstrings, 2))
    for i in axes(abstrings, 2)
        temp = Matrix{Int8}(undef, size(strings, 1), 0)
        j = 1
        while j <= size(strings, 2)
            if !(pauliprod(abstrings[:, i], strings[:, j])[3])
                temp = [temp strings[:, j]]
                strings = strings[:, 1:end .!= j]
            else
                j += 1
            end
        end
        symstrings[i] = temp
    end
    return symstrings
end

function cleangenerators!(symgenerators::Vector{Matrix{Int8}}, abstrings::Matrix{Int8})::Matrix{Int8}
    
    temp = Matrix{Int8}(undef, size(abstrings, 1), 0)
    while temp in symgenerators
        ind = findfirst(==(temp), symgenerators)
        popat!(symgenerators, ind)
        abstrings = abstrings[:, 1:end .!= ind]
    end
    return abstrings
end
