using Cartan
using Test

@testset "Cartan.jl" begin
    h = Int8[[1, 1, 2, 2] [1, 1, 3, 3] [1, 1, 4, 4] [2, 2, 1, 1] [3, 3, 1, 1] [4, 4, 1, 1] [2, 2, 3, 3] [2, 2, 4, 4] [3, 3, 2, 2] [3, 3, 4, 4] [4, 4, 2, 2] [4, 4, 3, 3]]
    b = Int8[[1, 1, 2, 2] [1, 1, 3, 3] [2, 2, 1, 1] [3, 3, 1, 1]]
    k = Int8[[1, 3, 2, 4] [1, 2, 3, 4] [1, 4, 2, 3] [1, 2, 4, 3] [1, 4, 3, 2] [1, 3, 4, 2] [3, 2, 4, 1] [2, 3, 4, 1] [4, 2, 3, 1] [2, 4, 3, 1] [4, 3, 2, 1] [3, 4, 2, 1] [3, 1, 4, 2] [3, 2, 1, 4] [4, 1, 3, 2] [4, 2, 1, 3] [2, 1, 4, 3] [2, 3, 1, 4] [4, 3, 1, 2] [4, 1, 2, 3] [2, 4, 1, 3] [2, 1, 3, 4] [3, 4, 1, 2] [3, 1, 2, 4]]
    symk = [Int8[[1, 3, 2, 4] [1, 4, 2, 3] [1, 4, 3, 2] [1, 3, 4, 2] [3, 2, 4, 1] [2, 3, 4, 1] [4, 2, 3, 1] [2, 4, 3, 1] [3, 1, 4, 2] [3, 2, 1, 4] [4, 1, 3, 2] [4, 2, 1, 3] [2, 3, 1, 4] [4, 1, 2, 3] [2, 4, 1, 3] [3, 1, 2, 4]],
    Int8[[1, 2, 3, 4] [1, 2, 4, 3] [4, 3, 2, 1] [3, 4, 2, 1] [2, 1, 4, 3] [4, 3, 1, 2] [2, 1, 3, 4] [3, 4, 1, 2]],
    Matrix{Int8}(undef, 4, 0), Matrix{Int8}(undef, 4, 0)]

    ab = subalgred(h)
    subalgred_test = false
    if size(ab, 2) != size(b, 2)
        @test subalgred_test
    else
        for string in eachcol(ab)
            string in eachcol(b) && continue
            @test subalgred_test
            break
        end
        @test !subalgred_test
    end

    kk = symsubspaces(k, ab)
    symsubspaces_test = false
    for i in eachindex(symk)
        if size(kk[i], 2) != size(symk[i], 2)
            @test symsubspaces_test
            break
        else
            for string in eachcol(kk[i])
                string in eachcol(symk[i]) && continue
                @test symsubspaces_test
                break
            end
            i == length(symk) ? (@test !symsubspaces_test) : nothing
        end
    end
end