using Cartan
using Test

@testset "Cartan.jl" begin
    tfimham = Int8[[2, 2, 1] [1, 2, 2] [4, 1, 1] [1, 4, 1] [1, 1, 4]]
    tfimdla = Int8[[2, 2, 1] [1, 2, 2] [4, 1, 1] [1, 4, 1] [1, 1, 4] [1, 2, 3] [1, 3, 2] [2, 3, 1] [3, 2, 1] [3, 4, 2] [3, 3, 1] [2, 4, 3] [2, 4, 2] [1, 3, 3] [3, 4, 3]]
    tfimk = Int8[[1, 2, 3] [1, 3, 2] [2, 3, 1] [3, 2, 1] [3, 4, 2] [2, 4, 3]]
    tfimm = Int8[[2, 2, 1] [1, 2, 2] [4, 1, 1] [1, 4, 1] [1, 1, 4] [3, 3, 1] [2, 4, 2] [1, 3, 3] [3, 4, 3]]
    tfimh = Int8[[4, 1, 1] [1, 4, 1] [1, 1, 4]]

    cd = involutionlessdecomp(tfimham)
    k = cd["k"]
    m = cd["m"]
    h = cd["h"]
    cartandecomp_test = false
    if size(k, 2) != size(tfimk, 2) || size(m, 2) != size(tfimm, 2) || size(h, 2) != size(tfimh, 2)
        @test cartandecomp_test
    else
        for string in eachcol(k)
            string in eachcol(tfimk) && continue
            @test cartandecomp_test
            break
        end
        for string in eachcol(m)
            string in eachcol(tfimm) && continue
            @test cartandecomp_test
            break
        end
        for string in eachcol(h)
            string in eachcol(tfimh) && continue
            @test cartandecomp_test
            break
        end
        @test !cartandecomp_test
    end
end