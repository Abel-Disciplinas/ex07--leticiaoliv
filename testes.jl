using Base.Test, BenchmarkTools

include("lu.jl")

function tests()
    n = 100
    A = rand(n, n) + n * I
    b = A * ones(n)

    LU = declu(copy(A))
    for rlu = [reslu, reslu2, reslu3, reslu4]
        x = rlu(LU, copy(b))
        @test norm(x - ones(n)) < 1e-12
        t = @btime $rlu($LU, $(copy(b))) samples = 10 evals = 1
    end
end

tests()
