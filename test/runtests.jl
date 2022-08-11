using ClusterEnumeration
using Test
using Spacey

@testset "ClusterEnumeration.jl" begin
    @test typeof(read_lattice_vectors()) == Matrix{Float64}
    A = read_lattice_vectors()
    #@test gen_points_in_supercell(A,3)[1]|>length == 10
    #ptsTest,_=gen_points_in_supercell(A,40)
    _, rots = pointGroup_robust(eachcol(A)...)
    #_,dgn,_=genSymEqvPoints(ptsTest,rots)
    #@test testDegeneracies(dgn)
end
