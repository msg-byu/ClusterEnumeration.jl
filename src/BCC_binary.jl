# Generate clusters for BCC lattice
using ClusterEnumeration 
using Plots
using Spacey
using LinearAlgebra

cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/data/")
default(ms=2,msw=0,legend=:none,st=:scatter)

# Define the BCC lattice
A = [-1 1 1;
     1 -1 1;
     1 1 -1]

# Get the symmetry operations of the lattice     
rots = pointGroup_robust(eachrow(A)...)



# Find all lattice points in specified sphere. These are used to construct clusters
# 5 seconds for up to 20, 50 secs for up to 30, 120 secs for up to 30
@time rpts, norms = genLatticePts(A,rots,35); length(rpts)
sphPts, degen, lengths =genSymEqvPoints(rpts,rots)
lengths
length(sphPts)

# Generate possible pairs from sphere of lattice points
begin
    n = 500
    sphPts, degen, lengths =genSymEqvPoints(rpts[1:n],rots); @show length(sphPts)
    figCands2 = makeFigureCandidates(sphPts,2); @show length(figCands2)
    rdFigs2 = symReduceFigList(figCands2,rots); @show length(rdFigs2)
 #   @time figCands3 = iterAugmentFigures(rdFigs2,sphPts); @show length(figCands3)
#    @time rdFigs3 = symReduceFigList(figCands3,rots); @show length(rdFigs3)
end