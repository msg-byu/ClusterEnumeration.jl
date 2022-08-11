using ClusterEnumeration
using Revise
using Spacey
names(Main)
names(ClusterEnumeration)

cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/develop")

A = read_lattice_vectors()
_,rots = pointGroup_robust(eachcol(A)...)
#plot!(1:1000,cbrt.(1:1000)*cbrt(16π/3)/2)
pts,_ = gen_points_in_supercell(A,2)
ptsTest,_=gen_points_in_supercell(A,30)
fgCnds, fgDgn, ptDists = genSymEqvPoints(pts,rots)