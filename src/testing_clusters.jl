using Plots
using LinearAlgebra
using DelimitedFiles
using Spacey

cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl")
include("find_clusters.jl")
using .clusters

cd("pairs")

A = clusters.read_lattice_vectors()
_, rots = pointGroup_robust(eachcol(A)...) 
Rpts = gen_points_in_sphere(A,100) # At 5000, takes a long time...10-15 mins

# make list of equivalent points from small set of Rpts
pairFigs = Rpts
scalefontsizes(1.3)
plot(norm.(pairFigs),xlabel="Pair figure number",ylabel="Norm of figure (a0)",legend=:none,lw=2) # Probably start to miss things around 10-20% of the total...
histogram(norm.(pairFigs),legend=:none,xlabel="Norm of figure (a0)",ylabel="Number of figures",color=:cyan,title="Pair Figures vs. Distance")
r = 0:0.5:40
num = [count(norm.(pairFigs).< i) for i in r]
plot(r,num,legend=:none,xlabel="Pair length (a0)",ylabel="Number of figures inside length",color=:red,lw=3)
# Should climb as the cube of the length. Breaks down around length = 10.0
nPairs = 1000
norm.(Rpts)[nPairs] 

pairFigs = [hcat([0.,0.,0.],i) for i ∈ Rpts[1:nPairs]]
nV = ones(Int,nPairs)*2
svecs = [[1,1] for i ∈ 1:nPairs]
dists = norm.(pairFigs[1:nPairs])
id = collect(1:nPairs)

# Add point cluster at the beginning of the list
clust = vcat([[0;0;0]],pairFigs)
svecs = pushfirst!(svecs,[1]) # should this be just [1]?
prepend!(nV,1)
prepend!(dists,0.0)
prepend!(id,0)
write_clusters(clust, svecs, nV, dists, id,"clusters.2b_"*string(nPairs))
clust, svecs, nV, dists, id = read_clusters_from_file("clusters.2b_"*string(nPairs))
#>  cp clusters.2b_496 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.2b_3500
m = readdlm("pimat.2b_1000")
rank(m) # still 137 when we go out to 3500 pairs (last l.i. pair is at 368)
# The last(?) pair occurs around ≈9.3.
s = svd(m);
plot(s.S[1:141],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="Size of singular values",xlabel="Singular values")
# I wonder if we should keep a few of the columns that have a low singular value. I'll try keeping the next three after 137
_,R=qr(m);
plot(abs.(diag(R)),st=:scatter,msw=0,color=:red,ms=2,yaxis=:log,legend=:none,title="R matrix diag entry sizes")
plot(sort(abs.(diag(R)))[1:141],st=:scatter,msw=0,color=:red,ms=2,yaxis=:log,legend=:none,title="R matrix diag entries sorted by size, truncated")
on = get_nonzero_index(R,2.5e-14)
count(on)
findlast(on)
# The longest one still is 137 in the sing val list
#on = vcat(on,falses(nPairs-size(m,1)+2))
idx = collect(1:size(m,2))
c =  [i ? :blue : :orange for i ∈ on[3:end]]
sz = [i ? 4 : 2 for i ∈ on[3:end]]
plot(idx[3:end],dists[2:end],label="Independent pairs",legend=(0.7,0.7),color=c,imagesize=(1600,1000),dpi=300,ms=sz,xlabel="Pair number",ylabel="Pair length (a0)",st=:scatter,msw=0,markershape=:diamond)
# 137 l.i. pairs (but kept three extra, just in case). Last one is #365. Length is ≈9.25
dists[365]
# Eliminate the dependent terms and write a new cluster file
popfirst!(on) # Drop empty cluster flag from mask
write_clusters(clust[on], svecs[on], nV[on], dists[on], id[on],"clusters.2b_li139")
#>  cp clusters.2b_li_137 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2b_li139")
rank(m) # Still have a 137 rank matrix with no extra columns included.
s = svd(m);
plot(s.S,msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="Size of singular values",xlabel="Singular values",title="After throwing away dependent pairs")

# Now, let's get some 3-body clusters
pts =  make_eqvPoints(Rpts[1:145],rots)
figs3bt = make_figure_candidates(pts,3)
@time fullfigs3b = reduce_figList(figs3bt[1:1_200_000],rots) # ~10mins for 200k
figs3b = fullfigs3b[1:end]
n3b = length(figs3b)
#figs3b, sv3, nV3, dists3, id3 = read_clusters_from_file("clusters.2b3b_7921")
#n3b = length(figs3b[nV3.==3])
clust3 = vcat(clust[on],figs3b)
sv3 = vcat(svecs[on],[[1,1,1] for i ∈ 1:n3b])
id3 =  vcat(id[on],1:n3b)
dists3 = vcat(dists[on],MADfromCOM.(figs3b))
nV3 = vcat(nV[on],ones(Int,n3b)*3)
write_clusters(clust3, sv3, nV3, dists3, id3,"clusters.2b3b_"*string(n3b))
#>  cp clusters.2b3b_4470 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.2b_li136
m = readdlm("pimat.2b3b_"*string(n3b))
rank(m) # Finally got to rank=415. Not entirely sure this is all of the clusters, but it's as many as I had on the last round. # wow! I got 418 with the new way of getting the cluster candidates (that solves the problem of missing 3 at 6-bodies case.)
# Now back to 415 with a longer list of 3-bodies. What is going on?!
# Ok!  I took the first 4470 three bodies by working in two parts. Now the rank is 429!
# On re-running, the rank fell to 424. Grr! Regenerated 3-bodies with slightly larger pool. Still 424
# Added sorting to the make_figure_candidates routine changed things a bit. Fewer unique clusters. Now must go farther out to find l.i. clusters
# Kept 139 pairs (extra 3). Rank 426 with ~7900 triplets
# with ~11,800 triplets, rank is 430


plot(dists3[140:end],st=:scatter,msw=0,ms=2,xlabel="Triplet number",ylabel="Triplet diameter",legend=:none)
rank(m[:,1:2346]) # rank 430 => 412 implies 18 of the l.i. 3-bodies are right of the diagonal
# _,R3 = qr(m);
# on3 = get_nonzero_index(R3,8e-14)
# count(on3) # This keeps a few we might not need but we'll see at n=6...
# plot(370:425,sort(abs.(diag(R3)),rev=true)[370:425],st=:scatter,msw=0,color=:red,ms=2,yaxis=:log,legend=:none,title="R3 matrix diag entries sorted by size, truncated")
# s = svd(m);
# plot(400:435,s.S[400:435],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="(2+3) Size of singular values",xlabel="Singular values")
#
t = get_leftmost_indep_columns(m)
rank(m[:,t]) # Still have 430 rank
mt = m[:,t]
s = svd(mt);
plot(400:455,s.S[400:455],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="(2+3) Size of singular values",xlabel="Singular values")
on = [i ∈ t for i ∈ 1:length(dists3)]
c =  [i ? :blue : :orange for i ∈ on[3:end]]
sz = [i ? 4 : 2 for i ∈ on[3:end]]
plot(dists3[140:end],xlabel="cluster Number", title="3-body Clusters",ylabel="Diameter",color=c,ms=sz,label="Independent 3-bodies",legend=:bottom)


popfirst!(t) # Don't need the first element -- it refers to the empty cluster, not listed in clusters.out
t .-= 1 # Subtract 1 from each entry because clusters array and pi matrix column index are shifted by 1
clust3li = clust3[t]
sv3li = sv3[t]
id3li =  id3[t]
dists3li = dists3[t]
nV3li = nV3[t]
write_clusters(clust3li, sv3li, nV3li, dists3li, id3li,"clusters.2b3b_li460")
#>  cp clusters.2b3b_li418 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2b3b_li460")
rank(m)

on3 = falses(length(clust3))
on3[t] .= 1 # set mask for l.i. clusters
idx3 = collect(1:length(clust3))
c =  [i ? :blue : :orange for i ∈ on3]
sz = [i ? 4 : 2 for i ∈ on3]
plot(idx3[140:end],dists3[140:end],label="Independent triplets",legend=(0.7,0.65),color=c,imagesize=(1600,1000),ms=sz,xlabel="Triplet number",ylabel="Triplet \"diameter\"",st=:scatter,msw=0,markershape=:diamond)
# It's hard to see how many lin.dep. triples there are at small diameters. Rethink this viz

# 4 body terms
pts = make_eqvPoints(Rpts[1:13],rots) 
figs4bt = make_figure_candidates(pts,4) # if 248 pts then 2.5 M candidates
@time figs4b = reduce_figList(figs4bt[1:800_000],rots) # Using the first 800K, took >1 hour to generate, generated  8202
# Read in from file if we don't want to regenerate
# figs4b, sv4, nV4, dists4, id4 = read_clusters_from_file("clusters.2-4b_4307")
# mask = nV4.==4
# n4b = count(mask)
n4b = length(figs4b)
clust4 = vcat(clust3li,figs4b)
sv4 = vcat(sv3li,[[1,1,1,1] for i ∈ 1:n4b])
id4 =  vcat(id3li,1:n4b)
dists4 = vcat(dists3li,MADfromCOM.(figs4b))
nV4 = vcat(nV3li,ones(Int,n4b)*4)

write_clusters(clust4, sv4, nV4, dists4, id4,"clusters.2-4b_"*string(n4b))
#>  cp clusters.2-4b_2796 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-4b_"*string(n4b))
rank(m)
rank(m[:,1:2346]) # 9 l.i. columns to the right of the diagonal
t = get_leftmost_indep_columns(m,)
rank(m[:,t]) # Still have 940 rank
mt = m[:,t]
s = svd(mt);
plot(925:945,s.S[925:945],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="(2+3) Size of singular values",xlabel="Singular values")
_,Rtest = qr(mt);
plot(abs.(diag(Rtest)),yaxis=:log,legend=:none)
plot(sort(abs.(diag(Rtest)),rev=true),yaxis=:log)
plot([960:1000],sort(abs.(diag(Rtest))[960:1000],rev=true),yaxis=:log,st=:scatter,ms=2,msw=0)

# The singular values and rank don't closely match the size of entries in the R'S
# So let's keep some extras from the R's, but not all:
on = get_nonzero_index(Rtest,2.2e-13); count(on)
s = copy(t) # Keep s, just in case.
t = s[on]

popfirst!(t) # Don't need the first element -- it refers to the empty cluster, not listed in clusters.out
t .-= 1 # Subtract 1 from each entry because clusters array and pi matrix column index are shifted by 1
clust4li = clust4[t]
sv4li = sv4[t]
id4li =  id4[t]
dists4li = dists4[t]
nV4li = nV4[t]
write_clusters(clust4li, sv4li, nV4li, dists4li, id4li,"clusters.2-4b_li"*string(length(t)))
#>  cp clusters.2b3b_li418 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-4b_li"*string(length(t)))
rank(m) # Monday: 934, 940

# 5 body terms
pts = make_eqvPoints(Rpts[1:7],rots) 
@time figs5bt = make_figure_candidates(pts,5) # up to 6 Rpts was just a few seconds, up to 7took 90 seconds and generated 12M candidates. about 50 mins to reduce the first 800K 
@time figs5bfull = reduce_figList(figs5bt[1:800_000],rots)
n5b = length(figs5b)
figs5b = figs5bfull[1:n5b]
#figs5b, sv5, nV5, dists5, id5 = read_clusters_from_file("clusters.2-5b_4500")
#mask=nV5.==5
#n5b = count(mask)
clust5 = vcat(clust4li,figs5b)
sv5 = vcat(sv4li,[[1,1,1,1,1] for i ∈ 1:n5b])
id5 =  vcat(id4li,1:n5b)
dists5 = vcat(dists4li,MADfromCOM.(figs5b))
nV5 = vcat(nV4li,ones(Int,n5b)*5)
write_clusters(clust5, sv5, nV5, dists5, id5,"clusters.2-5b_"*string(n5b))
#>  cp clusters.2-5b_2477 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-5b_"*string(n5b))
rank(m)
rank(m[:,1:2346]) # 37 l.i. columns to the right of the diagonal
t = get_leftmost_indep_columns(m)

r=rank(m[:,t]) # Still have 1364 rank
mt = m[:,t]
s = svd(mt);
plot(1340:1380,s.S[1340:1380],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="Size of singular values",xlabel="Singular values")
_,Rtest = qr(mt);
plot(abs.(diag(Rtest)),yaxis=:log)
plot(sort(abs.(diag(Rtest)),rev=true),yaxis=:log,st=:scatter,ms=2,msw=0,color=:purple)
plot([1370:1390],sort(abs.(diag(Rtest))[1370:1390],rev=true),yaxis=:log,st=:scatter,ms=2,msw=0)

# The singular values and rank don't closely match the size of entries in the R'S
# So let's keep some extras from the R's, but not all:
on = get_nonzero_index(Rtest,8e-13); count(on)
s = t[on]



popfirst!(s) # Don't need the first element -- it refers to the empty cluster, not listed in clusters.out
s .-= 1 # Subtract 1 from each entry because clusters array and pi matrix column index are shifted by 1
clust5li = clust5[s]
sv5li = sv5[s]
id5li =  id5[s]
dists5li = dists5[s]
nV5li = nV5[s]
write_clusters(clust5li, sv5li, nV5li, dists5li, id5li,"clusters.2-5b_li"*string(r))
#>  cp clusters.2b3b_li418 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-5b_li"*string(r))
rank(m)
rank(m[1:137,:])

# 6 body terms
pts = clusters.make_eqvPoints(Rpts[1:2],rots)
@time figs6bt = clusters.make_figure_candidates(pts,6) # with 5 rpts, 78 eqv points, 21 million cands, 15 mins
@time figs6btemp = reduce_figList(figs6bt[1:1_000_000],rots) # 200_000 took 5 mins
figs6b = figs6btemp[1:end] 
n6b = length(figs6b)
figs6b, sv6, nV6, dists6, id6 = read_clusters_from_file("clusters.2-6b_3232")
#mask = nV6.==6
#n6b = count(mask)
clust6 = vcat(clust5li,figs6b)
sv6 = vcat(sv5li,[[1,1,1,1,1,1] for i ∈ 1:n6b])
id6 =  vcat(id5li,1:n6b)
dists6 = vcat(dists5li,MADfromCOM.(figs6b))
nV6 = vcat(nV5li,ones(Int,n6b)*6)
# clust6=figs6b
# sv6 = [[1,1,1,1,1,1] for i ∈ 1:n6b]
# id6 =  1:n6b
# dists6 = MADfromCOM.(figs6b)
# nV6 = ones(Int,n6b)*6
write_clusters(clust6, sv6, nV6, dists6, id6,"clusters.2-6b_"*string(n6b))
##write_clusters(clust6, sv6, nV6, dists6, id6,"clusters.2-6b_"*string(n6b))
#>  cp clusters.2-6b_2600 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-6b_"*string(n6b))
rank(m)
rank(m[1:137,:]) #Still missing 2 columns, maybe missing lower order clusters...
rank(m[:,1:2346]) # 33 l.i. columns to the right of the diagonal

t = get_leftmost_indep_columns(m)

r=rank(m[:,t]) # Still have 1727 rank
mt = m[:,t]
s = svd(mt);
plot(1700:1750,s.S[1700:1750],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="Size of singular values",xlabel="Singular values")
_,Rtest = qr(mt);
plot(abs.(diag(Rtest)),yaxis=:log)
plot(sort(abs.(diag(Rtest)),rev=true),yaxis=:log,st=:scatter,ms=2,msw=0,color=:purple)
plot([1700:1750],sort(abs.(diag(Rtest))[1700:1750],rev=true),yaxis=:log,st=:scatter,ms=2,msw=0)

# The singular values and rank don't closely match the size of entries in the R'S
# So let's keep some extras from the R's, but not all:
on = get_nonzero_index(Rtest,1e-12); count(on)
s = copy(t)
s = t[on]

rank(m[:,s])

popfirst!(s) # Don't need the first element -- it refers to the empty cluster, not listed in clusters.out
s .-= 1 # Subtract 1 from each entry because clusters array and pi matrix column index are shifted by 1
clust6li = clust6[s]
sv6li = sv6[s]
id6li =  id6[s]
dists6li = dists6[s]
nV6li = nV6[s]
write_clusters(clust6li, sv6li, nV6li, dists6li, id6li,"clusters.2-6b_li"*string(r))
#>  cp clusters.2b3b_li418 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-6b_li"*string(r))
rank(m)
rank(m[1:137,:])


# Inserting new list of threes
clust6li, sv6li, nV6li, dists6li, id6li = read_clusters_from_file("clusters.2-6b_li1727")
mask = nV6li .!= 3
n3b = length(figs3b)
nV3p = nV6li[mask]
inst = findfirst(nV3p.==4)
temp = clust6li[mask] 
figs3b = [hcat(i) for i ∈ figs3b]
clust3p = vcat(temp[1:inst-1],figs3b,temp[inst:end])
sv6li = sv6li[mask]
svinst = [[1,1,1] for i ∈ 1:n3b]
sv3p = vcat(sv6li[1:inst-1],svinst,sv6li[inst:end])
nV6li = nV6li[mask]
nV3p = vcat(nV6li[1:inst-1],ones(Int,n3b)*3,nV6li[inst:end])
dists6li = dists6li[mask]
dists3p = vcat(dists6li[1:inst-1],MADfromCOM.(figs3b),dists6li[inst:end])
id6li = id6li[mask]
id3p = vcat(id6li[1:inst-1],1:n3b,id6li[inst:end])

write_clusters(clust3p, sv3p, nV3p, dists3p, id3p,"clusters.all3p_"*string(n3b))
#>  cp clusters.all3p_... clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix pimat.all3b_...
m = readdlm("pimat.all3b_"*string(n3b))
rank(m)
rank(m[1:137,:])

t = get_leftmost_indep_columns(m)

default(ms=2,msw=0,legend=:none,st=:scatter)
r=rank(m[:,t]) # Still have 1731 rank
mt = m[:,t]
s = svd(mt);
plot(1700:1750,s.S[1700:1750],msw=0,ms=2,st=:scatter,yaxis=:log,legend=:none,ylabel="Size of singular values",xlabel="Singular values")
_,Rtest = qr(mt);
plot(abs.(diag(Rtest)),yaxis=:log)
plot(sort(abs.(diag(Rtest)),rev=true),yaxis=:log,st=:scatter,ms=2,msw=0,color=:purple)
plot([1720:1783],sort(abs.(diag(Rtest)),rev=true)[1720:1783],yaxis=:log,st=:scatter,ms=2,msw=0)

# The singular values and rank don't closely match the size of entries in the R'S
# So let's keep some extras from the R's, but not all:
on = get_nonzero_index(Rtest,1e-12); count(on)
s = copy(t)
s = t[on]

rank(m[:,s])

popfirst!(s) # Don't need the first element -- it refers to the empty cluster, not listed in clusters.out
s .-= 1 # Subtract 1 from each entry because clusters array and pi matrix column index are shifted by 1
clust3pli = clust3p[s]
sv3pli = sv3p[s]
id3pli =  id3p[s]
dists3pli = dists3p[s]
nV3pli = nV3p[s]
write_clusters(clust3pli, sv3pli, nV3pli, dists3pli, id3pli,"clusters.2-6b_li"*string(r))
#>  cp clusters.2b3b_li418 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m6 = readdlm("pimat.2-6b_li"*string(r))
rank(m6)
rank(m6[1:137,:])
plot(m6[125,:],st=:scatter)
plot!(m6[126,:].+0.013,st=:scatter,color=:purple)
plot!(m6[127,:].+0.026,st=:scatter,color=:red)
plot!(m6[128,:].+0.039,st=:scatter,color=:green)