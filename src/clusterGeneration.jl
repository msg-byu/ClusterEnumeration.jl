# If I import a module locally, I need to "use" the packages that the module needs
using Printf
using Spacey
using LinearAlgebra
using Plots
using DelimitedFiles

default(ms=2,msw=0,legend=:none,st=:scatter)
include("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/find_clusters.jl")
using .clusters
cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/")

cd("pairs")

A = read_lattice_vectors()
_, rots = yowhGroup_robust(eachcol(A)...) 
@time Rpts = clusters.gen_points_in_supercell(A,1000) # At 5000, takes a long time...10-15 mins
writedlm("pairs_fcc_15_000.dat",Rpts)
R1 = readdlm("pairs_fcc_10_000.dat")
pts = make_eqvPoints(Rpts[1:6],rots)pfigs6b = clusters.make_figure_candidates(pts,6)
@time s1 = make_figure_candidates(make_eqvPoints(Rpts[1:4],rots),6) #3.1M 15 secs
@time s2 = make_figure_candidates(make_eqvPoints(Rpts[1:5],rots),6) #21M 150 secs
t1 = round.(MADfromCOM.(s1),digits=12)|>unique|>sort
t2 = round.(MADfromCOM.(s2),digits=12)|>unique|>sort
# t1 has 32K unique diameters

# Two long pair enumerations were identical matches to 14765 pairs. readdlm("pairs_fcc_14765.dat")
plot([norm(i) for i ∈ eachrow(R1)],title="Pairs",xlabel="Pair number",ylabel="Length")
plot!(norm.(Rpts[1:size(R1,1)]),title="Pairs",xlabel="Pair number",ylabel="Length",label="Longer enumeration")
plot([norm(i) for i ∈ eachrow(R1)] .- norm.(Rpts[1:size(R1,1)]) |>x-> plot([14000:15000],x[14000:15000]),title="Length difference")



plot(t1[1:1500]-t2[1:1500],title="Difference in diameters for 4NN and 5NN 6b sets",xlabel = "Diameter number",ylabel = "Difference")
findfirst(x->x>0,t1[1:1500]-t2[1:1500]) # First different diameter occurs at 341.

@time v1 = make_figure_candidates(make_eqvPoints(Rpts[1:6],rots),5) # 5: 1.4M 20 secs. 6: 2.1M, 28 secs
@time v2 = make_figure_candidates(make_eqvPoints(Rpts[1:7],rots),5) # 7: 12.8M 140 secs
s1 = round.(MADfromCOM.(v1),digits=12)|>unique|>sort
s2 = round.(MADfromCOM.(v2),digits=12)|>unique|>sort
plot(s1[1:1500]-s2[1:1500],title="Difference in diameters for 6NN and 7NN 5b sets",xlabel = "Diameter number",ylabel = "Difference")
findfirst(x->x>0,s1[1:1500]-s2[1:1500]) # First different diameter occurs at 725.

b2 = Vector{Matrix{Float64}}()
f2 = Vector{Matrix{Float64}}()
#for np = 2:4
np = 50 # Roughly the number of "shells", n-th neighbor
    pts = make_eqvPoints(Rpts[1:np],rots)
    b2 = make_figure_candidates(pts,2)
    f2 = reduce_figList(b2,rots)
    # Now add
#end
#b3 = Vector{Matrix{Float64}}()
b3 =  [sortslices(hcat(j,i),dims=2) for i in pts for j in f2 if !any(sum(j.-i.≈zeros(3),dims=1).==3)] # Add a third point to each pair, but don't add the same point twice   
@time f3 = reduce_figList(b3,rots);length(f3)
@time r3 = reduce_figList(make_figure_candidates(make_eqvPoints(Rpts[1:np],rots),3),rots);length(r3)

# up to reduced pools of 1196 (17 pairs), the two approaches give identical figures (in slightly different orders)

function augment_clusters(clust,shellpts)
    newclust = [sortslices(hcat(j,i),dims=2) for i in shellpts for j in clust if !any(sum(j.-i.≈zeros(3),dims=1).==3)]
    return newclust
end

b4 = augment_clusters(f3,pts)
f4 = reduce_figList(b4,rots)
b5 = augment_clusters(f4,pts)
b5p = b5[findall(x<2,MADfromCOM.(b5))]
f5 = reduce_figList(b5p,rots)
b6 = augment_clusters(f5,pts)
b6p = b6[findall(x->x<1.7,MADfromCOM.(b6))]
f6 = reduce_figList(b6p,rots)

n6b = length(f6)
clust, sv, nV, dists, id = read_clusters_from_file("clusters.2-5b_li1364")
clust6 = vcat(clust,f6)
sv6 = vcat(sv,[[1,1,1,1,1,1] for i ∈ 1:n6b])
nV6 = vcat(nV,ones(n6b)*6)
dists6 = vcat(dists,MADfromCOM.(f6))
id6 = vcat(id,1:n6b) 
write_clusters(clust6,sv6,nV6,dists6,id6,"clusters.out.new6")
# cp clusters.out.new6 clusters.out 
# rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.new6
m = readdlm("pimat.new6");
rank(m)
rank(m[1:137,:])
t = get_leftmost_indep_columns(m)
mt=m[:,t]
_,R = qr(mt);
plot([1700:1750],sort(abs.(diag(R)),rev=true)[1700:1750],yaxis=:log)
s = svd(mt);
plot([1700:1750],s.S[1700:1750],yaxis=:log,color=:red)
r = get_nonzero_index(R,1e-4)
u = t[r]
mt=m[:,u] 
rank(mt)

popfirst!(u) # First cluster is the empty cluster and not in the file, so don't include it explicitly
u .-= 1 # Shift the indexing
write_clusters(clust6[u],sv6[u],nV6[u],dists6[u],id6[u],"clusters.2-6b_li1736")
mr = readdlm("pimat.2-6b_li1736");
rank(mr)

# Try a bigger (and less "holey"?) 3-body pool (more than 27,000)
n3b = length(f3)
clust, sv, nV, dists, id = read_clusters_from_file("clusters.2-6b_li1736")
clust3 = vcat(clust,f3)
sv3 = vcat(sv,[[1,1,1] for i ∈ 1:n3b])
nV3 = vcat(nV,ones(n3b)*3)
dists3 = vcat(dists,MADfromCOM.(f3))
id3 = vcat(id,1:n3b) 
write_clusters(clust3,sv3,nV3,dists3,id3,"clusters.out.new3")
# cp clusters.out.new3 clusters.out 
# rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.new3
m = readdlm("pimat.new3");
rank(m)
rank(m[1:137,:])
t = get_leftmost_indep_columns(m)
mt=m[:,t]
_,R = qr(mt);
plot([1700:1750],sort(abs.(diag(R)),rev=true)[1700:1750],yaxis=:log)
s = svd(mt);
plot([1700:1750],s.S[1700:1750],yaxis=:log,color=:red)
r = get_nonzero_index(R,1e-4)
u = t[r]
mt=m[:,u] 
rank(mt)

popfirst!(u) # First cluster is the empty cluster and not in the file, so don't include it explicitly
u .-= 1 # Shift the indexing
write_clusters(clust3[u],sv3[u],nV3[u],dists3[u],id3[u],"clusters.2-6b_li1741")
mr = readdlm("pimat.2-6b_li1741");
rank(mr)

# Try 14765 pairs. This didn't add even 1 to the rank.
pairs = readdlm("pairs_fcc_14765.dat")
pairs = [hcat(i...) for i in eachrow(pairs) ]
pairs = [hcat([0.,0.,0.],vcat(ipair...)) for ipair ∈ pairs]
n2b = length(pairs)
clust, sv, nV, dists, id = read_clusters_from_file("clusters.2-6b_li1736")
clust2 = vcat(clust,pairs)
sv2 = vcat(sv,[[1,1] for i ∈ 1:n2b])
nV2 = vcat(nV,ones(n2b)*2)
dists2 = vcat(dists,MADfromCOM.(pairs))
id2 = vcat(id,1:n2b) 
write_clusters(clust2,sv2,nV2,dists2,id2,"clusters.out.new2")
# cp clusters.out.new2 clusters.out 
# rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.new2
m = readdlm("pimat.new2");
rank(m) # 1736 --- going out to 14765 pairs did NOT add anything to the rank. So the last pair probably really is 356 (137 pairs)
rank(m[1:137,:])

plot(m[129,1:1736])
# The 128 structure doesn't seem to have anything unusual about its correlations. the are mostly multiples of 1/3 instead of 1/6 but nothing that really stands out compared to the other structures in that block

np = 50 # Roughly the number of "shells", n-th neighbor
pts = make_eqvPoints(Rpts[1:np],rots)
b2 = make_figure_candidates(pts,2)
f2 = reduce_figList(b2,rots)
b3 =  [sortslices(hcat(j,i),dims=2) for i in pts for j in f2 if !any(sum(j.-i.≈zeros(3),dims=1).==3)] # Add a third point to each pair, but don't add the same point twice   
@time f3 = reduce_figList(b3,rots);length(f3)
# Made 14515 three bodies in 7 mins, using 50 shells 
plot(MADfromCOM.(f3))

# Find a big, but not too big pool of 4 bodies
MADfromCOM(f3[200]) # What is the distance of the n-th three body
# And how many pts to keep at that same distance
findfirst(norm.(eachrow(pts)) .> 2.738)

@time b4 = augment_clusters(f3[1:350],pts[1:368]);length(b4) #248 is the end of 13-th NN? 
@time f4_3 = reduce_figList(b4,rots); length(f4_3)
plot(MADfromCOM.(f4_3)) # Holes may appear around 5000, so enumerate some meore

clust, sv, nV, dists, id = read_clusters_from_file("clusters.2-6b_li1741")
n4b = 20_000
f4 = f4_3[1:n4b]
clust4 = vcat(clust,f4)
sv4 = vcat(sv,[[1,1,1,1] for i ∈ 1:n4b])
nV4 = vcat(nV,ones(n4b)*4)
dists4 = vcat(dists,MADfromCOM.(f4))
id4 = vcat(id,1:n4b) 
write_clusters(clust4,sv4,nV4,dists4,id4,"clusters.out.new4")
# cp clusters.out.new4 clust
# rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.new4
m = readdlm("pimat.new4");
rank(m) #1743 with 20_000 new 4-bodies, last l.i. 4-body at 11582
rank(m[1:137,:])

t = get_leftmost_indep_columns(m)
mt=m[:,t]
_,R = qr(mt);
plot([1700:1750],sort(abs.(diag(R)),rev=true)[1700:1750],yaxis=:log)
s = svd(mt);
plot([1700:1750],s.S[1700:1750],yaxis=:log,color=:red)
r = get_nonzero_index(R,1e-5);count(r)
u = t[r]
mt=m[:,u] 
rank(mt)

popfirst!(u) # First cluster is the empty cluster and not in the file, so don't include it explicitly
u .-= 1 # Shift the indexing
write_clusters(clust4[u],sv4[u],nV4[u],dists4[u],id4[u],"clusters.2-6b_li1743")
mr = readdlm("pimat.2-6b_li1743");
rank(mr)

@time b5 = augment_clusters(f4[1:350],pts[1:368]);length(b5) 
@time f5 = reduce_figList(b5,rots); length(f5)


clust, sv, nV, dists, id = read_clusters_from_file("clusters.2-6b_li1743")
n5b = 20_000
clust5 = vcat(clust,f5[1:n5b])
sv5 = vcat(sv,[[1,1,1,1,1] for i ∈ 1:n5b])
nV5 = vcat(nV,ones(n5b)*5)
dists5 = vcat(dists,MADfromCOM.(f5[1:n5b]))
id5 = vcat(id,1:n5b) 
write_clusters(clust5,sv5,nV5,dists5,id5,"clusters.out.new5")
# cp clusters.out.new5 clusters.out
# rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44; cp enum_PI_matrix.out pimat.new5
m = readdlm("pimat.new5");
rank(m) #1760 with 20_000 new 5-bodies, last l.i. 5-body at ?????
rank(m[1:137,:])

t = get_leftmost_indep_columns(m)
mt=m[:,t]
_,R = qr(mt);
plot([1740:1770],sort(abs.(diag(R)),rev=true)[1740:1770],yaxis=:log)
s = svd(mt);
plot([1740:1770],s.S[1740:1770],yaxis=:log,color=:red)
r = get_nonzero_index(R,1e-5);count(r)
u = t[r]
mt=m[:,u] 
rank(mt)

popfirst!(u) # First cluster is the empty cluster and not in the file, so don't include it explicitly
u .-= 1 # Shift the indexing
write_clusters(clust5[u],sv5[u],nV5[u],dists5[u],id5[u],"clusters.2-6b_li1760")
mr = readdlm("pimat.2-6b_li1760");
rank(mr)

# One more try at the three-body terms. Get a more "compact" set, no holes by taking more shells 
np = 110 # Roughly the number of "shells", n-th neighbor
pts = make_eqvPoints(Rpts[1:np],rots)
b2 = make_figure_candidates(pts,2)
f2 = reduce_figList(b2,rots)
b3 =  augment_clusters(f2,pts)
@time f3 = reduce_figList(b3,rots);length(f3)
# Made 14515 three bodies in 7 mins, using 50 shells 
# Made ????  3-bodies in ?? mins, using 110 shells
plot(MADfromCOM.(f3))





using StatsPlots
li2 = [2,5,10,19,29,48,69 ,100,135]
li3 = [0,3,10,19,39,65,120,189,288]
li4 = [0,0, 7,12,42,68,179,382,514]
li5 = [0,0, 0, 5,15,30,298,207,424]
li6 = [0,0, 0, 0, 9,14,258,151,378]
rd  = [0,0, 0, 0, 1,14,75 ,206,605]
groupedbar([li2 li3 li4 li5 li6 rd],bar_position=:stack)
groupedbar(replace!([li2 li3 li4 li5 li6 rd],0=>1),bar_position=:stack,yaxis=:log,ylims=(1,2400),xticks=(1:9,2:10))

plot(2:10,li2,yaxis=:log,label="pairs")
plot!(2:10,li3,yaxis=:log,ylims=(1,100))