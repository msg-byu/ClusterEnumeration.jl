# Trying to do everything with ints for the sake of speed
using ClusterEnumeration
using Revise
using Spacey
using LinearAlgebra
using DelimitedFiles
using Plots

cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/pairs/")
default(ms=2,msw=0,legend=:none,st=:scatter)
A = [0 1 1;
     1 0 1;
     1 1 0]

_,rots = pointGroup_robust(eachrow(A)...)
rots = [Int.(i) for i in rots]
@time rpts, norms = genLatticePts(A,rots,40); length(rpts)
#writedlm("rpts.upto40cells.dat",rpts) # Took several minutes to generate
rpts = readdlm("rpts.upto40cells.dat")
rpts =  [Int.(i) for i in eachrow(rpts)]
rpts
norms
sphPts, degen, lengths =genSymEqvPoints(rpts,rots)
lengths
length(sphPts)
rpts[end-2:end].|>norm

# figCands = makeFigureCandidates(sphPts,2)
# rdFigs = symReduceFigList(figCands,rots)
# for n ∈ 7:10
# @time figs = makeFullFigsFromShellPts(rpts[1:n],rots,3)
# display(length(figs))
# end

#@time figs = makeFullFigsFromShellPts(rpts[1:11],rots,3)

# This generated ??? 3b figs in ??? secs using the first ??? shells
begin
n = 500
sphPts, degen, lengths =genSymEqvPoints(rpts[1:n],rots); @show length(sphPts)
figCands2 = makeFigureCandidates(sphPts,2); @show length(figCands2)
rdFigs2 = symReduceFigList(figCands2,rots); @show length(rdFigs2)
@time figCands3 = iterAugmentFigures(rdFigs2,sphPts); @show length(figCands3)
@time rdFigs3 = symReduceFigList(figCands3,rots); @show length(rdFigs3)
end

f3 = rdFigs3./2
# Should have divided the fig components by two, to match the lattice parameter in lat.in!
cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/pairs/")
cl, sv, nv, dists, id = read_clusters_from_file("clusters.2-6b_li1760");
n3b = length(f3)
cl3 = vcat(f3,cl)
sv3 = vcat([[1,1,1] for i ∈ 1:n3b],sv)
nv3 = vcat(ones(n3b)*3,nv)
dists3 = vcat(diameter.(f3),dists)
id3 = vcat(1:n3b,id) 
write_clusters(cl3,sv3,nv3,dists3,id3,"clusters.radGen.3")
m = readdlm("pimat.radGen.3")
rank(m)
rank(m[1:137,:])
t = get_leftmost_indep_columns(m)

# Rearrange 1760 dists to be in ascending nV order

dp = sortperm(dists;alg=Base.Sort.DEFAULT_STABLE)
nvp = sortperm(nv[dp];alg=Base.Sort.DEFAULT_STABLE)
begin
plot(dists[nvp[dp]],xlabel="Cluster number",ylabel="Diameter of cluster",title="Independent Clusters, 2346 structures (size ≤ 10)",st=:scatter,msw=0,color=:red,ms=1,legend=:none,)
annotate!(60, 4, text("Pairs", :blue, :left,12))
annotate!(200, 1.6, text("Triplets", :blue, :left,12))
annotate!(550, 2, text("4-bodies", :blue, :left,12))
annotate!(1000, 1.9, text("5-bodies", :blue, :left,12))
p1 = annotate!(1450, 0.9, text("6-bodies", :blue, :left,12))
savefig(p1,"allclusters_upto10.pdf")
end
write_clusters(cl[nvp[dp]],sv[nvp[dp]],nv[nvp[dp]],dists[nvp[dp]],id[nvp[dp]],"clusters.2-6b_li1760")

# Start inside 8679-th shell, radius is ≈56.6
sphPts, degen, lengths =genSymEqvPoints(rpts,rots); @show length(sphPts)
figCands2 = makeFigureCandidates(sphPts,2); @show length(figCands2)
rdFigs2 = symReduceFigList(figCands2,rots); @show length(rdFigs2)
# radius is ≈47.1, about shell 5112
#writedlm("uqPairsTo5112.dat",rdFigs2) # Gave me a rank of 31
#rdFigs2 = readdlm("uqPairsTo5112.dat")
#rdFigs2 = [reshape(Int.(i),(3,2)) for i in eachrow(rdFigs2)]
@time figCands3 = iterAugmentFigures(rdFigs2[1:75],sphPts[1:3000]); @show length(figCands3)
plot(diameter.(rdFigs2[1:75]))
plot(norm.(sphPts[1:3000]))
# Makes new clusters inside a sphere of ≈10, about 75th shell, 113K candidates
# On big search, max l.i. 3-body had length ≈5
@time rdFigs3 = symReduceFigList(figCands3,rots); @show length(rdFigs3)
writedlm("rdFigs3from_75_3000.dat",rdFigs3)
# rdFigs3 = readdlm("rdFigs3from_75_3000.dat")
# rdFigs3 = [reshape(Int.(i),(3,2)) for i in eachrow(rdFigs2)]


@time figCands4 = iterAugmentFigures(rdFigs3[1:500],sphPts); @show length(figCands4)
@time rdFigs4 = symReduceFigList(figCands4,rots); @show length(rdFigs4)
@time figCands5 = iterAugmentFigures(rdFigs4[1:800],sphPts); @show length(figCands5)
@time rdFigs5 = symReduceFigList(figCands5,rots); @show length(rdFigs5)
@time figCands6 = iterAugmentFigures(rdFigs5[1:1000],sphPts); @show length(figCands6)
@time rdFigs6 = symReduceFigList(figCands6,rots); @show length(rdFigs6)
end

cl, sv, nv, dists, id = read_clusters_from_file("pointClusterOnly.out");
#write_clusters(cl[1:2], sv[1:2], nv[1:2], dists[1:2], id[1:2],  "pointClusterOnly.out")
f2 = rdFigs2./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n2b = length(f2)
cl2 = vcat(cl,f2)
sv2 = vcat(sv,[[1,1,1] for i ∈ 1:n2b])
nv2 = vcat(nv,ones(n2b)*2)
dists2 = vcat(dists,diameter.(f2))
id2 = vcat(id,1:n2b) 
write_clusters(cl2,sv2,nv2,dists2,id2,"clusters.radGen.2")
# re-run uncle 44
m = readdlm("pimat.radGen.2")
rank(m) # rank 31, 29 independent pairs
t = get_leftmost_indep_columns(m,1500)
begin
nP = 50; c = [i ∈ t ? :red : :blue for i ∈ 1:nP]; sz = [i ∈ t ? 5 : 6 for i ∈ 1:nP]
sh = [i∈t ? :circle : :star4 for i ∈ 1:nP]
p1 = plot(dists2[1:nP],ms=sz,shape=sh,xlabel="Pair number",ylabel="Pair length",color=c,title="Independent Pair Clusters",label="Independent",legend=:right)
annotate!(30,1.2,"Longest pair: ≈23.57\n Total pairs: 5128")
#savefig(p1,"pairs.pdf")
end
dists2[end] # ≈23.57
length(dists2)

rank(m[:,t])
t
popfirst!(t)
t.-= 1

write_clusters(cl2[t],sv2[t],nv2[t],dists2[t],id2[t],"clusters.radGen.2li_31")
# rerun uncle...
m = readdlm("pimat.radGen.2li_31")
rank(m)

cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.2li_31");
f3 = rdFigs3./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n3b = length(f3)
cl3 = vcat(cl,f3)
sv3 = vcat(sv,[[1,1,1] for i ∈ 1:n3b])
nv3 = vcat(nv,ones(n3b)*3)
dists3 = vcat(dists,diameter.(f3))
id3 = vcat(id,1:n3b) 
write_clusters(cl3,sv3,nv3,dists3,id3,"clusters.radGen.3")
# re-run uncle 44
m = readdlm("pimat.radGen.3")
rank(m) # rank 70, 70-31=39 three-bodies independent pairs

t = get_leftmost_indep_columns(m)
rank(m[:,t])

begin
     intv = 32:220
     c = [i ∈ t ? :red : :black for i ∈ intv]; sz = [i ∈ t ? 4 : 1 for i ∈ intv]
     sh = [i∈t ? :circle : :circle for i ∈ intv]
     p2 = plot(dists3[intv],ms=sz,shape=sh,xlabel="Triplet Number",ylabel="Triplet Diameter",color=c,title="Independent Triplet Clusters",label="Independent",legend=:right)
     annotate!(100,1.2,"Max. diameter searched: ≈25.50\n Total triplets searched: 15,270")
     savefig(p2,"triplets.pdf")
end
f3[end]/.2|>diameter
t[end]-32 # Last l.i. 3-body cluster is 167
dists3[199] # Diameter is 2.32065

popfirst!(t)
t.-= 1
write_clusters(cl3[t],sv3[t],nv3[t],dists3[t],id3[t],"clusters.radGen.3li_70")
# rerun uncle...
m = readdlm("pimat.radGen.3li_70")
rank(m)


cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.3li_70");
f3 = rdFigs3./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n3b = length(f3)
cl3 = vcat(cl,f3)
sv3 = vcat(sv,[[1,1,1] for i ∈ 1:n3b])
nv3 = vcat(nv,ones(n3b)*3)
dists3 = vcat(dists,diameter.(f3))
id3 = vcat(id,1:n3b) 
write_clusters(cl3,sv3,nv3,dists3,id3,"clusters.radGen.3")
# re-run uncle 44
m = readdlm("pimat.radGen.3")
rank(m) # rank 70, 70-31=39 three-bodies independent pairs

t = get_leftmost_indep_columns(m)
rank(m[:,t])

popfirst!(t)
t.-= 1
write_clusters(cl3[t],sv3[t],nv3[t],dists3[t],id3[t],"clusters.radGen.3li_70")
# rerun uncle...
m = readdlm("pimat.radGen.3li_70")
rank(m)


plot(diameter.(rdFigs3[1:300]))
plot(norm.(sphPts[1:320]))
@time figCands4 = iterAugmentFigures(rdFigs3[1:300],sphPts[1:320]); @show length(figCands4)
# Makes new clusters inside a sphere of ≈6, about 22th shell, 113K candidates
# On big search, max l.i. 4-body had length ⪅3
# Too big, try a cutoff of 5 => 44K cands
@time rdFigs4 = symReduceFigList(figCands4,rots); @show length(rdFigs4)
writedlm("rdFigs4from_300_320.dat",rdFigs4)
# rdFigs4 = readdlm("rdFigs3from_75_3000.dat")
# rdFigs4 = [reshape(Int.(i),(3,2)) for i in eachrow(rdFigs2)]

cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.3li_70");
f4 = rdFigs4./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n4b = length(f4)
cl4 = vcat(cl,f4)
sv4 = vcat(sv,[[1,1,1,1] for i ∈ 1:n4b])
nv4 = vcat(nv,ones(n4b)*4)
dists4 = vcat(dists,diameter.(f4))
id4 = vcat(id,1:n4b) 
write_clusters(cl4,sv4,nv4,dists4,id4,"clusters.radGen.4")
# re-run uncle 44
m = readdlm("pimat.radGen.4")
rank(m) # rank 112, 112-70=42 independent four-bodies 

t = get_leftmost_indep_columns(m,500)
rank(m[:,t])

begin
     intv = 72:750
     c = [i ∈ t ? :red : :black for i ∈ intv]; sz = [i ∈ t ? 4 : 1 for i ∈ intv]
     sh = [i∈t ? :circle : :circle for i ∈ intv]
     sw = [i∈t ? 0 : 1 for i ∈ intv]
     p3 = plot(dists4[intv],ms=sz,shape=sh,msw=sw,xlabel="4-body Number",ylabel="Diameter",color=c,title="Independent 4-body Clusters",label="Independent",legend=:right)
     annotate!(400,1.15,"Max. diameter searched: ≈2.56\n Total 4-bodies searched: 11,883")
     savefig(p3,"quads.pdf")
end
t[end]-72 # 632-nd quad is last l.i. one
dists4[704] # Diameter is 1.833
dists4[end] # 2.56 largest diameter searched

popfirst!(t)
t.-= 1
write_clusters(cl4[t],sv4[t],nv4[t],dists4[t],id4[t],"clusters.radGen.4li_112")
# rerun uncle...
m = readdlm("pimat.radGen.4li_112")
rank(m)

plot(diameter.(rdFigs4[1:1500]))
plot(norm.(sphPts[1:140]))
@time figCands5 = iterAugmentFigures(rdFigs4[1:1500],sphPts[1:140]); @show length(figCands5)
# Makes new clusters inside a sphere of ≈4, about 8th shell, 89K candidates
# On big search, max l.i. 5-body had length ⪅2
@time rdFigs5 = symReduceFigList(figCands5,rots); @show length(rdFigs5)
#writedlm("rdFigs5from_1500_140.dat",rdFigs5)
rdFigs5 = readdlm("rdFigs5from_1500_140.dat")
rdFigs5 = [reshape(Int.(i),(3,5)) for i in eachrow(rdFigs5)]

cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.4li_112");
f5 = rdFigs5./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n5b = length(f5)
cl5 = vcat(cl,f5)
sv5 = vcat(sv,[[1,1,1,1,1] for i ∈ 1:n5b])
nv5 = vcat(nv,ones(n5b)*5)
dists5 = vcat(dists,diameter.(f5))
id5 = vcat(id,1:n5b) 
write_clusters(cl5,sv5,nv5,dists5,id5,"clusters.radGen.5")
# re-run uncle 44
m = readdlm("pimat.radGen.5")
rank(m) # rank 127, 127-112=15 independent five-bodies 

t5 = get_leftmost_indep_columns(m,500)
rank(m[:,t])

begin
     intv = 114:5000
     c = [i ∈ t5 ? :red : :black for i ∈ intv]; sz = [i ∈ t5 ? 7 : 1 for i ∈ intv]
     sh = [i∈t5 ? :circle : :auto for i ∈ intv]
     sw = [i∈t5 ? 0 : 1 for i ∈ intv]
     p4 = plot(dists5[intv],ms=sz,shape=sh,msw=sw,xlabel="5-body number",ylabel="Diameter",color=c,title="Independent 5-body Clusters",label="Independent",legend=:right)
     annotate!(3000,1.1,"Max. diameter searched: ≈2.02\n Total 5-bodies searched: 23,066")
     savefig(p4,"quints.pdf")
end
dists5[end]
t5[end]
4599-112 # 4487-th quint is last l.i. one
dists5[4559] # diameter is ≈1.769

popfirst!(t)
t.-= 1
write_clusters(cl5[t],sv5[t],nv5[t],dists5[t],id5[t],"clusters.radGen.5li_127")
# rerun uncle...
m = readdlm("pimat.radGen.5li_127")
rank(m)

###
### 6 bodies!
###

plot(diameter.(rdFigs5),yticks=1:0.25:4,gridalpha=0.5,grid=:all)
plot(norm.(sphPts[1:480]),yticks=1:0.25:4.0)
#@time figCands6 = iterAugmentFigures(rdFigs5[1:9500],sphPts[1:134]); @show length(figCands6)
@time figCands6 = iterAugmentFigures(rdFigs5,sphPts[1:248]); @show length(figCands6)
norms = diameter.(figCands6./2)
count(norms.<1.85) #47874 #52830 #52968 #53042
idx = findall(norms.<1.60)
plot(norms)
histogram(norms)

@time rdFigs6 = symReduceByLengthFigList(figCands6[findall(norms.<1.85)],rots);
plot(diameter.(rdFigs6))
symReduceByLengthFigList(figCands6[findall(norms.<1.45)],rots);
# Makes new clusters inside a sphere of ≈3.75, 7th shell for shpPts, 126K candidates
# On big search, max l.i. 6-body had length ⪅1.7
@time rdFigs6 = symReduceFigList(figCands6[findall(norms.<1.6)],rots); @show length(rdFigs6)
#writedlm("rdFigs6from_3900_86.dat",rdFigs6)
#writedlm("quickDirty_figs6.dat",rdFigs6)
#rdFigs6 = readdlm("quickDirty_figs6.dat")
#rdFigs6 = [reshape(Int.(i),(3,6)) for i in eachrow(rdFigs6)]

cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.5li_127");
f6 = rdFigs6./2
# Should divide the fig components by two, to match the lattice parameter in lat.in!
n6b = length(f6)
cl6 = vcat(cl,f6)
sv6 = vcat(sv,[[1,1,1,1,1,1] for i ∈ 1:n6b])
nv6 = vcat(nv,ones(n6b)*6)
dists6 = vcat(dists,diameter.(f6))
id6 = vcat(id,1:n6b) 
write_clusters(cl6,sv6,nv6,dists6,id6,"clusters.radGen.qd6")
# re-run uncle 44
m = readdlm("pimat.radGen.qd6")
rank(m) # rank 137! Finally!

# Need a better algorithm for this now.
#t = get_leftmost_indep_columns(m,500)
r = rank(m)

failsafe = 1
next = size(m,2)
next =  8000
curr = 136
l6 = [58949]
while true
     _,R,pv1 = qr(m[:,1:next], ColumnNorm())
     next = maximum(pv1[1:137])-1
     r=rank(m[:,1:next])
     if r < curr push!(l6,next+1); curr=r;println(l6) end
     println("Next left: ", next,"   rank: ",r ,"   Iteration: ",failsafe)
     #if failsafe > 126 break end
     if r < 127 break end
     failsafe +=1
end
rank(m[:,1:7816])
rank(m[:,1:7815])
_,R,pv1 = qr(m[:,1:7816], ColumnNorm())

#t6 = collect(1:128)
t6 = append!(collect(1:128),[129,130,131,132,139,147,150,173,227,7817,58949])

rank(m[:,t6])
size(f6)
dists6[end]
begin
     intv = 129:60_000
     c = [i ∈ t6 ? :red : :black for i ∈ intv]; sz = [i ∈ t6 ? 5 : 1 for i ∈ intv]
     sh = [i∈t6 ? :circle : :circle for i ∈ intv]
     sw = [i∈t6 ? 0 : 1 for i ∈ intv]
     p5 = plot(dists6[intv],ms=sz,st=:scatter,shape=sh,msw=sw,xlabel="6-body Number",ylabel="Diameter",color=c,title="Independent 6-body Clusters",label="Independent",legend=:right)
     annotate!(35000,1.1,"Max. diameter searched: ≈1.85\n Total 6-bodies searched: 76,040")
     #savefig(p5,"sextuplets.pdf")
end

popfirst!(t6)
t6.-= 1
write_clusters(cl6[t6],sv6[t6],nv6[t6],dists6[t6],id6[t6],"clusters.radGen.6li_137")
# rerun uncle...
m = readdlm("pimat.radGen.6li_137")
rank(m)

cl, sv, nv, dists, id = read_clusters_from_file("clusters.radGen.6li_136");


begin
p6=heatmap(hcat(m[:,:]),aspect_ratio=1,xlim=(0.5,139.5),ylim=(0.5,137.5),yticks=130:-10:1,xticks=0:10:137,xlabel="Clusters",ylabel="Structure number",legend=:true,yflip=true,dpi=400,tick_direction=:out,framestyle=:box,xmirror=true)#,color=:viridis)
w=.3

# plot!([0,139],[11,11],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([0,139],[29.5,29.5],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([0,139],[57.5,57.5],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([31.5,31.5],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([70.3,70.3],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([129,129],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
# plot!([114.2,114.2],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
savefig(p6,"pimatrix.png")
end


Different lattices: bcc, sc, hcp, fct
sum rows, columns of pimatrix, look at distributions (not at sum, but at devations)
list rows in different order (by cell shape instead of concentration)
list rows entirely by concentration
QR decomposition
SVD decomposition
Fourier transformation?
Reorder columns by cluster length
Reorder rows by cell edge length
Cell "diameter" (orthogonality defect)

All the DFT data related questions
k-nary cases
s-vector degeneracies

