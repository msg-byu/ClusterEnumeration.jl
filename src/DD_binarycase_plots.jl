# Used this code to generate things for the DD for binary fcc
#module clusters
using Plots
using Combinatorics
using Spacey
using LinearAlgebra
using DelimitedFiles
using Printf
using MinkowskiReduction
using Random
#using SMatrix

"""
    read_clusters_from_file(filename="clusters.out")

 Read in clusters info from the 'clusters.out' file 
    clusters, svecs, Verts, avgDist, IDnum = read_clusters_from_file("[file]")
"""
function read_clusters_from_file(filename="clusters.out")
    f = readlines(filename)
    vNidx = findall(x-> occursin("# Number of vertices:",x),f).+1
    vCount = parse.(Int,f[vNidx])
    avgDist = parse.(Float64,f[vNidx.+2])
    id = parse.(Int,f[vNidx.-2])
    clusters = []
    svecs = []    
    vNidx = vNidx.+4
    for i ∈ 1:length(vNidx)
        lines = [split(f[vNidx[i]+j]) for j ∈ 0:vCount[i]-1]
        v = [parse.(Float64,i[1:3]) for i ∈ lines]
        s = [parse.(Int,i[5]) for i ∈ lines]
        push!(clusters,v)
        push!(svecs,s)
    end
    clusters = [hcat(i...) for i ∈ clusters]
    return clusters,svecs, vCount, avgDist, id
end


cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data")
cl, sv, nV, d, id = read_clusters_from_file("clusters.2-6b_li1761")

colors = replace(nV, 1=>:red, 2=>:orange,3=>:green,4=>:magenta,5=>:blue,6=>:red)
pushfirst!(colors,:black)

m = readdlm("pimat.2-6b_li1761")
rank(m)
en = readdlm("energiesPerAtom.fcc.upto10.urgr")
J = m\en
plot(J,color=colors,st=:scatter,msw=0,ms=2,title="Binary Case, Size 1--10",xlabel="Cluster Number",ylabel="Coefficient (ECI)",legend=:none)
plot(abs.(J),color=colors,st=:scatter,msw=0,ms=2,title="Binary Case, Size 1--10",xlabel="Cluster Number",ylabel="Coefficient (ECI)",legend=:none,yaxis=:log)

fitErr = map(1:size(m,2)) do ms
    println(ms)
    J = m[:,1:ms]\en 
    norm(m[:,1:ms]*J-en)/sqrt(2346)
end
plot(1:size(m,2),fitErr,yaxis=:log,xlabel="Model Size",ylabel="Fitting Error (eV/atom)",title="Binary Case, Size 1--10, Simple Fit",st=:scatter,msw=0,ms=2,color=colors,legend=:none,yrange=(0,.03),xticks=0:100:1800)

# Definitely we want to try ordering by radius, irrespective of vertex order
begin
Nits = 100
sizes = collect(1:2:650) 
data = map(sizes) do ms
    eFit = 0
    eVal =0
    rav = 0
    println(ms)
    for i = 1:Nits 
        t = randperm(size(m,1))
        nFit = 300
        fitIdx = t[1:nFit]
        valIdx = t[nFit+1:end]
        J = m[fitIdx,1:ms]\en[fitIdx]
        eFit += norm(m[fitIdx,1:ms]*J-en[fitIdx])/sqrt(nFit)
        rav += rank(m[fitIdx,1:ms])
        eVal += norm(m[valIdx,1:ms]*J-en[valIdx])/sqrt(size(m,1)-nFit)
    end
    eFit/Nits,rav/Nits, eVal/Nits
    end 
    errFit = [i[1] for i ∈ data]
    ranks = [i[2] for i ∈ data]
    errVal = [i[3] for i ∈ data]
end

plot(sizes,errVal,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Validation Error", title = "Fitting to 300 Random Structures (Pt-Cu, ≥ 10)",legend=:none,color=colors[sizes],yrange=(1e-2,8e-0)
#,xticks = 0:50:1100,xrange=(0,600)
)

plot(sizes,errFit,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Fitting Error", title = "Fitting to 300 Random Structures (Pt-Cu, ≥ 10)",legend=:none,color=colors[sizes],yrange=(1e-4,2e-2)
#,xticks = 0:50:1100,xrange=(0,600)
)

begin
plot(sizes,ranks,legend=:none,title="Rank vs Model Size (Pt-Cu, ≥ 10, 2346 structures)",xlabel="Model size (num parameters)",ylabel="Rank of Designm Matrix",xrange=(0,400))
plot!([0; 300],[0; 300])
end 

clIdx = [1:10; 137:144; 431:439; 11:50; 145:430; 51:136;  440:500; 947:1000; 501:946; 1388:1400; 1001:1387; 1401:1761 ]

fitErr = map(1:size(m,2)) do ms
    println(ms)
    J = m[:,clIdx][:,1:ms]\en 
    norm(m[:,clIdx][:,1:ms]*J-en)/sqrt(2346)
end
plot(1:size(m,2),fitErr,yaxis=:identity,xlabel="Model Size",ylabel="Fitting Error (eV/atom)",title="Binary Case, Size 1--10, Simple Fit, Rearranged Clusters",st=:scatter,msw=0,ms=2,color=colors[clIdx],legend=:none,yrange=(3e-5,.003),xticks=0:100:1761,xrange=(1,1761))

mr = m[:,clIdx]
begin
Nits = 100
sizes = collect(1:2:650) 
data = map(sizes) do ms
    eFit = 0
    eVal =0
    rav = 0
    println(ms)
    for i = 1:Nits 
        t = randperm(size(m,1))
        nFit = 300
        fitIdx = t[1:nFit]
        valIdx = t[nFit+1:end]
        J = mr[fitIdx,1:ms]\en[fitIdx]
        eFit += norm(mr[fitIdx,1:ms]*J-en[fitIdx])/sqrt(nFit)
        rav += rank(mr[fitIdx,1:ms])
        eVal += norm(mr[valIdx,1:ms]*J-en[valIdx])/sqrt(size(m,1)-nFit)
    end
    eFit/Nits,rav/Nits, eVal/Nits
    end 
    errFit = [i[1] for i ∈ data]
    ranks = [i[2] for i ∈ data]
    errVal = [i[3] for i ∈ data]
end
    
plot(sizes,errVal,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Validation Error", title = "Fitting to 300 Random Structures (Pt-Cu, ≥ 10)",legend=:none,yrange=(1e-3,8e-0),color=colors[clIdx][sizes]
#,xticks = 0:50:1100,xrange=(0,600)
)

plot(sizes,errFit,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Fitting Error", title = "Fitting to 300 Random Structures (Pt-Cu, ≥ 10)",legend=:none,color=colors[sizes],yrange=(1e-4,2e-2)
#,xticks = 0:50:1100,xrange=(0,600)
)

begin
    plot(sizes,ranks,legend=:none,title="Rank vs Model Size (Pt-Cu, ≥ 10, 2346 structures)",xlabel="Model size (num parameters)",ylabel="Rank of Designm Matrix",xrange=(0,400))
    plot!([0; 300],[0; 300])
    end 
       
2==2






















 