# Try and re-sort the clusters in the pi matrix to maximize design matrix rank as clusters are added
using DelimitedFiles
using LinearAlgebra
using Plots
using Random
using ClusterEnumeration

### Moved these into ClusterEnumeration
### begin
""" Shift some elements of `list` (identified by `index`) to the left end, otherwise preserving order of `list`."""
function shift_elements_right(list,index)
    m = length(list)
    left = [list[i] for i ∈ 1:m if i ∉ index] # Gather specified elements to the left
    right = [list[i] for i ∈ 1:m if i ∈ index]
    return append!(left,right) # Append unspecified elements on the right
end

""" ord = reorder_for_linear_independence(m)
Find the column reordering of matrix m that gathers the left-most linearly independent columns of the matrix to the left """
function reorder_for_linear_independence(m)
    ms, nc = size(m) 
    r = rank(m)
    idx = collect(1:nc)
    #println("Starting rank: ",r)
    it = 1
    while true
        _,R = qr(m[:,idx]) 
        rR = rank(R[:,1:ms])
        ondiag = get_index_for_zeros(R,rR)
        #println("num zero terms: ",length(ondiag),"   num its: ",it, "    Rank of R:",rank(R[:,1:250]))
        if rank(R[:,1:ms])≠r-length(ondiag)
            println("rank of R: ",rank(R[:,1:ms]))
            println(ondiag)
            println("r-length(ondiag): ",r-length(ondiag))
            s = r-length(ondiag)
            println("smallest diag R elements",sort(abs.(diag(R)))[s-3:min(ms,s+3)])
            error("ranks not consistent")
        end
        idx = shift_elements_right(idx,ondiag)
        if rank(R[:,1:ms])==r 
            return idx[1:r]
        end
        it += 1
        if it > nc
            println("Initial rank: ",r)
            println("Final rank: ",rank(R[:,1:ms]))
            println("Length ondiag: ",length(ondiag))
            println("idx: ",idx[1:ms])
            error("Didn't get l.i. set")
        end
    end 
end


""" Return index of list elements with absval ≈ 0 """
function get_index_for_zeros(m, r)
    reps = sort(abs.(diag(m)),rev=:true)[r] # Get the r largest elements by finding the value of the smallest in the set
    idx = findall(x -> x < reps, abs.(diag(m)))
    return idx
end
### end

cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/")
m = readdlm("pimat.li_1-6.ternary")
_,_,nV,_,_ = read_clusters_from_file("clusters.fcc_ternary.1-6_1081")
colors = replace(nV, 1=>:red, 2=>:orange,3=>:green,4=>:magenta,5=>:blue,6=>:red)
pushfirst!(colors,:black)
ms = 200
Ncl = size(m,2)

n = readdlm("pimat.li_1-6.ternary")
p = Int.(readdlm("reorderCl_1081.dat.2"))
colors=colors[p]
n = n[:,p]
Ncl = size(m,2)
Nits = 1000
clProb = zeros(Int,Ncl);
for it ∈ 1:Nits
    println("iteration: ", it)
    #p = randperm(size(n,2))
    #m = n[randperm(1081)[1:ms],p]
    m = n[randperm(1081)[1:ms],:]
    idx = reorder_for_linear_independence(m)  
    #clProb .+= [count(i->i==j,p[idx]) for j ∈ 1:Ncl]  
    clProb .+= [count(i->i==j,idx) for j ∈ 1:Ncl]  
end
plot(clProb./Nits,st=:scatter,msw=0,ms=2,legend=:none,xrange=(1,400),yaxis=:identity,yrange=(0,1),xlabel="Cluster Number",ylabel="Probability",title="Cluster Probability",color=colors)
clProb1 = deepcopy(clProb)
clProb2 = deepcopy(clProb)
plot!(clProb,st=:scatter,msw=0,ms=2,legend=:none,color=:red)

plot(clProb,st=:scatter,msw=0,ms=2,legend=:none,xrange=(1,600),yaxis=:log10,yrange=(0,Nits+1),xlabel="Cluster Number",ylabel="Probability",title="Cluster Probability",
color=colors)

# Find the natural order for the clusters
ρ = clProb/maximum(clProb)
s = sortperm(ρ,rev=:true)
plot(ρ[s],st=:scatter,msw=0,ms=3,color=colors[s],xlabel="Cluster Number",legend=:none,ylabel="Cluster Probability",xrange=(1,400))
writedlm("reorderCl_1081.dat.2",p[s])





# Need some stats on diagonals of R.
# Sort them, calculate all the gaps, take the biggest log10 gap (divides zeros from non-zeros). The zero is the top of that gap. What is the biggest zero?
Nits = 10
ms = 250
for it ∈ 1:Nits
    println("iteration: ", it)
    m = n[randperm(1081)[1:ms],:]
    _,R = qr(m)
    d = abs.(diag(R))
    p = sortperm(d)
    d = d[p]
    gaps = [d[i]-d[i-1] for i ∈ 2:length(d)]  
end
