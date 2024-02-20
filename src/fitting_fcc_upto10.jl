using Plots
using DelimitedFiles
using LinearAlgebra

default(ms=2,msw=0,legend=:none,st=:scatter)
cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/")

#fullEn = readdlm("energiesPerAtom.fcc.upto10")
fullEn = readdlm("relaxed_En.dat")
conc = readdlm("relaxed_conc.dat")

#sm = readdlm("struct_enum_mat.fcc.upto10")
#f1 = Figure(resolution=(1600,900));
#scatter(f1[1,1],conc,fullEn) 
#ax1 = Axis(f1[1,1],xlabel="Ni Concentration")
#p1 = scatter!(ax1,conc,fullEn)
p1 = scatter(conc,fullEn,color=:blue,markersize=5, marker=:diamond; xlabel="Ni Concentration",ylabel="Energy/atom (eV)",title="Convex Hull (Al-Ni)")

m6 = readdlm("pimat.radGen.6li_137")
m = readdlm("pimat.2-6b_li1761")
rank(m6)
rank(m)
En = fullEn[1:137]
testEn = fullEn[138:end]
rank(m6)
rank(m[1:137,:])

s=svd(m6)
plot(s.S,yaxis=:log,title="Singular values of 137 Pi Matrix")
_,R= qr(m6)
plot(abs.(diag(R)),yaxis=:log,title="QR diagonal values of 137 Pi Matrix")
#_,Rt=qr(transpose(m6))
#plot(abs.(diag(Rt)),yaxis=:log)

J = m6\En
col = vcat(:black,:red, [:orange for i ∈ 1:29],[:blue for i in 1:39],[:magenta for i in 1:42],[:green for i in 1:15],[:purple for i in 1:10])
p2 = plot(J,color=col,ms=3,xlabel="Coefficient number",ylabel="Coefficient value")
plot(p1,p2,layout=(2,2),size=(1600,900))
testErr = zeros(136)
err = zeros(136)
for i ∈ 1:136
    J = m6[:,1:i]\En
    err[i] = sum(abs.(m6[:,1:i]*J .- En))/137    
    testErr[i] = sum(abs.(m[138:end,1:i]*J-testEn))/length(testEn)
end
err
plot(1:136,err,title="Error vs terms",xlabel="Number of terms",ylabel="Mean absolute error",color=col,yaxis=:log)
plot!(1:136,testErr,yticks=[1e-4,1e-3,1e-2,1e-1,1,10])
begin
terms = collect(1:137)
keep = []
errList = []
testErr = []
for i ∈ 1:length(terms)-1
    println("i: ",i)
    minerr = 1e300
    best = 0
    for j ∈ terms
        if j ∈ keep continue end
        idx = vcat(keep,j)
        testm = m6[:,idx]
        J = testm\En
        err = sum(abs.(testm*J .- En))/137
        if err < minerr
            minerr = err
            best = j 
        end
    end
    push!(keep,best)
    push!(errList,err)
    J = m6[:,keep]\En 
    push!(testErr,sum(abs.(m[138:end,keep]*J-testEn))/length(testEn))
end
end
plot(1:136,errList,color=col[keep],yaxis=:log,ms=3,xlabel="Number of terms",ylabel="MAE (eV/atom)")
plot!(1:136,testErr,yticks=[1e-4,1e-3,1e-2,1e-1,1,10])
errList

begin
    p6=heatmap(m6,aspect_ratio=1,xlim=(0.5,137.5),ylim=(0.5,137.5),yticks=130:-10:1,xticks=0:10:137,xlabel="Clusters",ylabel="Structure number",legend=:true,yflip=true,dpi=400,tick_direction=:out,framestyle=:box,xmirror=true,grid=false)#,color=:viridis)
    w=.5
    default(st=:line)
    plot!([0,139],[10.5,10.5],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([0,139],[29.5,29.5],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([0,139],[57.5,57.5],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([31.5,31.5],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([70.4,70.4],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([127.4,127.4],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
    plot!([112.4,112.4],[0,139],color=:black,yflip=true,aspect_ratio=1,lw=w,legend=:none)
#    savefig(p6,"pimatrix.png")
end

#m = readdlm("pimat.2-6.reordered")

J10 = m\fullEn
col = vcat(:black,:red, [:orange for i ∈ 1:137],[:blue for i in 138:431],[:magenta for i in 432:947],[:green for i in 948:1388],[:purple for i in 1389:1761])
plot(abs.(J10),yaxis=:log,color=col,ms=3,xlabel="Coefficient number",ylabel="Coefficient value",st=:scatter)
plot(m[1,:],st=:scatter)

# cl, sv, nv, dists, id = read_clusters_from_file("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/clusters.2-6b_li1761")
# tmp = hcat(nv,dists)
# p = sortperm(collect(eachslice(tmp,dims=1)))
# write_clusters(cl[p],sv[p],nv[p],dists[p],id[p],"/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/clusters.2-6.reordered")
# m = readdlm("pimat.2-6.reordered")
# rank(m)

norms = ones(size(m,2))
testErr = ones(size(m,2))
err = ones(size(m,2))
for i ∈ 1:size(m,2)-1
    if i == 137 continue end
    J = m[1:137,1:i]\En
    err[i] = sum(abs.(m[1:137,1:i]*J .- En))/137    
    testErr[i] = sum(abs.(m[138:end,1:i]*J-testEn))/length(testEn)
    norms[i] = norm(J)
end
err
plot(1:size(m,2)-2,err[1:end-1],title="Error vs terms",xlabel="Number of terms",ylabel="Mean absolute error",color=col,yaxis=:log,st=:scatter)
plot!(1:size(m,2)-2,testErr[1:end-1],yticks=[1e-4,1e-3,1e-2,1e-1,1,10],st=:scatter)


plot(norms,st=:scatter,color=col)