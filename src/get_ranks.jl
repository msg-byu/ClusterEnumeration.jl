using DelimitedFiles
using LinearAlgebra
using Random
using Plots
using LaTeXStrings

cd("/Users/glh43/home/fortranCodes/uncle/cluster_tests/")
println("Put all the files in a directory and then change the preceding line.")

""" Parse input files to get number of atoms/cell and energies """
function get_cell_sizes_and_energies()
    # Get number of atoms for each structure
    atoms = readlines("structures.in")
    idx = findall(occursin.("Direct",atoms)).-1
    atoms = atoms[idx]
    atoms = [parse.(Int,i) for i ∈ split.(atoms)].|>sum
 #   grep -B 1 Energy structures.in|grep -v En|grep -v "\-\-" > energies.dat
    # Get energy for each structure, convert to "per atom"
    f=readlines("energies.dat")
    en = parse.(Float64,f)./atoms
    return atoms, en
end

nAt, energies = get_cell_sizes_and_energies()
plot(energies,xlabel="Structure number",ylabel="Energy (eV/atom)",legend=:none)
savefig("energies.pdf")

# Read in the descriptors (clusters)
m = readdlm("enum_PI_matrix.out")


""" Calculate the rank of matrix as more and more clusters are included. Takes a while to run. Don't run it unless you need to """
function get_ranks(m)
    # Calculate rank for increasing number of clusters in PI file
    ranks = zeros(Int,size(m,2))
    Threads.@threads for icol ∈ 1:size(m,2)
        if mod(icol,10)==0 println(icol) end
        ranks[icol] = rank(m[:,1:icol])
    end
    writedlm("ranks.ternary.upto4",ranks)
    return ranks
end
# ranks = get_ranks(m)
# writedlm("ranks.ternary.upto5",ranks)

# Read the ranks in from file (do this instead of running 'get_ranks' unless "structures.in" file changes)
ranks = readdlm("ranks.ternary.upto6")

""" Read cluster file to find vertex order for each cluster """
function parse_clusterfile(filename="clusters.out")
    lines=[]
    vOrder=[0]
    f=readlines(filename)
    for (i,s) in enumerate(f)
        if occursin("Cluster number (",s)
            append!(lines,i+1)
            append!(vOrder,parse(Int,f[i+3]))
        end
    end 
    x = 1:length(vOrder)
    breaks=[findlast(vOrder.<i) for i in 3:6]
    return x, breaks
end
x, breaks = parse_clusterfile()

# Plot the rank of design matrix as a function of adding descriptors (clusters)
begin
default(guidefont=16,tickfont=12,legendfont=12)
plot(x[1:breaks[1]],ranks[1:breaks[1]],color=:orange,lw=4,label="pairs")
plot!(x[breaks[1]+1:breaks[2]],ranks[breaks[1]+1:breaks[2]],color=:green,lw=3,label="triplets")
plot!(x[breaks[2]+1:breaks[3]],ranks[breaks[2]+1:breaks[3]],color=:magenta,lw=3,label="4-bodies")
plot!(x[breaks[3]+1:breaks[4]],ranks[breaks[3]+1:breaks[4]],color=:blue,lw=3,label="5-bodies")
plot!(x[breaks[4]+1:end],ranks[breaks[4]+1:end],color=:red,lw=3,label="6-bodies",ylabel="Matrix Rank",xlabel="Number of basis functions",legend=:bottomright,
ylim = (0,1100), xlim=(0,2200),
title="Ternary FCC Case (n=1-6)")
plot!([0,2220],[1081,1081],ls=:dash,label="",lc=:black)
annotate!(2000,1110,text("Full rank",9))
savefig("~/home/projects/double_descent/paper/figures/rankVSclusters.pdf")
end


ranks = readdlm("ranks.ternary.upto5")
x, breaks = parse_clusterfile("clusters.upto6_2220")
# Plot the rank of design matrix as a function of adding descriptors (clusters)
begin
    default(guidefont=16,tickfont=12,legendfont=10)
    m = readdlm("enum_PI_matrix.upto5")
    ranks = readdlm("ranks.ternary.upto5")
    plot(x[1:breaks[1]],ranks[1:breaks[1]],color=:orange,lw=4,label="pairs")
    plot!(x[breaks[1]+1:breaks[2]],ranks[breaks[1]+1:breaks[2]],color=:green,lw=3,label="triplets")
    plot!(x[breaks[2]+1:breaks[3]],ranks[breaks[2]+1:breaks[3]],color=:magenta,lw=3,label="4-bodies")
    plot!(x[breaks[3]+1:breaks[4]],ranks[breaks[3]+1:breaks[4]],color=:blue,lw=3,label="5-bodies")
    plot!(x[breaks[4]+1:end],ranks[breaks[4]+1:end],color=:red,lw=3,label="6-bodies",ylabel="Matrix Rank",xlabel="Number of basis functions",legend=:bottomright,
    ylim = (0,300), xlim=(0,2200),
    title="Ternary FCC Case (n=1-6)")
    plot!([0,2220],[291,291],ls=:dash,label="",lc=:black)
    annotate!(2000,300,text("Full rank",9))

    m = readdlm("enum_PI_matrix.upto4")
#ranks = get_ranks(m)g
#writedlm("ranks.ternary.upto5",ranks)
    ranks = readdlm("ranks.ternary.upto4")
    plot!(x[1:breaks[1]],ranks[1:breaks[1]],color=:orange,lw=4,label="")
    plot!(x[breaks[1]+1:breaks[2]],ranks[breaks[1]+1:breaks[2]],color=:green,lw=3,label="")
    plot!(x[breaks[2]+1:breaks[3]],ranks[breaks[2]+1:breaks[3]],color=:magenta,lw=3,label="")
    plot!(x[breaks[3]+1:breaks[4]],ranks[breaks[3]+1:breaks[4]],color=:blue,lw=3,label="")
    plot!(x[breaks[4]+1:end],ranks[breaks[4]+1:end],color=:red,lw=3,label="",ylabel="Matrix Rank",xlabel="Number of basis functions",legend=:bottomright,
    #ylim = (0,130), xlim=(0,2200),
    title="Ternary FCC")
    plot!([0,2220],[126,126],ls=:dash,label="",lc=:black)
    annotate!(2000,137,text("Full rank",9))
    annotate!(650,250,text(L"n=1\cdots 5",))
    annotate!(1030,147,text(L"n=1\cdots4",))
        savefig("~/home/projects/double_descent/paper/figuresrankVSclusters.upto5.pdf")
    end
    

println()

# # This alternate approach seem to be a bit faster than using 'rank'
# m2 = m[:,1]
# idx = [1]
# diff = Float64[]
# for icol ∈ 1:size(m,2)-3300
#     mod(icol,10)==0 ? println(icol) : true
#     Δ = norm(m2*(m2\m[:,icol+1])-m[:,icol+1])
#     append!(diff,Δ)
#     if Δ > 1.0e-10
#         m2 = hcat(m2,m[:,icol+1])
#         push!(idx,icol)
#     end
# end

nSt = length(energies)
rndShuffle = shuffle(1:nSt)
enTrain = energies[rndShuffle[1:800]]
enTest = energies[rndShuffle[801:end]]

# RMSE errors for different numbers of descriptors
mTrain = m[rndShuffle[1:800],:]
mTest = m[rndShuffle[801:end],:]

# Plot test and train errors as design matrix is expanded (takes ~15 mins to run)
begin
errTrain = Float64[]
errTest = Float64[]
for n ∈ 1:size(mTrain,2)
    mod(n,50)==0 ? println(n) : true
    c = mTrain[:,1:n]\enTrain # Get coefficients of model
    Δ = sum(abs.(mTrain[:,1:n]*c-enTrain))/size(mTrain,1) # Compute train error
    δ = sum(abs.(mTest[:,1:n]*c-enTest))/size(mTest,1) # test error (MAE)
    append!(errTrain,Δ)
    append!(errTest,δ)
end
end

# Plot train/test errors as a function of adding descriptors (clusters)
begin
plot(x[1:breaks[1]-1],errTrain[1:breaks[1]-1],color=:orange,lw=3,label="pairs")
# plot!(x[breaks[1]:breaks[2]+1],errTrain[breaks[1]:breaks[2]+1],lw=3,color=:green,label="triplets")
# plot!(x[breaks[2]:breaks[3]+1],errTrain[breaks[2]:breaks[3]+1],lw=3,color=:magenta,label="4-bodies")
# plot!(x[breaks[3]:breaks[4]+1],errTrain[breaks[3]:breaks[4]+1],lw=3,color=:blue,label="5-bodies")
# plot!(x[breaks[4]:end],      errTrain[breaks[4]:end],lw=3,      color=:red,label="6-bodies",yaxis=:log)
plot!(errTest,xlabel="Number of descriptors",ylabel="MAE per structure",label="Test set",ylim=(7e-4,0.1))
savefig("errors.pdf")
end

plot(x[1:end],errTrain,yaxis=:log)


""" Compute a "reduced" design matrix that doesn't include any linear dependent terms """
function get_reduced_matrix(m)
    mRed = Array{Float64}(undef,size(m,1),0)
    liList = Array{Int}(undef,0)
    for iR ∈ 1:size(m,2)
        if mod(iR,25)==0 println(iR,' ',size(mRed,2)) end
        testM = hcat(mRed,m[:,iR])
        if rank(testM)==size(testM,2)
            mRed=testM
            append!(liList,iR)
        end
    end
    writedlm("reducedM.dat",mRed)
    return liList, mRed
end
#mRed = get_reduced_matrix(m)
mRed = readdlm("reducedM.dat")

# Plot train/test error for reduced (full rank) design matrix
begin
mRedTrain = mRed[rndShuffle[1:800],:]
mRedTest = mRed[rndShuffle[801:end],:]
errTrain2 = Float64[]
errTest2 = Float64[]
for n ∈ 1:size(mRed,2)
    mod(n,50)==0 ? println(n) : true
    c = mRedTrain[:,1:n]\enTrain
    Δ = sum(abs.(mRedTrain[:,1:n]*c-enTrain))/size(mRedTrain,1)
    δ = sum(abs.(mRedTest[:,1:n]*c-enTest))/size(mRedTest,1)
    append!(errTrain2,Δ)
    append!(errTest2,δ)
end
end 

writedlm("mRed_errTrain.dat",errTrain2)
writedlm("mRed_errTest.dat",errTest2)

begin
plot(errTrain2,label="Train")
plot!(errTest,label="Test (orig)",xlim=(1,841))
plot!(errTest2,ylim=(7e-4,.1),yaxis=:log,label="Test",xlabel="Number of descriptors",ylabel="MAE prediction (eV/atom)")
savefig("errors_reduced_matrix.pdf")
end

# Examine condition numbers. Keep only clusters that don't increase the condition number by more than a factor of 5
#conds = [cond(mRed[:,1:i]) for i ∈ 1:size(mRed,2)-641]
conds = []
mCond = Array{Float64}(undef,size(m,1),0)
r = 1
for i ∈ 1:size(mRed,2)
    mod(i,50)==0 ? println(i) : true
    s = cond(mRed[:,1:i])
    if s/r < 1.5
        mCond = hcat(mCond,mRed[:,1:i])
        append!(conds,s)
        r = s
    else
        println("skipped: ",s/r," ",i)
    end
end
plot(conds,yaxis=:log,xlabel="Number of terms",ylabel="Condition number",legend=:none)
plot(conds[2:end]./conds[1:end-1],xlabel="Number of terms",ylabel="% Increase",legend=:none)
# The condition number doesn't seem to be a good metric for eliminating clusters

plot(errTrain[2:end]./errTrain[1:end-1],yaxis=:log,ylim=(1e-1,1.2),ylabel="Train Error Reduction",legend=:none)
savefig("UsingClustersThatImproveTrainErr.pdf")
plot(errTest[2:end]./errTest[1:end-1],yaxis=:log,ylim=(-.2,1.2),xlim=(1,100))
plot(errTest.-errTrain,ylim=(-.01,.12),ylabel="Test error increase over train error",xlabel="Number of terms",legend=:none)

# If we eliminate clusters that don't reduce train error much, will that reduce test error?
impv = []
cIdx = [1]
mImpTr = mRedTrain[:,1]
mImpTe = mRedTest[:,1]
c = mImpTr\enTrain
r = sum(abs.(mImpTr[:,1]*c-enTrain))/size(mImpTr,1)
errTrImp = [r]
errTeImp = [sum(abs.(mImpTe[:,1]*c-enTest))/size(mImpTe,1)]
for i ∈ 2:size(mRed,2)
    mod(i,50)==0 ? println(i) : true
    test = hcat(mImpTr,mRedTrain[:,i])
    c = test\enTrain
    err = sum(abs.(test*c-enTrain))/size(mImpTr,1)
    imp = r/err 
    if imp > 1.005
        mImpTr = test
        append!(impv,imp)
        r = err
        append!(cIdx,i)
        append!(errTrImp,err)
        mImpTe = hcat(mImpTe,mTest[:,i])
        append!(errTeImp,sum(abs.(mImpTe*c-enTest))/size(mImpTe,1))
    else
        println("skipped: ",i,"  ",imp)
    end
end

end

plot(errTeImp,label="Test")
plot!(errTrImp,yaxis=:log,xlabel="Number of terms",ylabel="Train error",label="Train")
# Doesn't help

