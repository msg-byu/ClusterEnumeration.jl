module clusters
using Plots
using Combinatorics
using Spacey
using LinearAlgebra
using DelimitedFiles
using Printf
using MinkowskiReduction
#using SMatrix

export get_dsites_cell_sizes, gen_points_in_sphere, read_lattice_vectors, make_figure_candidates, make_eqvPoints, reduce_figList, write_clusters, get_nonzero_index, MADfromCOM, read_clusters_from_file, get_leftmost_indep_columns, rd_list, get_MAD, isMatrixEqvBySymmetry

""" Use QR decomposition iteratively to pull out left-most independent columms """
function get_leftmost_indep_columns(m)
    r = rank(m)
    nr,nc = size(m)
    idx = collect(1:nc)
    it = 1
    mask = trues(nc)
    while true
        println("it: ",it)
        println("len idx: ",length(idx))
        if it > 80 break end
        _,R = qr(m[:,idx])
        ondiag = get_nonzero_index(R)
        r1 = rank(m[:,idx[1:length(ondiag)][ondiag]])
        println("rank: ",r1)
        if r1==r 
            idx = idx[1:length(ondiag)][ondiag]
            break 
        end
        pad = max(length(idx)-nr,0)
        mask[idx] = vcat(ondiag,trues(pad))
        idx = collect(1:nc)[mask]
        it += 1
    end
    rank(m[:,mask])
    return idx
end


""" Check and see if the all the entries in each column of m are equal """
function isCommonShift(m)
    return all(y->y≈m[:,1],eachcol(m)) 
end

""" Remove translationally equivalent figures from a list """
function remove_transEqv_figures!(figs)
    uqFigs = []
    for jFig ∈ figs
        uq = true
        for iFig ∈ uqFigs
            if length(iFig) == length(jFig)
                r = iFig - jFig # Turn r into a matrix
                # All all elements of each row equal?
                if isCommonShift(r)
                    println("Common shift") 
                    display(iFig)
                    display(jFig)                
                    uq = false
                    break
                end
            end
        end
        if uq
            push!(uqFigs,jFig)
        end
    end
    figs = uqFigs
    return figs
end


""" Get HNF index and poingroup of each structure in struct_enum file  """
function cell_index_and_pointgroup(filename = "struct_enum.out")
    f = readlines(filename)[17:end]
    hnfIdx = [parse(Int,split(i)[2]) for i in f]
    pgSize = [parse(Int,split(i)[8]) for i in f]
    return hnfIdx, pgSize
end

""" Get the atom positions of all structures from 'structures.in'.
    Return cell sizes as well, and longest lattice vector """
function get_dsites_cell_sizes(filename="structures.in")
    f = readlines(filename)
    s = findall(x->occursin("Direct",x),f).+1
    e = findall(x->occursin("Energy",x),f).-2
    idx = hcat(s,e)
    cell_sizes = [j[2]-j[1]+1 for j in eachrow(idx)]
    # Collect all atom positions for each structure
    dsites = []; maxLen = Float64[]
    for j in 1:size(idx,1)
        nextSites = [parse.(Float64,split(i)[1:3]) for i in f[idx[j,1]:idx[j,2]]]
        # Get lattice vectors and convert from direct to Cartesian coordinates
        S = f[idx[j,1].-collect(5:-1:3)]
        S = hcat([parse.(Float64,split(s)) for s ∈ S]...)
        longest = maximum(norm.(eachcol(S)))
        nextSites = S*hcat(nextSites...)
        nextSites = [nextSites[:,i] for i ∈ 1:size(nextSites,2)]
        push!(dsites,nextSites)
        push!(maxLen,longest)
    end
    return dsites, cell_sizes, maxLen
end

""" Read in the lattice vectors from lat.in """
function read_lattice_vectors()
    g = readlines("lat.in")
    A = hcat([parse.(Float64,split(i)) for i ∈ g[6:8]]...)
    return A
end

""" Expand a list of positions/structure (from 'structures.in') to cluster combinations """
function generateFiguresFromSites(dsites)
    # For each structure, generate all candidate clusters of all vertex orders
    candClust = [collect(combinations(jSt)) for jSt ∈ dsites]
    # Undo one layer of nesting (structure-by-structure nesting not needed now)
    # candClust will now be a vector of clusters (a cluster is a vector of vectors)
    candClust = vcat(candClust...)
    # Toss out single-vertex cluster candidates
    candClust = filter(x->length(x)> 1,candClust) 
    # Convert from vectors of vectors to matrices 
    candClust = [hcat(i...) for i ∈ candClust]
    # Sort the columns, put vertices in a canonical order (handy for comparisons)
    candClust = [sortslices(i,dims=2) for i ∈ candClust]
    # Shift all clusters so that vertex 1 is the origin
    candClust = [i.-i[:,1] for i ∈ candClust]
    # Toss out exact duplicates (tossing out other kinds of duplicates will happen later)
    candClust = unique(candClust)
    # Sort/group clusters by number of vertices (not necessary but aesthetically pleasing) 
    candClust = candClust[sortperm(size.(candClust,2))]
    return candClust
end

""" Attach s-vectors to figures """
function attachSvectors(clusters)
    sV = [permutations(vcat([i*ones(Int,j) for i in 1:2]...),j)|>collect|>unique for j in 2:6]
    svec = []
    scl = []
    for iCl ∈ clusters
        l = size(iCl,2)
        for iSvec ∈ 1:2^l # Loop over all the s-vectors for this vertex order 
            push!(scl,iCl) # Just another copy of the geometric figure
            push!(svec,sV[l-1][iSvec])
        end
    end
    return scl, svec
end

""" Generate list of symmetry-equivalent points from list of points """
function make_eqvPoints(pointList,rots)
    pts = []
    for i ∈ pointList
        push!(pts,unique([irot*i for irot ∈ rots])...)
    end
    return pts
end

function MADfromCOM(fig)
    nV = size(fig,2)
    com = sum(fig,dims=2)./nV
    diffs = [2*norm(i.-com) for i ∈ eachcol(fig)]
    return sum(diffs)./nV
end

""" make_figure_candidates(point_list, n)

Make figures (with 'n' vertices, n > 1) from list of points. Figures always include the origin. This routine does not remove translationally or rotationally equivalent figures.
          
"""
function make_figure_candidates(pts, n)
    figs = collect(combinations(pts,n-1))
    [pushfirst!(ifig,[0.,0.,0.]) for ifig ∈ figs]
    figs = [hcat(ifig...) for ifig ∈ figs] # Convert vecs of vecs to matrices 
    figs = [sortslices(i,dims=2) for i ∈ figs] # Sort the vertices
    return figs[sortperm(MADfromCOM.(figs))]
end

""" reduce_figList(figures_list,symmetry_operations)

Remove translationally and symmetrically equivalent figures.
"""
function reduce_figList(figList,rots)
    # s = size(figList[1])
    # println("s: ",s)
    # l = length(figList)
    # println("l: ",l)
    radii = MADfromCOM.(figList)
    p = sortperm(radii)
    figList = figList[p]
    radii = radii[p]
    uqFigs = Vector{Matrix{Float64}}()
    for (i,iFig) ∈ enumerate(figList)
        uqFlag = true
        if mod(i,10_000)==0 println(i," ",length(uqFigs)) end
        for jFig ∈ uqFigs
            #print(".")m
            rj = MADfromCOM(jFig)
            if radii[i]≉rj && radii[i] > rj 
                continue
            elseif radii[i]≉ rj && radii[i] < rj
                break
            end
            #print("x")
            for iRot ∈ rots
                t = sortslices(replace!(iRot*iFig,-0.0=>0.0),dims=2;lt=<)
                r = t - jFig # Compute translation difference between the vertices
                if t ≈ jFig # then iFig is a rotation duplicate of jFig
                    uqFlag = false
                    break
                elseif all(y->y≈r[:,1],eachcol(r)) # Then iFig is a rotated translation duplicate
                    uqFlag = false
                    break
                end
            end
            if uqFlag==false break end
        end
        if uqFlag
            #print("uq")
            push!(uqFigs,iFig)
        end
        #println()
    end
    return uqFigs
end

""" isMatrixEqvBySymmetry(iFig,jFig,rots) """
function isMatrixEqvBySymmetry(iFig,jFig,rots)
    uqFlag = false
    for iRot ∈ rots
        t = sortslices(replace!(iRot*iFig,-0.0=>0.0),dims=2) # The 'lt=<' avoids sorting problems with negative zeros
        r = t - jFig # Compute translation difference between the vertices
        if t ≈ jFig # then iFig is a rotation duplicate of jFig
            uqFlag = true
            break
        elseif all(y->y≈r[:,1],eachcol(r)) # Then iFig is a rotated translation duplicate
            uqFlag = true
            break
        end
    end
    return uqFlag
end

function rd_list2(figList,rots)
    s = size(figList[1])
    println("s: ",s)
    l = length(figList)
    println("l: ",l)
    radii = get_MAD.(figList)
    p = sortperm(radii)
    figList = figList[p]
    radii = radii[p]
    mask = falses(length(figList))
    for (i,iFig) ∈ enumerate(figList)
        uqFlag = true
        if mod(i,10_000)==0 println(i," ",count(mask)) end
        for jFig ∈ figList[mask]
            #print(".")
            rj = get_MAD(jFig)
            if radii[i]≉rj && radii[i] > rj 
                continue
            elseif radii[i]≉ rj && radii[i] < rj
                break
            end
            uqFlag = isMatrixEqvBySymmetry(iFig,jFig,rots)
            if uqFlag==false break end
        end
        if uqFlag
            mask[i] = true
        end
    end
    return figList[mask]
end

""" Read in clusters info from the 'clusters.out' file 
    clusters, svecs
    , Verts, avgDist = read_clusters_from_file("[file]")
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

""" gen_points_in_supercell(A,n)
Generate the shortest 'n' the points on the lattice defined by A 
This doesn't work robustly yet. Will provide good results only in ideal cases.
"""
function gen_points_in_supercell(A,npts)
    u,v,w = eachcol(A)
    _, rots = pointGroup_robust(u,v,w)
    s = ceil(Int,∛npts)
    l,m,n = [s,s,s]
    pts =[A*[i, j, k] for i ∈ -l:l for j ∈ -m:m for k ∈ -n:n]
    pts = pts[sortperm([norm(i) for i ∈ pts])]
    uqPts = []
    for iPt ∈ pts 
        l = length(uqPts)
        if mod(l,10000)==0 println("pts so far: ",l) end
        uq = true
        for jPt ∈ uqPts
            if norm(iPt) > norm(uqPts[end]) break end # should speed up things
            if norm(iPt) > 45 break end
            if norm(iPt)≉norm(jPt) continue end # This speeds things up immensely
            for iRot ∈ rots
                if iRot*iPt≈jPt
                    uq = false
                    break
                end
            end
            if uq==false
                break
            end
        end
        if uq
            push!(uqPts,iPt)
        end
    end
    return uqPts[2:end] # Skip point at origin
end

""" Write out clusters & s-vectors to clusters.out-type file """
function write_clusters(clusters, svectors, numVerts, avgDist, idLabel,filename="clusters.out")
    str = 
"""#--------------------------------------------------------------------------------
# Cluster number: 0
# This is the empty cluster = constant term for the cluster expansion.
#
# - The empty cluster is not explicitly listed in this file. 
#   (This comment here only reminds the user of the fact that
#   the empty cluster does exist in UNCLE.)
# - The empty cluster does not contain any vertices.
# - The empty cluster is not read in from this file.
# - The empty cluster is always used in the fitting procedure.
# - The empty cluster counts as a separate cluster for the fitting procedure.
#   (E.g., a simple fit with 1 cluster is a fit with the empty cluster only.)
"""
    nCl = 0
    for iFig ∈ 1:length(clusters)
        println("iFig ",iFig)
        nV = size(clusters[iFig],2)
        if nV!=numVerts[iFig]
            println("numVerts: ",numVerts[iFig])
            println("nV:     : ",nV)
            error("Num Verts doesn't agree")
        end
        thiskind = count(size.(clusters[1:iFig],2).==nV)
        fig = clusters[iFig] # Grab the vertex coordinates
        #com = sum(fig,dims=2)/nV # Compute center of mass
        #radius = √(sum((fig.-com).^2)/nV)
        radius = avgDist[iFig]
        str *= "#--------------------------------------------------------------------------------\n"
        str *= @sprintf "# Cluster number (total/of this kind): %5d / %6d\n" iFig thiskind
        str *= @sprintf "%6d\n" idLabel[iFig]
        str *= "# Number of vertices:\n   $nV\n# Avg. Distance:\n"        
        str *= @sprintf "%10.6f\n" radius
        str *= "# Vertices: (x,y,z) | d-vector label | s-vector\n"
        for iV ∈ 1:nV
            a,b,c = fig[:,iV]
            sVlabel = svectors[iFig][iV]
            str *= @sprintf "  %12.6f" a
            str *= @sprintf "    %12.6f" b
            str *= @sprintf "    %12.6f" c
            str *= @sprintf "         1       %d\n" sVlabel
        end
        str *= "#\n"
    end
    println(pwd())
    open(filename,"w") do io 
        println(io,str)
    end
end

function get_nonzero_index(m,reps=1e-13)
    mask = abs.(diag(m)).>reps
    return mask
end


end
# d, cellSize, len = get_dsites_cell_sizes()
# plot(cellSize,len,st=:scatter,title="Maximum cell vector length vs # of atoms",xlabel="Number of atoms",ylabel="Max cell edge length (a0%)",legend=:none,msw=0)
# plot(cellSize,title="Cell sizes",legend=:none,xlabel="Structure number",ylabel="Volume factor")

# #m = readdlm("enum_PI_matrix.out.2n3b")

# m = m[:,3:end]
# uq = [length(unique(round.(i,digits=14))) for i in eachrow(m)]
# plot(range(1,241),uq[1:241],title="Count of unique π's vs. structure",xlabel="Structure number",ylabel="Count of unique correlation values",label="Size < 8")
# plot!(range(242,631),uq[242:631],color=:red,label="Size 8")
# plot!(range(632,1135),uq[632:1135],color=:purple,label="Size 9")
# plot!(range(1136,length(uq)),uq[1136:end],color=:magenta,label="Size 10")

# id,pg = cell_index_and_pointgroup()
# plot(id)
# plot!(range(1136,length(id)),id[1136:end],color=:magenta)
# plot!(range(242,631),id[242:631],color=:red,label="Size 8")

# p = sortperm(id)
# id = id[p]
# plot(id,uq[p],st=:scatter,msw=0)

# ratio=[maximum(uq[id.==i]) for i ∈ 1:87]./[unique(cellSize[id.==i])[1] for i in 1:87]
# plot(ratio,legend=:none,
# title="ratio of unique π's to cell size",xlabel="Cell size-shape index",ylabel="Ratio",color=:red,msw=0,st=:scatter)

# # The last unique π in row, can come after many repeats. E.g.,
# plot([findfirst(x->x≈i,m[900,:]) for i ∈ unique(m[900,:])])
# plot(m[900,:])

# # It seems that the last unique value in every case is 1. This may be a bug in the correlation calculations, or it might be a consquence of...?
# [unique(i)[end] for i ∈ eachrow(m)]

# # How far out is last unique π over all structures?
# plot([findfirst(x->x≈1,i) for i ∈ eachrow(m)])
# argmax([findfirst(x->x≈1,i) for i ∈ eachrow(m)])
# findfirst(x->x≈1,m[1136,:])
# rank(m[:,1:177])
# # Going out as far as the last unique π, in any row, does not yeld the maximum ranks...
# rank(m[:,1:end])

# # There doesn't seem to be a problem with small singular values in the rank computation. Big break at 136
# s=svd(m)
# plot(s.S,yaxis=:log)

# # The 369rd column is the last lin. ind. one
# rank(m[:,1:368])
# plot(m[:,369])

# # Plot the impact of including that last pair on train and test error

# #pending

# #
# clen = readdlm("clusters.lengths")
# plot(clen,legend=:none,xlabel="Cluster number",ylabel="Pair length")

# # Structures 8,9,10 can't be distinguished from previous structures with ANY pair cluster. (must need a three body)
# # Ditto for 16-29
# [rank(m[1:i,:]) for i ∈ 1:35]|>plot

# # Really large pair cluster pool
# #m = readdlm("enum_PI_matrix.out")
# mfull = readdlm("PImatrix.2n3Full")
# rank(mfull)
# # pivoting picks up the 3 columns beyond the end of the diagonal, but pivoting does not pick the leftmost vector when it has a choice.
# _, R, p = qr(mfull,ColumnNorm());
# rank(mfull[:,p[1:415]])
# _, R = qr(mfull);
# rank(R) #412 l.i. columns before the end of the diagonal
# cl,sv,nv,dist,id = read_clusters_from_file("clusters.out.2n3")
# write_clusters(cl,sv,nv,dist,id,"clusters.2n3Full_write")
# # cp clusters to clusters.out, run mode 44,
# rank(readdlm("enum_PI_matrix.out")[:,1:2346])
# # So reading in the clusters, writing them out, regenerating the PI file works.
# # What if we delete the zero columns?
# mask = get_nonzero_index(R,1e-12)
# mask = vcat(mask[2:end],falses(length(cl)-size(R,2)+2))
# size(mfull)
# count(mask)
# write_clusters(cl[mask[2:end]],sv[mask[2:end]],"clusters.nonzeros")
# # cp clusters to clusters.out, run mode 44,
# rank(readdlm("enum_PI_matrix.out"))
# plot(mask[1:200])



# ##############################################################
# # Really large pair + 3body cluster pool
# cd("/Users/glh43/home/fortranCodes/ClusterGeneration.jl/pairs")
# m = readdlm("PImatrix.2n3Full")
# rank(m)
# #_, R, p = qr(m,ColumnNorm());
# #heatmap(log.(abs.(R[end:-1:1,:])),size=(1600,1000),yticks=:none)
# _, R = qr(m);
# mask=vcat(abs.(diag(R)).>1e-13,falses(3289-2346));
# count(mask) #415 l.i. columns
# plot(collect(1:length(mask))[mask])
# rank(m[:,mask],eps==1e-13) # without the eps option, we come up 3 short (412 instead of 415)

# cl,sv,nv,dist,id = read_clusters_from_file("clusters.2n3Full")
# licl = cl[mask[2:end]]
# lisv = sv[mask[2:end]]
# linv = nv[mask[2:end]]
# lidist=dist[mask[2:end]]
# liid=id[mask[2:end]]
# write_clusters(licl,lisv,linv,lidist,id,"clusters.2n3.li")

# # Now write out the new PI file
# mli = readdlm("PImatrix.2n3.li.only")
# rank(mli,eps==1e-13) # 415, as expected
# # YAY!
# m4 = readdlm("PImatrix.2thru4Full")
# rank(m4) #899 when 2008 4b are used, for n<=4 cell sizes, full rank (so we got them all)
# _, R4 = qr(m4)
# plot(abs.(diag(R4)),yaxis=:log)
# tmask=vcat(abs.(diag(R4)).>1e-13,trues(size(R4,2)-size(R4,1)));
# idx = collect(1:size(R4,2))
# rank(m4[:,tmask])
# idx[tmask]
# _, R4p = qr(m4[:,tmask])
# tmask = idx[tmask][get_nonzero_index(R4p)]
# mask = falses(size(R4,2))
# for i in tmask
#     mask[i]=1
# end
# plot(mask,st=:scatter,legend=:none,msw=0,ms=2,color=:purple) 

# cl,sv,nv,dist,id = read_clusters_from_file("clusters.2thru4Full")
# licl = cl[mask[2:end]]
# lisv = sv[mask[2:end]]
# linv = nv[mask[2:end]]
# lidist=dist[mask[2:end]]
# liid=id[mask[2:end]]
# write_clusters(licl,lisv,linv,lidist,id,"clusters.2thru4.li")

# # 5 bodies
# m5 = readdlm("PImatrix.2thru5Full")
# rank(m5)#1302 
# rank(m5[1:58,:]) #58 So full rank for 5bodies and 5cells
# rank(m5[1:137,:]) #126
# _, R5 = qr(m5);
# plot(abs.(diag(R5)),yaxis=:log,st=:scatter,ms=2,msw=0,legend=:none,color=:green,ylabel="R diag component size",xlabel="Cluster number")

# tmask=vcat(abs.(diag(R5)).>1e-13,trues(size(R5,2)-size(R5,1)));
# idx = collect(1:size(R5,2))
# rank(m5[:,tmask])
# idx[tmask]
# _, R5p = qr(m5[:,tmask])
# tmask = idx[tmask][get_nonzero_index(R5p)]
# mask = falses(size(R5,2))
# for i in tmask
#     mask[i]=1
# end
# plot(mask,st=:scatter,legend=:none,msw=0,ms=2,color=:purple) 

# cl,sv,nv,dist,id = read_clusters_from_file("clusters.2thru5Full")
# licl = cl[mask[2:end]]
# lisv = sv[mask[2:end]]
# linv = nv[mask[2:end]]
# lidist=dist[mask[2:end]]
# liid=id[mask[2:end]]
# write_clusters(licl,lisv,linv,lidist,id,"clusters.2thru5.li")
# m=readdlm("PImatrix.2thru5.li")
# rank(m) 

# # 6 bodies
# m6 = readdlm("PImatrix.2thru5Full")
# rank(m6)#1621 
# rank(m6[1:137,:]) #134 So *not* full rank for 6bodies and 6cells, incomplete list of clusters?
# rank(m6[1:241,:]) #220
# temp = readdlm("PImatrix.6only")
# m = hcat(m6,temp)
# rank(m[1:137,:])

# #Break cluster file for 6b into two big sections and try combining after generating PI files to avoid line length restriction


# _, R6 = qr(m6);
# plot(abs.(diag(R6)),yaxis=:log)
# mask=vcat(abs.(diag(R5)).>1e-13,trues(size(R5,2)-size(R5,1)));
# cl,sv,nv,dist,id = read_clusters_from_file("clusters.2thru5Full")
# licl = cl[mask[2:end]]
# lisv = sv[mask[2:end]]
# linv = nv[mask[2:end]]
# lidist=dist[mask[2:end]]
# liid=id[mask[2:end]]
# write_clusters(licl,lisv,linv,lidist,id,"clusters.2thru5.li")
# rank


# ## Example using only 1-4 size
# h = mfull[1:29,:]
# _, rh = qr(h)
# heatmap(rh[:,450:-1:1])
# rank(rh[:,1:15])

# h = readdlm("PImatrix.2n3small")[1:29,:]
# rank(h)
# _, rh = qr(h)
# heatmap(rh[27:-1:1,:])
# plot(abs.(diag(rh)),yaxis=:log)
# count(abs.(diag(rh)).>1e-5)
# mask = abs.(diag(rh)).>1e-5
# cl,sv = read_clusters_from_file("clusters.small2n3_27")

# write_clusters(cl[mask[2:end]],sv[mask[2:end]],"clusters.smalltest")
# rank(readdlm("enum_PI_matrix.out")[1:29,:])
# ####################################
# # Generate 6-body clusters from first three neighbors

# A = read_lattice_vectors()
# u,v,w = eachcol(A)
# _, rops =  pointGroup_robust(u,v,w)


# pool = [zeros(3)]
# a1 = A[:,1]
# pool = append!(pool,unique([i*a1 for i in rops]))
# a2 = [1.;0;0]
# push!(pool,unique([i*a2 for i in rops])...)
# a3 = [1.0;.5;.5]
# push!(pool,unique([i*a3 for i in rops])...)




# #####################################################
# mask = trues(size(candClust,1))
# for (i,iCl) ∈ enumerate(candClust)
#     for (j,jCl) ∈ enumerate(candClust[1:i-1])
#         if length(iCl) == length(jCl) && iCl≈jCl
#             mask[i] = false
#         end
#     end
# end
# # 638 is reduced to 408 when 'isapprox' is used to compare
# candClust = candClust[mask]
# plot(length.(candClust),title="Unique Candidate figures",legend=:none)


# # Finally we need to eliminate clusters that are equivalent under rotation (and rotation+lattice shift)
# mask = trues(size(uqClust,1))
# for (i,iCl) ∈ enumerate(uqClust)
#     for (j,jCl) ∈ enumerate(uqClust[1:i-1])
#         if length(iCl) != length(jCl) 
#             continue
#         end
#         for iRot ∈ ops
#             t = sort(transpose(transpose(iCl)*iRot),dims=2) # This is a faulty concept. Can't rotate a vector with the integer rots
#             r = t - jCl # Compute translation difference between the vertices
#             if t ≈ jCl # then iCl is a rotation duplicate
#                 mask[i] = false
#                 break
#             #All all elements of each row equal?
#             elseif isCommonShift(r) # Then iCl is a rotated translation duplicate
#                 mask[i] = false
#                 println("shift dup")
#                 break
#             end
#         end
#         if !mask[i] break end
#     end
# end
# rotuqClust = uqClust[mask]
# # 277 of these found July 11 8:50 pm
# # 267 at 9:40 (added sorting of vertices to 't')
# # 108 at 9:30 July 12 (shifted all clusters to have vertex 1 at origin)...but now they are unsorted...
# # 314 at 2 am Juli 13

# nFigPerOrder = [sum(length.(rotuqClust).==i*3) for i in 2:6]
# nClustPerOrder = [2^(i+1)*nFigPerOrder[i] for i ∈ 1:length(nFigPerOrder)]
# println("Number of figures and clusters at each order (2-6): \n",nFigPerOrder,"\n",nClustPerOrder)
# println("Total number of clusters (including s-vectors): ",sum(nClustPerOrder))



# ## Now that we have the s-vectors, we need to eliminate s-clusters that are symmetrically equivalent
# mask = trues(size(sCl,1))
# for (i,iCl) ∈ enumerate(sCl[1:end])
#     println(i)
#     for (j,jCl) ∈ enumerate(sCl[1:i-1])
#         if length(iCl) != length(jCl) 
#             continue
#         end
#         for iRot ∈ ops
#             t = transpose(transpose(iCl)*iRot) # Rotate the i-th sCl
#             p = sortperm(collect(eachslice(t,dims=2))) # how t is permuted by the sort
#             t = t[:,p] # sort t
#             svec = sVlist[i][p] # Permute the labels on the vertices
#             if t ≈ jCl && svec==sVlist[j]
#                 mask[i] = false
#                 break
#             end
#             r = t - jCl
#             if isCommonShift(r) && svec==sVlist[j]
#                 mask[i] = false
#                 println("breaking")
#                 break
#             end
#         end
#         if !mask[i] break end
#     end
# end

# end
