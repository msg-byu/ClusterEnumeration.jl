# Used this file to find cluster sets for 8-atom fcc configs. Also a bunch of cluster-related function but the intent is that all of these are superceded by the functions in ClusterEnumeration.jl
#
# Finally, used this code to get the full-s-vector clusters for n=1-6 ternary case
module clusters
using Plots
using Combinatorics
using Spacey
using LinearAlgebra
using DelimitedFiles
using Printf
using MinkowskiReduction
using Random
#using SMatrix

export get_dsites_cell_sizes, gen_points_in_supercell, read_lattice_vectors, make_figure_candidates, make_eqvPoints, reduce_figList, write_clusters, get_nonzero_index, MADfromCOM, read_clusters_from_file, get_leftmost_indep_columns, rd_list, get_MAD, isMatrixEqvBySymmetry

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
        #if it > 10*ceil(nc/nr) break end
        if length(idx) ≤ nr break end
        #if it > 25 break end 
        _,R = qr(m[:,idx])
        ondiag = get_nonzero_index(R;reps=1e-8)
        r1 = rank(m[:,idx[1:length(ondiag)][ondiag]])
        println("rank: ",r1)
        if r1==r 
            idx = idx[1:length(ondiag)][ondiag]
            break 
        end
        pad = max(length(idx)-nr,0)
        mask[idx] = vcat(ondiag,trues(pad))
        idx = collect(1:nc)[mask]
        println("last: ",idx[end])
        println("rank idx: ",rank(m[:,idx]),"\n")
        it += 1
        if length(idx) < nr break end
    end
    #rank(m[:,mask])
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

""" Get the atom positions of all structures from 'to '.
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

# This function is needed any more since the "clusters from internal coords" idea didn't pan out.
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

""" attachSvectors(clList)

Attach s-vectors to figures 

(Actually, this only works for ternaries as currently written)
"""
function attachSvectors(clusters)
    sV = [permutations(vcat([i*ones(Int,j) for i in 1:2]...),j)|>collect|>unique for j in 1:6]
    svec = []
    scl = []
    for iCl ∈ clusters
        println(iCl)
        l = size(iCl,2)
        for iSvec ∈ 1:2^l # Loop over all the s-vectors for this vertex order 
            push!(scl,iCl) # Just another copy of the geometric figure
            push!(svec,sV[l][iSvec])
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

""" get_nonzero_index(m,reps=1e-13) """

function get_nonzero_index(m; reps=1e-13)
    mask = abs.(diag(m)).>reps
    return mask
end


end # End of module "clusters"

# d, cellSize, len = get_dsites_cell_sizes()
# plot(cellSize,len,st=:scatter,title="Maximum cell vector length vs # of atoms",xlabel="Number of atoms",ylabel="Max cell edge length (a0%)",legend=:none,msw=0)
# plot(cellSize,title="Cell sizes",legend=:none,xlabel="Structure number",ylabel="Volume factor")

cd("/Users/glh43/home/fortranCodes/uncle/tests/david/bulk_CE")
m = readdlm("enum_PI_matrix_full.out")

rank(m)

size(m)
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


### Working in a single cell shape Start with doubled fcc (8 atoms, 16 cells)
# # Really large pair cluster pool
using .clusters
m = readdlm("enum_PI_matrix.2b_300")
rank(m) # 4 - empty, on-site, and two pairs
_,R2 = qr(m)
mask = get_leftmost_indep_columns(m)
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1
cl,sv,nv,dist,id = read_clusters_from_file("clusters.2b_300")
licl = cl[mask]
lisv = sv[mask]
linv = nv[mask]
lidist=dist[mask]
liid=id[mask]
#Includes the empty and on-site in number of pairs
write_clusters(licl,lisv,linv,lidist,id,"clusters.2b_li4")

A = read_lattice_vectors()
_, rots = pointGroup_robust(eachcol(A)...) 

rots = [A*i for i in pointGroup_fast(eachcol(A)...)] # This wasn't returning the Cartesion operators...
rots = [[1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 -1.0 0.0], [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], [-1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 -1.0 0.0], [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 -1.0], [0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0], [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0 1.0 0.0; 0.0 0.0 -1.0; -1.0 0.0 0.0], [0.0 0.0 1.0; 0.0 -1.0 0.0; -1.0 0.0 0.0], [0.0 1.0 0.0; 0.0 0.0 -1.0; 1.0 0.0 0.0], [0.0 0.0 1.0; 0.0 -1.0 0.0; 1.0 0.0 0.0], [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 -1.0], [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0], [0.0 0.0 -1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0 0.0 -1.0; 1.0 0.0 0.0; 0.0 -1.0 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0], [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0], [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0], [-1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0], [0.0 -1.0 0.0; 0.0 0.0 -1.0; -1.0 0.0 0.0], [0.0 0.0 -1.0; 0.0 -1.0 0.0; -1.0 0.0 0.0], [0.0 -1.0 0.0; 0.0 0.0 -1.0; 1.0 0.0 0.0], [0.0 0.0 -1.0; 0.0 -1.0 0.0; 1.0 0.0 0.0], [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], [0.0 1.0 0.0; 0.0 0.0 1.0; -1.0 0.0 0.0], [0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0], [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0], [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0], [-1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0], [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0], [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0], [0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0], [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.0 0.0 -1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], [0.0 -1.0 0.0; 0.0 0.0 1.0; -1.0 0.0 0.0], [0.0 0.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0], [0.0 -1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0], [0.0 0.0 -1.0; -1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0], [0.0 0.0 -1.0; 1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [-1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0], [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]]
# begin 
 # modified this for 7-body cluster enumeration
@time Rpts = gen_points_in_supercell(A,30)
Rpts = Rpts[1] # gen_points returned lengths too
@time pts = make_eqvPoints(Rpts[1:10000],rots)
figs7bt = make_figure_candidates(pts[1:50],7)
figs7b = symReduceByLengthFigList(figs7bt,rots)
n7b = length(figs7b)
cl,sv,nv,dist,id = read_clusters_from_file("clusters.2-6b_li1761")
clust7b = vcat(cl,figs7b)
sv7 = vcat(sv,[[1,1,1,1,1,1,1] for i ∈ 1:n7b])
id7 = vcat(id,1:n7b)
dists7 = vcat(dist,MADfromCOM.(figs7b))
nV7 = vcat(nv,ones(Int,n7b)*7)
write_clusters(clust7b, sv7, nV7, dists7, id7, "clusters.2-7b_"*string(n7b))
#>  cp clusters.2-7b_234000 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-7b_"*string(n7b))
rank(m) # 6 (so added two 3-body clusters)
_,R7 = qr(m)
mask = get_leftmost_indep_columns(m)
rank(m) #6
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl7= clust7b[mask]
lisv7 = sv7[mask]
linV7 = nV7[mask]
lidist7=dists7[mask]
liid7=id7[mask]
write_clusters(licl7,lisv7,linV7,lidist7,liid7,"clusters.2-7b_li1791")
# cp clusters.2-7b_li1791 clusters.out 
# Run mode 44 uncle to make pi file
m = readdlm("pimat.2-7b_li1791")
rank(m)
#figs7b = rd_list2(figs7bt[1:100_000],rots)
# end
@time pts = make_eqvPoints(Rpts[1:40],rots)
figs3bt = make_figure_candidates(pts,3)
figs3b = reduce_figList(figs3bt[1:100_000],rots)
n3b = length(figs3b)
clust3b = vcat(licl,figs3b)
sv3 = vcat(lisv,[[1,1,1] for i ∈ 1:n3b])
id3 = vcat(liid,1:n3b)
dists3 = vcat(lidist,MADfromCOM.(figs3b))
nV3 = vcat(linv,ones(Int,n3b)*3)
write_clusters(clust3b, sv3, nV3, dists3, id3, "clusters.2-3b_"*string(n3b))
#>  cp clusters.2-3b_6111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-3b_"*string(n3b))
rank(m) # 6 (so added two 3-body clusters)
_,R3 = qr(m)
mask = get_leftmost_indep_columns(m)
rank(m) #6

popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl3= clust3b[mask]
lisv3 = sv3[mask]
linV3 = nV3[mask]
lidist3=dists3[mask]
liid3=id3[mask]
write_clusters(licl3,lisv3,linV3,lidist3,liid3,"clusters.2-3b_li6")

# Get 4 bodies
#licl3,lisv3,linV3,lidists3,liid3 = read_clusters_from_file("clusters.2-3b_li6")
pts = make_eqvPoints(Rpts[1:12],rots)
figs4bt = make_figure_candidates(pts,4)
figs4b = reduce_figList(figs4bt[1:200_000],rots)
n4b = length(figs4b)
clust4b = vcat(licl3,figs4b)
sv4 = vcat(lisv3,[[1,1,1,1] for i ∈ 1:n4b])
id4 = vcat(liid3,1:n4b)
dists4 = vcat(lidists3,MADfromCOM.(figs4b))
nV4 = vcat(linV3,ones(Int,n4b)*4)
write_clusters(clust4b, sv4, nV4, dists4, id4, "clusters.2-4b_"*string(n4b))
#>  cp clusters.2-4b_6111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-4b_"*string(n4b))
rank(m) # 11 (so added five 4-body clusters)
_,R4 = qr(m)
mask = get_leftmost_indep_columns(m)
rank(m) # 11 but it really seems there is an epsilon issue and its really 10
# meaning that four 4-body clusters were added to the list
mask = mask[1:10]
rank(m[:,mask])
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl4= clust4b[mask]
lisv4 = sv4[mask]
linV4 = nV4[mask]
lidist4=dists4[mask]
liid4=id4[mask]
write_clusters(licl4,lisv4,linV4,lidist4,liid4,"clusters.2-4b_li10")

# Get 5 bodies
#licl3,lisv3,linV3,lidists3,liid3 = read_clusters_from_file("clusters.2-3b_li6")
pts = make_eqvPoints(Rpts[1:6],rots)
figs5bt = make_figure_candidates(pts,5)
figs5b = reduce_figList(figs5bt[1:150_000],rots)
n5b = length(figs5b)
clust5b = vcat(licl4,figs5b)
sv5 = vcat(lisv4,[[1,1,1,1,1] for i ∈ 1:n5b])
id5 = vcat(liid4,1:n5b)
dists5 = vcat(lidist4,MADfromCOM.(figs5b))
nV5 = vcat(linV4,ones(Int,n5b)*5)
write_clusters(clust5b, sv5, nV5, dists5, id5, "clusters.2-5b_"*string(n5b))
#>  cp clusters.2-5b_6111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 55
m = readdlm("pimat.2-5b_"*string(n5b))
rank(m) # 12 (so added two 5-body clusters)
_,R5 = qr(m)
mask = get_leftmost_indep_columns(m)
rank(m) #12, adding two 5-body clusters
rank(m[:,mask])
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl5= clust5b[mask]
lisv5 = sv5[mask]
linV5 = nV5[mask]
lidist5=dists5[mask]
liid5=id5[mask]
write_clusters(licl5,lisv5,linV5,lidist5,liid5,"clusters.2-5b_li12")

# Get 6 bodies
#licl3,lisv3,linV3,lidists3,liid3 = read_clusters_from_file("clusters.2-3b_li6")
pts = make_eqvPoints(Rpts[1:5],rots)
figs6bt = make_figure_candidates(pts,6)
figs6b = reduce_figList(figs6bt[1:200_000],rots)
n6b = length(figs6b)
clust6b = vcat(licl5,figs6b)
sv6 = vcat(lisv5,[[1,1,1,1,1,1] for i ∈ 1:n6b])
id6 = vcat(liid5,1:n6b)
dists6 = vcat(lidist5,MADfromCOM.(figs6b))
nV6 = vcat(linV5,ones(Int,n6b)*6)
write_clusters(clust6b, sv6, nV6, dists6, id6, "clusters.2-6b_"*string(n6b))
#>  cp clusters.2-6b_6111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-6b_"*string(n6b))
rank(m) # 14 (so added two 6-body clusters)
_,R6 = qr(m)
mask = get_leftmost_indep_columns(m)
# something went wrong and the routine returned ~6600 elements. 
rank(m) #14
rank(m[:,mask])
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl6= clust6b[mask]
lisv6 = sv6[mask]
linV6 = nV6[mask]
lidist6=dists6[mask]
liid6=id6[mask]
write_clusters(licl6,lisv6,linV6,lidist6,liid6,"clusters.2-6b_li14")

# Get 7 bodies
#licl3,lisv3,linV3,lidists3,liid3 = read_clusters_from_file("clusters.2-3b_li6")
pts = make_eqvPoints(Rpts[1:3],rots)
figs7bt = make_figure_candidates(pts,7)
figs7b = reduce_figList(figs7bt[1:150_000],rots)
n7b = length(figs7b)
clust7b = vcat(licl6,figs7b)
sv7 = vcat(lisv6,[[1,1,1,1,1,1,1] for i ∈ 1:n7b])
id7 = vcat(liid6,1:n7b)
dists7 = vcat(lidist6,MADfromCOM.(figs7b))
nV7 = vcat(linV6,ones(Int,n7b)*7)
write_clusters(clust7b, sv7, nV7, dists7, id7, "clusters.2-7b_"*string(n7b))
#>  cp clusters.2-7b_7111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-7b_"*string(n7b))
rank(m) # 15 (so added one 7-body clusters)
_,R7 = qr(m)
mask = get_leftmost_indep_columns(m)
# something went wrong and the routine returned ~7700 elements. 
rank(m) #15
rank(m[:,mask])
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl7= clust7b[mask]
lisv7 = sv7[mask]
linV7 = nV7[mask]
lidist7=dists7[mask]
liid7=id7[mask]
write_clusters(licl7,lisv7,linV7,lidist7,liid7,"clusters.2-7b_li15")

# Get 8 bodies
#licl3,lisv3,linV3,lidists3,liid3 = read_clusters_from_file("clusters.2-3b_li7")
pts = make_eqvPoints(Rpts[1:2],rots)
figs8bt = make_figure_candidates(pts,8)
figs8b = reduce_figList(figs8bt[1:100],rots)
n8b = length(figs8b)
clust8b = vcat(licl7,figs8b)
sv8 = vcat(lisv7,[[1,1,1,1,1,1,1,1] for i ∈ 1:n8b])
id8 = vcat(liid7,1:n8b)
dists8 = vcat(lidist7,MADfromCOM.(figs8b))
nV8 = vcat(linV7,ones(Int,n8b)*8)
write_clusters(clust8b, sv8, nV8, dists8, id8, "clusters.2-8b_"*string(n8b))
#>  cp clusters.2-8b_8111 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-8b_"*string(n8b))
rank(m) # 12 (so added two 8-body clusters)
t_,R8 = qr(m)
mask = get_leftmost_indep_columns(m)
# something went wrong and the routine returned ~8800 elements. 
rank(m) #16
rank(m[:,mask])
popfirst!(mask) # Get rid of empty-cluster; it isn't in the clusters file
mask .-= 1 #The clusters index and pi matrix index are shifted by 1

licl8= clust8b[mask]
lisv8 = sv8[mask]
linV8 = nV8[mask]
lidist8=dists8[mask]
liid8=id8[mask]
write_clusters(licl8,lisv8,linV8,lidist8,liid8,"clusters.2-8b_li16")
#>  cp clusters.2-8b_li16 clusters.out 
#>  rm enum_PI_matrix.out; ../../uncle/src/uncle.x 44
m = readdlm("pimat.2-8b_li16")
rank(m) # 12 (so added two 8-body clusters)




pwd()

2

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
# # Generate 7-body clusters from first three neighbors

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


## Get a full set of clusters, with s-vectors for fcc n=1-6 *ternary* case
cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data")
m=readdlm("pimat.1-6.ternary")
idx = get_leftmost_indep_columns(m)
idx .-= 1
idx = idx[2:end]
cl,sv,nV,rad,id = read_clusters_from_file("clusters.1-6.ternary_overcomplete")
write_clusters(cl[idx],sv[idx],nV[idx],rad[idx],id[idx],"clusters.fcc_ternary.1-6_1081")
cl,sv,nV,rad,id = read_clusters_from_file("clusters.fcc_ternary.1-6_1081")
m = readdlm("pimat.li_1-6.ternary")
en = readdlm("energiesPerAtom_ternary1-6.dat")
J = m\en 
colors = replace(nV, 1=>:red, 2=>:orange,3=>:green,4=>:magenta,5=>:blue,6=>:red)
pushfirst!(colors,:black)
plot(J,st=:scatter,msw=0,color=colors,ms=2,legend=:none,xlabel="Cluster Number",ylabel="J Coefficient Size (ECIs)")
plot(abs.(J),st=:scatter,msw=0,color=colors,ms=2,legend=:none,xlabel="Cluster Number",ylabel="J Coefficient Size (ECIs)",yaxis=:log)
plot(abs.(m*J-en))
norm(m*J-en)
c2 = colors[idx]


# Should be able to start from scratch from here
m=readdlm("pimat.li_1-6.ternary")
p = Int.(readdlm("reorderCl_1081.dat.2"))
#m = m[:,[p...]]
enpa = readdlm("energiesPerAtom_ternary1-6.dat")
s = [1:22;63:80;150:1081;23:62;81:149]
p = p[s]
m = m[:,[p...]]

# This function now in ClusterEnumeration.jl, "sweep_model_sizes"
begin
Nits = 20
sizes = collect(1:1:400) 
#sizes = collect(1:1:1081) 
data = map(sizes) do ms
eFit = 0
eVal =0
rav = 0
println(ms)
for i = 1:Nits 
t = randperm(size(m,1))
nFit = 350
fitIdx = t[1:nFit]
valIdx = t[nFit+1:end]
J = m[fitIdx,1:ms]\enpa[fitIdx]
eFit += norm(m[fitIdx,1:ms]*J-enpa[fitIdx])/sqrt(nFit)
rav += rank(m[fitIdx,1:ms])
eVal += norm(m[valIdx,1:ms]*J-enpa[valIdx])/sqrt(size(m,1)-nFit)
end
eFit/Nits,rav/Nits, eVal/Nits
end 
errFit = [i[1] for i ∈ data]
ranks = [i[2] for i ∈ data]
errVal = [i[3] for i ∈ data]
#plot(sizes,errFit,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Validation Error",legend=:none,color=colors)
end

plot(sizes,errFit,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Fit Error",legend=:none,color=colors[p][sizes],yrange=(1e-5,3e-1),xticks = 0:50:1100,xrange=(0,600))
plot(sizes, errVal,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Val Error",legend=:none,color=colors[p][sizes],yrange=(1e-2,1e0),#ytick=[3e-2,5e-2,.0,1,1e1,1e2]
#xrange=(1,80) 
)
plot!(sizes,[i[2] for i in errFit])
plot(sizes,ranks,st=:scatter,msw=0,ms=2,ylabel="Rank of Design Matrix",xlabel="Model Size",xrange=(0,450),legend=:none)
plot!([0; 450],[0; 450],color=:red)

plot(sizes,errVal,yaxis=:log,st=:scatter,msw=0,ms=2,xlabel="Model size (num parameters)",ylabel="Avg Validation Error",legend=:none,color=colors[sizes],yrange=(4e-2,8e-0),xticks = 0:50:1100,xrange=(0,600))

# Argh, energies are not per atom
#
# tail -n +17 struct_enum.out.ternary|awk '{print $7}'
perAtom = readdlm("perAtom.dat")
writedlm("energiesPerAtom_ternary1-6.dat",en./perAtom)
enpa = en./perAtom


fitErr = map(1:1080) do x
    J = m[:,1:x]\enpa
    sum(abs.(m[:,1:x]*J-enpa))/x
    #norm(m[:,1:x]*J-enpa)/x
end 
plot(1:1080,valErr,st=:scatter,color=colors,xlabel="Model size",ylabel="Avg validation error",legend=:none,msw=0,ms=2,yaxis=:log)
savefig(ans,"fittingErrors.pdf")
end
