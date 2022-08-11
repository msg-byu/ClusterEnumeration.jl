module ClusterEnumeration
using Spacey
using LinearAlgebra
using Combinatorics
using Printf
using DataStructures

export gen_points_in_supercell, read_lattice_vectors, genSymEqvPoints, testDegeneracies, genLatticePts, makeFigureCandidates, diameter, symReduceFigList, makeFullFigsFromShellPts, augmentFigures, read_clusters_from_file, write_clusters, iterAugmentFigures, get_leftmost_indep_columns, get_nonzero_index, symReduceByLengthFigList, isMatrixEqvBySymmetry

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

function symReduceByLengthFigList(figlist,rots)
    norms = round.(diameter.(figlist),digits=12)
    p = sortperm(norms)
    figlist = figlist[p]
    norms = norms[p]
    #println.(norms)
    uqNorms = unique(norms)
    nN = length(uqNorms)
    d = [findfirst(i.==norms) for i ∈ uqNorms]
    sd = [findlast(i.==norms) for i ∈ uqNorms]
    #println.(sd)
    uqFigs = []
    println("Number of norms to search: ",nN)
    println("Total number to search: ",length(figlist))
    for iNorm ∈ 1:nN
        #println("iNorm Norm: ",uqNorms[iNorm])
        if mod(iNorm,500)==0 println("Norm #: ",iNorm) end
        uqAtLen = []
        for id ∈ d[iNorm]:sd[iNorm]
            #println("id:",id)
            uq = true
            iFig = figlist[id]
            #println("id diameter: ",diameter(iFig))
            for jFig ∈ uqAtLen
                if isMatrixEqvBySymmetry(iFig,jFig,rots)
                    uq = false
                    #println("Found eqv")
                    break
                end
            end
            if uq 
                #println("** saving **")
                push!(uqAtLen,iFig)
            end
        end
        append!(uqFigs,uqAtLen)
    end
    println("Found ",length(uqFigs))
    return uqFigs
end

""" leftmost_QR(m)
Get the leftmost columns of m that are linearly independent. Do this using QR with pivoting and then throwing away the the rightmost column and searching all of m to the left of that. Repeat until the rightmost column thrown away has no dependent copies to its left.
"""
function leftmost_QR(m)
    
end


""" Find the non-zero elements on the diagonal of a matrix """
function get_nonzero_index(m,reps=1e-13)
    mask = abs.(diag(m)).>reps
    return mask
end

""" Use QR decomposition iteratively to pull out left-most independent columms """
function get_leftmost_indep_columns(m,maxits=80)
    r = rank(m)
    nr,nc = size(m)
    idx = collect(1:nc)
    it = 1
    mask = trues(nc)
    while true
        println("it: ",it)
        println("len idx: ",length(idx))
        if it > maxits break end
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

""" augmentFigures(figs,sphPts)
Generate all (N+1)-vertex clusters from N-vertex clusters by adding each sphere point to each cluster. Doesn't eliminate duplicates
"""
function augmentFigures(figs,sphPts)
    # The if statement prevents a point being added that is already one of the vertices.
    newclust = [sortslices(hcat(j,i),dims=2) for i in sphPts for j in figs if !any(sum(j.-i.≈zeros(3),dims=1).==3)]
    return newclust
end

""" iterAugmentFigures(figs,sphPts)
Generate all (N+1)-vertex clusters from N-vertex clusters by adding each sphere point to each cluster. Doesn't eliminate duplicates
"""
function iterAugmentFigures(figs,sphPts)
    rcut = maximum(diameter.(figs))
    newFigs = []
    c = 0
    t = round(Int,length(figs)*length(sphPts)/100)
    println("t: ",t)
    for j ∈ figs
        for i ∈ sphPts
            if any(sum(j.-i.==zeros(Int,3),dims=1).==3) continue end
            newfig = hcat(j,i)
            c +=1
            if diameter(newfig) > rcut continue end
            push!(newFigs,sortslices(newfig,dims=2))
            if mod(c,t)==0 print(100*round(c/t,digits=2),"%  ") end
        end
    end
    println("count: ",c)
    return newFigs
end

""" makeFullFigsFromShellPts(shellPts,rotations,vertexorder)
Generate symmetrically distinct figures from a list of representative points of the lattice. This can be really slow for large numbers of pts or vertex orders bigger than 2. Start small.
"""
function makeFullFigsFromShellPts(Rpts,rots,nV)
    sphPts,_,_ = genSymEqvPoints(Rpts,rots)
    CandFigs = makeFigureCandidates(sphPts,nV)
    reducedFigs = symReduceFigList(CandFigs,rots)
    return reducedFigs
end

""" shiftFigures(figs_list)
Shift the origin of figures so that the first vertex is [0,0,0]. Subtract the first vertex from all the rest.
"""
function(figs_list)
    shiftList =  [m .- m[:,1] for m in figs_list]
    return shiftedList
end

""" symReduceFigList(fig_list,rots)

Remove redundant figures from the list using translations and rotations
"""
function symReduceFigList(figList,rots)
    radii = diameter.(figList)
    p = sortperm(radii)
    figList = figList[p]
    radii = radii[p]
    println("Uq radii: ", length(unique(round.(radii,digits=12))))
    uqFigs = Vector{Matrix{Int}}()
    for (i,iFig) ∈ enumerate(figList)
        uqFlag = true
        if mod(i,10_000)==0 println(i," ",length(uqFigs)) end
        for jFig ∈ uqFigs
            rj = diameter(jFig)
            if radii[i]≉rj && radii[i] > rj
                continue
            elseif radii[i]≉ rj && radii[i] < rj
                break
            end
            if isMatrixEqvBySymmetry(iFig,jFig,rots)
                uqFlag = false
                break
            end
            # for iRot ∈ rots
            #     t = sortslices(replace!(iRot*iFig,-0.0=>0.0),dims=2)
            #     #t = sortslices(iRot*iFig,dims=2)
            #     r = t - jFig # Compute translation difference between the vertices
            #     if t ≈ jFig # then iFig is a rotation duplicate of jFig
            #         uqFlag = false
            #         break
            #     elseif all(y->y≈r[:,1],eachcol(r)) # Then iFig is a rotated translation duplicate
            #         uqFlag = false
            #         break
            #     end
            # end
            # if uqFlag==false break end
        end
        if uqFlag
            push!(uqFigs,iFig)
        end
    end
    return uqFigs
end


""" diameter(figure)
Calculate diameter of a figure """
function diameter(fig)
    nV = size(fig,2) # Number of vertices
    com = sum(fig,dims=2)./nV # Center of mass
    diffs = [2*norm(i.-com) for i ∈ eachcol(fig)] # Dists from center of mass, x2
    return sum(diffs)./nV # Average of dists
end

""" make_figure_candidates(sph_point_list, n)

Make figures (with 'n' vertices, n > 1) from list of points. This routine does not remove translationally or rotationally equivalent figures.
"""
function makeFigureCandidates(pts, n)
    figs = collect(combinations(pts,n-1))
    [pushfirst!(ifig,[0,0,0]) for ifig ∈ figs]
    figs = [hcat(ifig...) for ifig ∈ figs] # Convert vecs of vecs to matrices
    figs = [sortslices(i,dims=2) for i ∈ figs] # Sort the vertices
    # We don't want to keep any candidate with a diameter bigger than the radius of the sphere
    d = diameter.(figs)
    rcut = maximum(d)/1.2
    idx = findall(x-> x<rcut, d)
    figs = figs[idx] # Return a list, sorted by diameters.
    figs = figs[sortperm(diameter.(figs))]
end

# """ (iter version) Make figures (with 'n' vertices, n > 1) from list of points. This routine does not remove translationally or rotationally equivalent figures.
# """
# function iterMakeFigureCandidates(pts, n)
#     rcut = maximum(norm.(pts))/2
#     println("cutoff:",rcut)
#     figs = combinations(pts,n-1)
#     candfigs = []
#     for ifig ∈ figs
#         pushfirst!(ifig,[0,0,0])
#         ifig = hcat(ifig...)
#         ifig = sortslices(ifig,dims=2)
#         if diameter(ifig) > rcut continue end
#         push!(candfigs,ifig)
#     end
#     p = sortperm(diameter.(candfigs))
#     candfigs = candfigs[p]
# end


""" Test that all degeneracies are 48, 12, 6, or 24 """
function testDegeneracies(list)
    mask = mod.(list,48).==0 .|| falses(length(list))
    mask =  mod.(list,24).==0 .|| mask
    mask =  mod.(list,6).==0 .|| mask
    mask =  mod.(list,12).==0 .|| mask
    mask =  mod.(list,8).==0 .|| mask
    if !all(mask) error("Bad degeneracies") end
    return true
end

""" Generate list of symmetry-equivalent points from list of points """
function genSymEqvPoints(pointList,rots)
    pts = []
    degen = []
    lengths = []
    for i ∈ pointList
        uqPts = unique([irot*i for irot ∈ rots])
        push!(pts,uqPts...)
        push!(degen,length(uqPts))
        push!(lengths,norm(i))
    end
    if !testDegeneracies(degen) error("Degeneracy error") end
    return pts,degen,lengths
end

""" Read in the lattice vectors from lat.in """
function read_lattice_vectors()
    g = readlines("lat.in")
    A = hcat([parse.(Float64,split(i)) for i ∈ g[6:8]]...)
    return A
end

""" genLatticePts(A,Rots,multiple)
Generate lattice points in a sphere. Return points and their distances
"""
function genLatticePts(A,rots,l)
    vol = abs(det(A))
    pts =[A*[i, j, k] for i ∈ -l:l for j ∈ -l:l for k ∈ -l:l]
    # Sort the points to speed up the elimination of duplicates
    pts = pts[sortperm([norm(i) for i ∈ pts])]
    uqPts = []
    pcount = 1000
    for iPt ∈ pts
        len = length(uqPts)
        uq = true
        if norm(iPt) > √2*l continue end # This works for fcc...the extra factor of two integer fcc
        if mod(len,pcount)==0 println("pts so far: ",len); pcount+=1000 end
        for jPt ∈ uqPts
            if norm(iPt) > norm(uqPts[end]) break end # should speed up things
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
    return uqPts[2:end],norm.(uqPts[2:end]) # Skip point at origin
end


""" gen_points_in_supercell(A,rcut)
Generate points on the lattice defined by A that are shorter than rcut.
This doesn't work robustly yet. Will provide good results only in ideal cases.
"""
function gen_points_in_supercell(A,l)
    u,v,w = eachcol(A)
    vol = abs(u×v⋅w)
    println("Generating approximately ",round(Int,(l/√2)^3/3))
    _, rots = pointGroup_robust(u,v,w)
    #s = ceil(Int,∛npts)
    #l,m,n = [s,s,s]
    pts =[A*[i, j, k] for i ∈ -l:l for j ∈ -l:l for k ∈ -l:l]
    pts = pts[sortperm([norm(i) for i ∈ pts])]
    uqPts = []
    for iPt ∈ pts
        pcount = 1000
        len = length(uqPts)
        uq = true
        if norm(iPt) > √2/2*l continue end # This works for fcc...
        if mod(len,pcount)==0 println("pts so far: ",len); pcount+=1000 end
        for jPt ∈ uqPts
            if norm(iPt) > norm(uqPts[end]) break end # should speed up things
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
    return uqPts[2:end],norm.(uqPts[2:end]) # Skip point at origin
end

end # module
