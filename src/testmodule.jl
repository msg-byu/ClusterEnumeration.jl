module testmodule
using Printf 

export write_clusters, read_lattice_vectors


""" Read in the lattice vectors from lat.in """
function read_lattice_vectors()
    g = readlines("lat.in")
    A = hcat([parse.(Float64,split(i)) for i ∈ g[6:8]]...)
    return A
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

end
