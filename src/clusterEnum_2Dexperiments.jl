nShells = 200
rpts = [i for i in eachrow(R1[1:nShells,:])]
rads = unique(round.(norm.(rpts),digits=12))
npts = cumsum([count(r.≈norm.(make_eqvPoints(rpts,rots))) for r in rads])
N = 0:npts[end]
begin
plot(N,cbrt.(3N/16π),lw=2,color=:red,st=:line,xlabel="Number of points",ylabel="Radius",label="Continuous",legend=:topleft)
plot!(npts,rads,ms=3,color=:blue,label="Discrete")
end