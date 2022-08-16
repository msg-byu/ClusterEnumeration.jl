using Plots

default(ms=2,msw=0,legend=:none,st=:scatter)
cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/")
# data from supercomputer
#  for i in dir_POSCAR*;do echo $i; grep -i "free  energy" $i/OUTCAR; done > rawout.txt
reN = r"dir_POSCAR([0-9]+)"
lines = readlines("upto10raw.data")
strN = [parse(Int,i[1])+1 for i in match.(reN,lines) if i!==nothing]
reE = r"TOTEN  = \s+(-[0-9.]+)"
energies = [parse(Float64,i[1]) for i in match.(reE,lines) if i!==nothing]
p = sortperm(strN)
strN = strN[p]
energies = energies[p]
plot(energies,ylabel="Total energy",xlabel="Structure number",color=:blue)

strlines = readlines("struct_enum.out")