using Plots
using DelimitedFiles

default(ms=2,msw=0,legend=:none,st=:scatter)
cd("/Users/glh43/home/juliaCodes/ClusterEnumeration.jl/data/")
# data from supercomputer
#  for i in dir_POSCAR*;do echo $i; grep -i "free  energy" $i/OUTCAR; done > rawout.txt
reN = r"dir_POSCAR([0-9]+)"
lines = readlines("rawout.txt.urmp")
strN = [parse(Int,i[1])+1 for i in match.(reN,lines) if i!==nothing]
reE = r"TOTEN  = \s+(-[0-9.]+)"
energies = [parse(Float64,i[1]) for i in match.(reE,lines) if i!==nothing]

#  Order the energies by structure number
p = sortperm(strN)
strN = strN[p]
energies = energies[p]
plot(energies,ylabel="Total energy",xlabel="Structure number",color=:blue)

# Get the struct_enum.out "matrix"
str = readlines("struct_enum.fcc.upto10")
strlines = str[17:end]
strmat = vcat([parse.(Int,split(i))' for i in strlines]...)
conc =  count.("1",string.(strmat[:,27]))./strmat[:,7] # Number of b atoms/total atoms
enpat = energies./strmat[:,7]
plot(conc,enpat,xlabel="Pt Concentration",ylabel="Energy/atom",color=:purple,title="Convex Hull (Cu-Pt)")
#m6 = readdlm("pimat.radGen.6li_137")
#m = readdlm("pimat.2-6b_li1761")

writedlm("energiesPerAtom.fcc.upto10.urmp",enpat)
writedlm("concentration.fcc.upto10.urmp",conc)
#writedlm("struct_enum_mat.fcc.upto10",strmat)
#rank(m6)
#rank(m)
engr = readdlm("energiesPerAtom.fcc.upto10.urgr")
enmp = readdlm("energiesPerAtom.fcc.upto10.urmp")
diffs = engr-enmp
histogram(diffs,xlabel="Energy difference (eV/atom)",ylabel="Number of structures",title="autoGR vs Monkhorst-Pack (Cu-Pt, unrelaxed)",yaxis=:log)

# data from supercomputer
#  for i in dir_POSCAR*;do echo $i; grep -i "irreduc" $i/OUTCAR; done > kpts.raw.urmp
reN = r"dir_POSCAR([0-9]+)"
lines = readlines("kpts.raw.urmp")
strN = [parse(Int,i[1])+1 for i in match.(reN,lines) if i!==nothing]
reK = r" Found \s+([0-9]+)"
kptsMP = [parse(Int,i[1]) for i in match.(reK,lines) if i!==nothing]

reN = r"dir_POSCAR([0-9]+)"
lines = readlines("rawout.txt.urgr")
strN = [parse(Int,i[1])+1 for i in match.(reN,lines) if i!==nothing]
p = sortperm(strN) 
lines = readlines("kpts.raw.urgr")
reK = r"KPOINTS: autoGR kpoint generation: \s+([0-9]+)"
kptsGR = [parse(Int,i[1]) for i in match.(reK,lines) if i!==nothing][p]
