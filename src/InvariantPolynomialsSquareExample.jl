using LinearAlgebra
function empty(x)
    return 1
end

function onsite(x)
    return sum(x)/length(x)
end

function nn1(x)
    idx = vcat(length(x),1:length(x)-1)
    return sum(x.*x[idx])/length(x)
end

function nn2(x)
    return (x[1]*x[3]+x[2]*x[4])/2
end

function siteSquared(x)
    return sum(x.*x)/length(x)
end
 
function siteCubed(x)
    return sum(x.^3)/length(x)
end

function threeNN(x)
    idx1 = [4,1,2,3]
    idx2 = [3,4,1,2]
    return sum(x.*x[idx1].*x[idx2])/length(x)
end

function three21(x) 
    
colorings = [
    -1 -1 -1 -1;
    -1 -1 -1  0;
    -1 -1 -1  1;
    -1 -1  0  0;
    -1 -1  0  1;
    -1 -1  1  1;
    -1  0 -1  0;
    -1  0 -1  1;
    -1  0  0  0;
    -1  0  0  1;
    -1  0  1  0;
    -1  0  1  1;
    -1  1 -1  1;
    -1  1  0  1;
    -1  1  1  1;
     0  0  0  0;
     0  0  0  1;
     0  0  1  1;
     0  1  0  1;
     0  1  1  1;
     1  1  1  1
]    

m = hcat(
    [empty(i) for i in eachrow(colorings)],
    [onsite(i) for i in eachrow(colorings)],
    [nn1(i) for i in eachrow(colorings)],
    [siteSquared(i) for i in eachrow(colorings)],
    [nn2(i) for i in eachrow(colorings)],
#    [siteCubed(i) for i in eachrow(colorings)]
    [threeNN(i) for i in eachrow(colorings)]
)
# siteCubed term will be identical to single site term


rank(m)
s = svd(m); s.S[1:end]./s.S[1]
