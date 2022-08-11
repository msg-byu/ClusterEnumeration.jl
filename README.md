# Getting a complete cluster set
In the cluster expansion method, the clusters (geometric "figures" + basis function labels, "s-vectors") form an orthonormal basis. 

But when a finite number of clusters are used and the training set is finite, the cluster basis will no longer be orthonormal or even linearly independent in the subspace of training set structures. Clusters may be degenerate (like higher order fourier terms that are "aliased").

The purpose of this code is to generate a set of clusters that form a full-rank correlation matrix, that is every structure in the training set is linearly dependent under the choice of the cluster pool. 

The general strategy is, for each input structure, to generate all possible clusters (of every vertex order) using the lattice sites in the unit cell of the input structure. This list is concantenated for every input structure. The resulting set will generally have more clusters than necessary but will not be missing any clusters needed for a full-rank correlation matrix. (Historically, when clusters are generated from the lattice, using heuristic cut-off's etc., the correlation matrix is always rank deficient.)

The cluster pool can be "trimmed" after the full rank is verified by using only those clusters that are linearly independent.

