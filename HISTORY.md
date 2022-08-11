*14 Jul 2022*
Started a proper repo for the code.  
Moving each block of code (each "method") into functions, so ideas are easier to explore.  

<< My current example is fcc ternary superstructures up to n=6. There are 1081 of them. >>

My initial attempts generated a PI matrix with a higher rank than the default uncle approach, but they still weren't full rank. As a debugging strategy, I'm going to run tests without trying to remove any translation or  rotation duplicates, even for identical figures with equivalent but permuted s-vectors.

** The case of 1-, 2-atom cells **
I generated 7 clusters: (1) empty cluster, (2-3) two on-site clusters, (4-7) four pair clusters (same figure but 4 s-vectors). The 1-2 and 2-1 svectors make equivalent clusters but I'm avoiding any reduction right now so that reduction bugs don't become red herrings for testing the hypothesis that I can use the input structures internal sites to generate a complete set of clusters (and get full-rank pi matrix).

When I do this, I get a full-rank 1081x6 matrix. The 6x6 matrix is also full rank.

---
** The case of 1-, 2-, 3-atom cells **  
Generated 44 2- and 3-body clusters + svecs.
There are 30 n<=3 structures so the rank should be 30. But it's only 24 in the 1081x44 matrix. Also confusing, the rank is only 20 when only the 30x44 submatrix is used. That means that the "degenerate" clusters are differentiating a few of the larger structures. That's surprising (but maybe not wrong.)
