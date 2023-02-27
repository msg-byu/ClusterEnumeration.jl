vibup to 3-body clusters in pool
| max cell size | Num structs | PI rank |
| ---           |          ---| ---     |
| 2             | 4           | 4       |
| 3             | 10          | 10      | round3:
| 4             | 29          | 22      |  23!
| 5             | 57          | 40      |  41!
| 6             | 137         | 70      |  70
| 7             | 241         | 115     |  115
| 8             | 631         | 191     |  191
| 9             | 1135        | 286     |  291!
| 10            | 2346        | 415(418)|  429! (424 when I re-ran. Numerical issues?) Mon: 426

up to 4-body clusters in pool
| max cell size | Num structs | PI rank |
| ---           |          ---| ---     |
| 2             | 4           | 4       |
| 3             | 10          | 10      |
| 4             | 29          | 29      |
| 5             | 57          | 53      | (52)
| 6             | 137         | 112     | (112)
| 7             | 241         | 180     | (183)
| 8             | 631         | 356     | (364)
| 9             | 1135        | 550     | (559)
| 10            | 2346        | 899     | (924) (930 on rerun) Mon: 934, 940

up to 5-body clusters in pool
| max cell size | Num structs | PI rank |
| ---           |          ---| ---     |
| 2             | 4           | 4       |
| 3             | 10          | 10      |
| 4             | 29          | 29      |
| 5             | 57          | 57      |
| 6             | 137         | 126     | (126)  (127) Monday: 127, 127
| 7             | 241         | 207     | (211)  (213)
| 8             | 631         | 459     | (468)  (471)
| 9             | 1135        | 750     | (758)  (770)
| 10            | 2346        | 1300    | (1322) (1349) Monday: 1354, 1364

up to 6-body clusters in pool
| max cell size | Num structs | PI rank |
| ---           |          ---| ---     |
| 2             | 4           | 4       |
| 3             | 10          | 10      |
| 4             | 29          | 29      |
| 5             | 57          | 57      | (57)
| 6             | 137         | 134!?   | (because three bodies were missing? That turned up on round 2. Still 134! Rats)
Which *rows* are degenerate? Look at those structures and see if that gives any insight.
What linear combination are those rows? What combination of the other nows make them?
A_indpendent * combi = dependent structure (use Pi matrix transpose) combare the DFT energy of the raw structure to the linear combo of the energies of those that make it up.

Could also do the same with the columns. And check whether the dependent column is "close" to an important column (short pair, e.g.) or something unimportant.
                                                 Monday: 135 (got one more!), Now 136!
| 7             | 241         |         |  (224) Tuesday 227
| 8             | 631         |         |  (546)         556
| 9             | 1135        |         |  (900)         926
| 10            | 2346        |         |  (1655) (1666) (1689) Monday: 1694, 1727, Tues: 1731

Other ideas:
Keep a few of the "small" diags at each order. (Maybe look at the singular values or size of R diag entries)

Plot dependent columns with close singular values
Look at distance vs singular value size


Gus Hart
8 month base: $167,186 ($20898/month)
10 month base: $197,289