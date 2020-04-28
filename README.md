# Primal-Dual Interior-Point Solver (PDIPS)

[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/sinha-abhi/PDIPS.jl/blob/master/LICENSE)

`PDIPS` is an implementation of the homogeneous self-dual interior-point
algorithm for solving linear optimization problems using interior-point methods.

This implementation solves linear programs of the form
```
min     c'x
 x
s.t    Ax = b
    lo <= x <= hi
```
by transforming it into standard form:
```
min     c'x
 x
s.t    Ax = b
       0 <= x.
```


#### References

For a detailed explanation of the algorithm see
[Tulip.jl](https://github.com/ds4dm/Tulip.jl).
