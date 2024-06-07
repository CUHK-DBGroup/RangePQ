<p align="center">
  <img src="RangePQ3.png" width="320">
</p>

# Efficient Dynamic Indexing for Range Filtered Approximate Nearest Neighbor Search

This repository implements a dynamic indexing (RangePQ+) for range filtered approximate nearest neighbor search. Given an object set, each object represents a vector with an associated attribute value. Range filtered ANN search contains a query range $Q=[l,r]$. The query can be viewed as, ANN search only for those vectors whose attribute values are in the particular range $Q$.

## Features of RangePQ+
- Efficient range filtered ANN.
- RangePQ+ support efficiently insert and delete object.
- RangePQ+ is a lightweight index scheme.


## Compiling
To compile the C++ file, you can run 
```
./build.sh
```

 After that, you can directly use our RangePQ index in Python by

 ```python
import segr

```

## Experiments details
For our RangePQ index, we set $L_base$ = 10000 to test on the Sift and WIT, and $L_base$ = 30000 to test on the Gist.

For Rii, we set $L$ = 1000 to test on the Sift and WIT, and $L$ = 3000 to test on the Gist.

For Milvus, we set $nlist$ = 6 to test on the Sift,  $nlist$ = 15 to test on the Gist, and $nlist$ = 6 to test on the WIT.

For Vbase, we set $nlist$ = 5 to test on the Sift,  $nlist$ = 20 to test on the Gist, and $nlist$ = 8 to test on the WIT.
