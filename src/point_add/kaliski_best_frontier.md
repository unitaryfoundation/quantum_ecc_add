# Current strict best frontier

Using the explicit specialized backward (not the failed generalized reverse),
the current best strict `main.rs`-validated passing settings found so far are:

| k | avg executed Toffoli |
|---|---:|
| 255 | 4,261,946 |
| 272 | 4,247,076 |
| 288 | 4,238,868 |
| 297 | 4,236,726 |
| 299 | 4,236,492 |
| 300 | 4,236,408 |
| 302 | 4,236,306 |
| 304 | **4,236,292** |
| 305 | 4,236,318 |
| 310 | 4,236,778 |
| 312 | 4,237,116 |
| 313 | 4,237,318 |
| 318 | 4,238,658 |
| 319 | 4,238,992 |
| 323 | 4,240,548 |
| 324 | 4,240,992 |
| 326 | 4,241,946 |

## Best strict result so far
- with `KAL_BULK3_EXPERIMENT=1`
- and default `BULK_PREFIX_SAFE_ITERS = 304`
- avg executed Toffoli = **4,236,292**
- savings vs baseline `4,394,546` = **158,254**

## Shape of the frontier
The frontier remains highly nonmonotone.
There is a strong passing cluster around `k ≈ 300`, then failures quickly blow
up beyond the low 330s.
