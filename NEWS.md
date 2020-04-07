# LMN 1.1

## Changes to interface

* All `lmn.*()` functions converted to `lmn_*()`.

# LMN 0.0.1.9001

## Changes to interface

* `lmn.suff()` specified variance via `V` and `Vtype`.  So e.g., `lmn.suff(acf = acf)` becomes `lmn.suff(V = acf, Vtype = "acf")`.
* All `lmn.*` functions except `lmn.suff` now only accept a `suff` argument instead of `Y`, `X`, `V`, etc.
* `lmn.suff()` return element `Beta.hat` renamed `Bhat`.

# LMN 0.0.1.9000

* Added `Hospitals` dataset.
