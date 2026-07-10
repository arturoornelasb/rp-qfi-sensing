# Q4 -- operator-structure scope law (RESULTS)

| gamma | structure | alpha | enh = a - a_dense | 1 - D2 | class |
|---|---|---|---|---|---|
| 1.4 | diag | +1.395 | +0.392 | 0.376 | INHERITS |
| 1.4 | banded | +1.219 | +0.216 | 0.376 | PARTIAL (58%) |
| 1.4 | lowrank | +1.000 | -0.002 | 0.376 | BLIND |
| 1.4 | dense | +1.003 | +0.000 | 0.376 | BLIND |
| 1.8 | diag | +1.896 | +0.921 | 0.861 | INHERITS |
| 1.8 | banded | +1.794 | +0.819 | 0.861 | INHERITS |
| 1.8 | lowrank | +0.957 | -0.018 | 0.861 | BLIND |
| 1.8 | dense | +0.975 | +0.000 | 0.861 | BLIND |

**Classification (consistent across gammas):**
- diag: **INHERITS**
- banded: **MIXED ['PARTIAL (58%)', 'INHERITS']**
- lowrank: **BLIND**
- dense: **BLIND**

Wall: 204 s.
