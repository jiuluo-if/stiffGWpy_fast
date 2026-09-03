# Extended matched fast-vs-reference spots (parameter edges + interiors)

Production fast engine vs continuous-sigma DOP853 reference on the SAME grid-independent frequency grid, SAME `z_tail=8`, rtol=1e-9, freq_res=1.0. Points cover parameter-axis edges (u ~ 0.02/0.98 of the sampled box) and transition-sensitive interiors.

| label | r | n_t | cr | T_re | DN_re | kappa10 | sig rel max | DN rel | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| edge_dnre_hi | 1.000e-02 | 0.00 | 0 | 2.000e+03 | 29.4 | 1.000e-02 | 5.892e-04 | +1.116e-04 | ok |
| edge_dnre_lo | 1.000e-02 | 0.00 | 0 | 2.000e+03 | 0.6 | 1.000e-02 | 6.432e-04 | -2.069e-04 | ok |
| edge_kap_hi | 1.000e-02 | - | 1 | 2.000e+03 | - | 7.600e-01 | 6.747e-04 | -9.339e-04 | ok |
| edge_kap_lo | 1.000e-02 | - | 1 | 2.000e+03 | - | 1.300e-06 | 6.139e-04 | -2.724e-05 | ok |
| edge_nt_red | 1.000e-02 | -0.48 | 0 | 2.000e+03 | - | 1.000e-02 | 6.997e-04 | +1.918e-02 | ok |
| edge_r_hi | 7.940e-02 | - | 1 | 2.000e+03 | - | 1.000e-02 | 1.641e-03 | -4.458e-04 | ok |
| edge_r_lo | 1.260e-06 | - | 1 | 2.000e+03 | - | 1.000e-02 | 6.991e-04 | -2.863e-04 | ok |
| edge_tre_lo | 1.000e-02 | - | 1 | 1.260e+01 | - | 1.000e-02 | 6.515e-04 | +1.064e-03 | ok |
| interior_dnre20 | 3.000e-02 | 0.00 | 0 | 2.000e+03 | 20.0 | 1.000e-02 | 6.655e-04 | -6.883e-04 | ok |
| interior_kap03 | 2.000e-02 | - | 1 | 1.000e+03 | - | 3.000e-01 | 6.883e-04 | -6.002e-04 | ok |
| interior_ntred | 1.000e-02 | -0.30 | 0 | 2.000e+03 | - | 1.000e-02 | 6.997e-04 | +4.618e-03 | ok |
| interior_r05 | 5.000e-02 | - | 1 | 2.000e+03 | - | 1.000e-02 | 6.772e-04 | +7.672e-05 | ok |
| interior_tre1000_r5e3 | 5.000e-03 | - | 1 | 1.000e+03 | - | 5.000e-02 | 6.825e-04 | -7.918e-04 | ok |
| interior_tre300 | 3.000e-02 | - | 1 | 3.000e+02 | - | 1.000e-01 | 6.859e-04 | -2.675e-04 | ok |

Aggregates over 14 ok points: signal rel max **1.641e-03** (gate <1e-3: FAIL), transition rel max **1.641e-03** (gate <1e-3: FAIL), DN rel abs median **5.230e-04** / max **1.918e-02** (gate <1e-4: FAIL).

Non-ok records: 2 total (2 explicit shared_Neff_guard, others exception/rejected recorded per point; never silent).
