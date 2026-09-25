# Pulsar-timing-array likelihoods

[中文说明](README_zh.md) | English

This package contains IPTA and NANOGrav likelihood adapters. Their class
options and default resource paths are defined in `IPTA.yaml` and
`NANOGrav.yaml`.

| Adapter | Default inputs | Evaluation |
|---|---|---|
| `IPTA` | `EPTAdr2/EPTA_dr2new_mock.dat`, `EPTAdr2/freqs_dr2new.txt` | Uses the configured frequency bins, the predicted spectrum, and SMBHB nuisance parameters with per-bin sample-density evaluation. |
| `NANOGrav` | `NANOGrav15yr/density_mock.npy`, `NANOGrav15yr/log10rhogrid.npy`, `NANOGrav15yr/freqs.npy` | Combines the predicted spectrum with the SMBHB component and interpolates the packaged per-bin likelihood grid. |

`Nfreqs` controls the number of bins consumed by each adapter. The files are
packaged resources; use the YAML paths rather than assuming a repository-root
working directory. See the [likelihood index](../README.md) and
[`docs/cobaya.md`](../../../../docs/cobaya.md).
