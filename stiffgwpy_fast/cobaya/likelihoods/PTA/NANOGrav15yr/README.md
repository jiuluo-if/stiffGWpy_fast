# NANOGrav 15-year mock resources

[中文说明](README_zh.md) | English

`../NANOGrav.yaml` selects three NumPy resources in this folder:

- `density_mock.npy`: per-frequency likelihood-density samples.
- `log10rhogrid.npy`: the `log10(rho)` grid for those likelihood densities.
- `freqs.npy`: PTA frequency bins in Hz.

The adapter loads these arrays from its configured paths and uses `Nfreqs` to
select the bins for evaluation. Preserve the file names and their pairing with
the YAML configuration. This guide makes no new claim about data provenance.
