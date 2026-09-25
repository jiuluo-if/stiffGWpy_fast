# Packaged Cobaya likelihoods

[中文说明](README_zh.md) | English

This package contains the LVK and PTA likelihood adapters and the resources
referenced by their YAML configuration. The parent
[Cobaya adapter guide](../README.md) documents the theory interface; this page
maps likelihood classes to their configuration and data folders.

| Folder | Adapter/configuration | Packaged inputs |
|---|---|---|
| `LIGO_SGWB/` | `LVK_SGWB_CC` / `LVK_SGWB_CC.yaml` | `C_O1_O2_O3.dat` |
| `PTA/` | `IPTA` / `IPTA.yaml` | EPTA DR2 mock samples and frequency bins under `EPTAdr2/` |
| `PTA/` | `NANOGrav` / `NANOGrav.yaml` | Density and frequency grids under `NANOGrav15yr/` |

The YAML files define the resource paths, frequency counts, aliases, and
likelihood parameters. Use the folder README files to identify each input.
Data contents and provenance are not changed by this index.
