# LVK cross-correlation likelihood

[中文说明](README_zh.md) | English

`LVK_SGWB_CC.py` provides the `LVK_SGWB_CC` Cobaya likelihood. Defaults are
defined in `LVK_SGWB_CC.yaml`; `CC_file` points to `C_O1_O2_O3.dat`.

The adapter reads three columns: frequency, measured cross-correlation, and
uncertainty. It converts frequency to `log10(Hz)`, interpolates the theory
spectrum on supported bins, and evaluates a Gaussian chi-square log-likelihood.
The likelihood alias is set in the YAML file.

See the [likelihood package index](../README.md), the
[Cobaya theory guide](../../../../docs/cobaya.md), and the
[root README](../../../../README.md). This page describes the code and file
layout; it makes no separate claim about data provenance or model fit.
