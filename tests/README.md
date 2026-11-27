# Core tests

The test suite focuses on package behavior: Fast and reference solver contracts, physics limits, frequency integration, compatibility, and the two fixed numerical reference regressions in `data/`.

The compact helpers in `support/` are imported only by core solver tests; they are not command-line research drivers and are excluded from distributions.

Run the default suite with:

```bash
python -m pytest -q
```
