# Registration and approximation paper access audit — 2026-11-11

This update corrects older source-status summaries and records the auditable limits of the registration papers relevant to Fast. It does not change Fast or Reference code, claim a Fast validation, or replace the full scope of the research task.

## #9 — Barral, Taddei & Tifouti, registration-based model reduction

**Sources actually read:** the accessible 30-page arXiv v1 [HTML](https://arxiv.org/html/2308.01773) and [PDF](https://arxiv.org/pdf/2308.01773), plus the [publisher record](https://www.sciencedirect.com/science/article/pii/S0021999123008227). The [IRIS accepted-postprint record](https://iris.uniroma1.it/handle/11573/1745040) lists a 37-page version, but its manuscript route was blocked; that version was not read.

- §1, Eqs. (1), (3), and Algorithm 1 describe steady conservation-law snapshots, mesh adaptation/registration, and a reduced model.
- §4.1, Eqs. (26)–(27), uses a geometric map represented by a polynomial displacement, `id + Σ a_i φ_i`, with boundary constraints. The diffeomorphism facts cited there concern a restricted model class; fitted coefficients alone do not provide a certificate for every predicted parameter.
- §4.2, Eqs. (28)–(30), adds regularization for map smoothness/mesh quality/Jacobian behavior. §4.3 uses nonconvex optimization and parameter regression, so initialization and out-of-training predictions remain relevant.
- Eqs. (31)–(32) use a soft target for a detected event location. §5, Eqs. (35)–(42), describes hyper-reduced LSPG/Gauss–Newton. §5.4, Eq. (43), explicitly says its residual indicator is not a rigorous error bound and requires a high-fidelity-sized linear system.
- §7 reports steady nozzle and Gaussian-bump examples, 20 random test parameters, and a 15×15 training mesh. The reported offline training takes hours; the listed online runtimes and average HF relative `L2` errors are specific to those examples and do not establish Fast cost or root accuracy.

### Derived soft-alignment counterexample

This scalar construction is our mathematical falsifier, not a paper experiment. On `[0,1]`, let

`Φ_a(s) = s + a s(1−s)`, `|a|<1`,

so `Φ_a` is boundary preserving and strictly increasing. Let the reference event be `s*=1/2`, the desired displacement be `d>0`, and penalize the soft event target together with an `H²`-type regularizer:

`J(a) = (a/4 − d)² + 4 ξ a²`.

Its minimizer is `a = 4d/(1+64ξ)`, leaving event mismatch `−64ξd/(1+64ξ)`. At `d=0.1`, `ξ=0.001`, the map remains bijective (`min Φ′≈0.6241`) but misses the target by about `0.006015`. Thus a smooth, admissible registered map plus a soft target does not imply exact event alignment.

**Fast assessment (inference):** the method is an offline registered ROM for steady conservation laws. It does not derive the exact transformed Fast IVP, preserve both reheating-side data, or certify the self-consistent `DN_gw` root. Keep exact physical event coordinates and matching unless a separate exact change-of-variables proof is supplied.

## #11 — Tifouti, Barral & Taddei, local reduced bases

**Full-text status: partial, not complete.** The [official IRIS record](https://iris.uniroma1.it/handle/11573/1745041) lists DOI [10.2514/6.2025-1571](https://doi.org/10.2514/6.2025-1571) and a 9.49 MB postprint. Its [author PDF route](https://iris.uniroma1.it/retrieve/c79e3f36-8bb4-4dad-b7f4-d30e0f6e3bd8/Tifouti_Registration-based_2025.pdf) currently returns Cloudflare 403. An alternative [HAL record](https://inria.hal.science/hal-05065511) and [PDF](https://inria.hal.science/hal-05065511v1/file/AIAA_BTT24_draft.pdf) were found but returned an access challenge/SSL timeout. Crossref metadata uses a slightly different published title: “Local Registration-Based Model Reduction of Parameterized PDEs With Spatio-Parameter Adaptivity.”

The historical tool record from 2026-11-06 reports that a 22-page PDF was detected, but the extracted text actually available in that run ends during page 7; later pages were not opened. The recovered text is only an excerpt, not the complete paper. The verifiable excerpt includes:

- Eq. (1): a **steady** parameterized conservation law.
- Eq. (2): hard indicator-based selection among local bases; Eq. (3) uses a bijective registration; Eq. (4) gives separate local bases/maps for separate parameter subdomains.
- The authors state that the piecewise classifier is parameter-discontinuous and may be unsuitable for unsteady problems.
- Eqs. (5)–(6) describe a deformed mesh with unchanged connectivity. Algorithms 1–2 cover iterative snapshot generation, clustering, and random sampling.

Sections 3–6, numerical error/cost tables, residual guarantees, and full conclusions were not available; this note makes no claims about them.

### Derived local-map/root counterexample

This is our standalone scalar counterexample, not a paper result or Fast experiment. Let the true contraction be `G(x)=(1−δ)x+δa`, whose unique fixed point is `a`. Define a hard-classified local surrogate

`G̃(x)=G(x)−ε` for `x<a`, and `G̃(x)=G(x)+ε` for `x≥a`.

Each local formula differs from `G` by only `ε`, yet the surrogate has two fixed points `a−ε/δ` and `a+ε/δ`, both with zero surrogate residual, and a `2ε` jump at the classifier boundary. For `a=1`, `δ=0.01`, `ε=1e−5`, root error is `1e−3`, above Fast's `2e−4` root target. This shows why local fit/residual quality does not establish a continuous self-consistent root guarantee.

**Fast assessment (inference):** partitioning by known mathematical event topology can be useful. A learned classifier that switches local basis as `DN_eff` varies can instead create artificial jumps in the outer map unless cross-region continuity and error are bounded. The verified excerpt supplies neither the Fast transformed equation nor a root certificate.

## Access-status corrections for other listed works

- **#3, derivative-free uniformly accurate schemes:** the older “full author PDF unavailable” row is superseded by [`derivative_free_uniformly_accurate_audit_20261106.md`](derivative_free_uniformly_accurate_audit_20261106.md), which records a full 20-page author manuscript and theorem/algorithm/example coverage.
- **#14, Manolopoulos log-derivative method:** AIP article/PDF requests returned 403. Oxford ORA returns HTTP 200 metadata HTML with DOI links but no manuscript download. Formulas, error analysis, and experiments remain unread; applicability is undetermined from the abstract.
- **#15, Algorithm 1062:** the official [JuSER record](https://juser.fz-juelich.de/record/1055814) marks the paper OpenAccess/CC BY 4.0 and provides a [PDF route](https://juser.fz-juelich.de/record/1055814/files/3805698.pdf), but current requests time out or return a 248-byte anti-bot challenge page rather than a PDF. The [ACM full page](https://dl.acm.org/doi/full/10.1145/3805698) returns 403. However, search on the exact official DOI surfaced extensive article text from the official index. Core text was read for §§1–3, 4.1–4.4, 5–7 and part of §4.5; indexed Algorithm 1/2 and figures were only `Image` placeholders, and Eq. (46) was not surfaced. This supports a core-text audit, not a complete PDF/figure audit.

### #15 — Algorithm 1062 core-text findings and scope

- **Problem and target algorithm:** a slow, high-precision real scalar reference function is approximated by piecewise polynomials, with a domain partition chosen for fast, rounding-safe indexing and argument reduction. The target polynomial is evaluated in fixed Horner form; subrange coefficient alignment is designed around cache-line loads. The worked functions are from the complex error function library, not ODE transfer maps.
- **Approximation assumptions:** Lipschitz continuity supports convergence of nested Chebyshev interpolation, but the sharper geometric-tail estimate uses analyticity and a Bernstein-ellipse bound. Relative error also requires a domain excluding zeros (or a separate treatment near zeros).
- **Error formulas:** the indexed text splits total error into truncation and floating-point rounding. Under the coefficient envelope `|p_n| ≤ s r^n`, the truncation tail is bounded by `T = s r^(N+1)/(1−r)`. The total relative bound combines this with a Horner/coefficients rounding term and a lower bound on `|g|`; the article's worked target is close to `2ε`. It states that local argument `t` is exact and omits `O(ε²)` in the rounding analysis.
- **Critical qualification:** §3.3 says the actual `s,r` envelope estimation is heuristic: it fits candidate geometric bounds from finitely many computed coefficients. The formal infinite-tail inequality is therefore conditional on that envelope; finite coefficient observations do not by themselves certify every omitted coefficient.
- **Cost and reported tests:** the indexed text describes offline coefficient generation, a piecewise target with `2^M` subranges and fixed-degree Horner evaluation, cache alignment, and tests with 120k equidistant and 10M random points for the worked function. Those tests are paper-reported, not rerun here. The figure and Algorithm 1/2 images were not visible in the index output, and no PDF was retrieved.
- **Fast assessment (inference):** this is a possible scalar-background evaluation accelerator only after the particular function's remainder and rounding bounds are independently certified. It does not address mode propagation, frequency integration, or amplification into the self-consistent root. Its coefficient-envelope heuristic cannot substitute for the continuous parameter-domain and root-error evidence required here.

The supplied attachment itself lists 16 papers and 8 repositories and contains no 37-paper appendix. This statement is limited to the attachment that was provided.
