# Derivative-free uniformly accurate schemes: source audit for Fast

## Scope and source access

This audit covers Philippe Chartier, Mohammed Lemou, Florian Méhats, and Xiaofei Zhao, “Derivative-free high-order uniformly accurate schemes for highly-oscillatory systems,” *IMA Journal of Numerical Analysis* 42(2) (2022), 1623–1644, DOI [10.1093/imanum/drab014](https://doi.org/10.1093/imanum/drab014). The user-supplied title says “evolution equations”; the published title says “systems.” The [Oxford Academic record](https://academic.oup.com/imajna/article-abstract/42/2/1623/6260882) confirms the published title, authors, issue, pages, DOI, date history, and abstract.

Context7 resolution was attempted for the title/DOI. Context7 returned unrelated software libraries and no paper record, so it was not used as evidence. The publisher’s article page is marked “Available for Purchase”; the HAL author-manuscript URL listed by MaRDI (`hal-03141156`) returned an Anubis access-denied page in the web reader, and direct HTTPS retrieval timed out. I then opened the [20-page full-text PDF surfaced by SciSpace](https://scispace.com/pdf/derivative-free-high-order-uniformly-accurate-schemes-for-14utse9nfe.pdf). Its first pages identify HAL Id `hal-03141156`, submission date 15 February 2021, title, authors, and DOI, consistent with the official publisher record. I read the full body and appendix, including §§1–5, Theorems 2.1, 2.2, 3.1, and 4.1, the iterative algorithms, numerical examples, and the appendix proof of Theorem 3.1. The SciSpace-hosted author manuscript is the full-text source actually inspected; the publisher PDF itself was not accessible in this session.

## Paper’s problem class and hypotheses

The paper studies

\[
\dot u(t)=f_{t/\varepsilon}(u(t)),\qquad 0\le t\le1,
\]

on a Banach space. The phase argument of `f` is periodic with fixed period (2\pi); `ε` scales the phase rate. The solution is assumed to exist and remain in a bounded set (K) for all (0<\varepsilon\le1). The vector field is (p)-times continuously differentiable in the phase, and it plus its phase derivatives has a uniformly bounded analytic extension in the state variable to a complex neighborhood of (K). These hypotheses give uniform bounds needed for the averaging transform and its derivatives (Introduction and Theorem 2.1).

The important structural point for this audit is that the fast phase is a *fixed periodic carrier* (t/\varepsilon), while the transformed slow vector field and all required derivatives stay uniformly bounded. The theory is not stated for an arbitrary variable instantaneous frequency, a real-to-imaginary eigenvalue transition, or a piecewise-smooth event coefficient.

## Derivative-free transform and error statements

The paper builds a near-identity periodic map (\Phi_\theta^{[n]}) and averaged vector field (F^{[n]}) recursively. The derivative-based recurrence uses (\partial_u\Phi^{[k]}F^{[k]}); the contribution here is to replace that directional derivative with

\[
D_\eta g(u)v=\frac{g(u+\eta v)-g(u)}{\eta},\qquad \eta=\varepsilon^k,
\]

inside the same iteration. The implementation maintains (\Phi=\mathrm{id}+\varepsilon B), evaluates the original vector field at shifted states, averages over phase, and subtracts the phase mean to preserve periodicity (Eqs. (2.1), (2.4), (2.8)–(2.9), Table 1). It avoids symbolic or explicit Jacobian/tensor derivatives, but it does not eliminate the repeated vector-field evaluations or phase integration.

Theorem 2.1 first gives the defect estimate for the exact derivative-based iteration. In addition to the bounded-solution and state analyticity conditions, (f_\theta) must have the stated phase smoothness and bounded phase derivatives. For a chosen iteration order (n), the theorem requires ((n+1)\varepsilon\le\varepsilon_0); the defect and its phase derivatives through order (p) are bounded by a geometric-in-(n) expression of the form \((Q_p(n+1)\varepsilon/\varepsilon_0)^n\). Theorem 2.2 proves that, for fixed (n), the finite-difference recurrence retains an (O(\varepsilon^n)) defect, with constants depending on the order and the smoothness bounds. Its proof first bounds the derivative-free and derivative-based maps’ difference by (O(\varepsilon^{k+1})) at iterate (k), then transfers the exact-iteration defect estimate.

For the two-scale formulation, the prepared initial datum is

\[
U(0,\theta)=\Phi_\theta^{[n]}(\tilde u_0),\qquad
\tilde u_0=(\Phi_0^{[n]})^{-1}(u_0).
\]

Theorem 3.1 states uniform bounds on mixed time and phase derivatives of (U) up to the specified orders. The Fourier phase discretization followed by Duhamel integration uses exact exponential oscillatory weights (p_{\ell,k}^{[n]}) and an order-(n) interpolation of the forcing history; the local truncation error is (O(\Delta t^{n+2})), giving order (n+1) with an order-(n) preparation. Starting values are supplied by a lower-order method. The two-scale implementation costs (O(N_\theta\log N_\theta)) per time step for the phase FFT and stores/evaluates a phase grid (Eqs. (3.3)–(3.7)).

For the micro–macro formulation, (u=\Phi_{t/\varepsilon}^{[n]}(v)+w), with slow macro equation \(\dot v=F^{[n]}(v)\) and a micro equation driven by the transform defect. Theorem 4.1 gives order (r=\min(p,q,n)) with an error constant independent of (\varepsilon), assuming the numerical method has order (q) on the nonstiff (\varepsilon=1) problem and the stated small-step condition. The derivative-free Adams–Bashforth micro–macro implementation also has per-step (O(N_\theta\log N_\theta)) cost. The authors explicitly report it is less efficient than two-scale AB because it recomputes the transform at every step; AB-TSF computes it once, or a few times when restarted (Eqs. (4.1)–(4.8), discussion after Eq. (4.8)).

## What the numerical section establishes

The experiments are three distinct periodic/near-periodic applications, not a varying-frequency oscillator benchmark:

- **Hénon–Heiles:** after filtering its stiff harmonic component, the system is put in the periodic (f_{t/\varepsilon}) form. Third- and fourth-order AB-TSF and AB-MM are tested with (N_\theta=128) for (\varepsilon\) from order one down to (10^{-4}). The comparison solution is a fourth-order micro–macro calculation with (\Delta t=10^{-5}), (N_\theta=256). The plots show the expected third/fourth order slopes across sampled (\varepsilon), with phase-discretization error decreasing spectrally as (N_\theta) grows. A restarted AB-TSF experiment plots energy error over a long interval; the paper notes linear drift is expected because the method is neither symmetric nor symplectic.
- **Nonlinear Klein–Gordon in the nonrelativistic scaling:** a filtered complex field produces the periodic high-phase form; the spatial operator is spectrally discretized and phase uses (N_\theta=128). Errors in (H^1) are plotted against timestep and (\varepsilon), for third- and fourth-order AB-TSF, with corresponding AB-MM plots described as similar and omitted. The phase-grid error is separately shown against (N_\theta). This demonstrates the proposed construction on a PDE after spatial discretization; it does not validate a turning-point crossing.
- **Vlasov with a strong magnetic field:** the stated regime has (B=B_0+O(\varepsilon)), with constant leading-order magnitude, so the stiff characteristic flow is periodic. A single characteristic is used to illustrate AB-MM at (\varepsilon=2^{-5}), (\Delta t=1/4), then quadratic interpolation recovers a fine trajectory. This is an illustration, not a parameter-uniform convergence study for vanishing frequency.

The paper gives no tensor-mode benchmark, event-crossing experiment, final-spectrum error, self-consistent outer-root error, or end-to-end cost against a batched transfer kernel. Its claimed uniformity is supported by the theory’s stated model and assumptions, while its plots show sampled convergence behavior for those examples.

## Mapping to Fast and the horizon turn

Fast’s public/reference equations evolve (z=\ln(k/aH)), (x), and (y):

\[
z'=1.5\sigma-1,\quad
x'=(-3+1.5\sigma)x-e^z y,\quad
y'=(-1+1.5\sigma)y+e^z x.
\]

The production kernel factors the background and advances a scaled two-state mode with frozen-(z) matrix

\[
A(q)=\begin{pmatrix}-1&-q\\q&1\end{pmatrix},\qquad q=e^z=k/(aH),
\]

whose eigenvalues are (\pm\sqrt{1-q^2}). Thus (q<1) is a real-eigenvalue, nonoscillatory region; at (q=1), the eigenvalues coalesce at zero and the frozen matrix is nilpotent; for (q>1), the pair is (\pm i\sqrt{q^2-1}). This is the mode’s horizon turning structure. Fast’s exact reheating event is also a kink in (\sigma(N)), and the production path treats it as an exact segment boundary.

**Conclusion from the model comparison:** the paper’s UA theorem does not apply uniformly across Fast’s full event-to-present chart. The natural oscillatory phase on the (q>1) side is generated by (\sqrt{q^2-1}), which vanishes at the turn. Reparameterizing time by its accumulated phase makes the inverse phase rate singular as (q\downarrow1); beyond the turn the generator has real rather than purely imaginary eigenvalues. There is therefore no single nonvanishing periodic carrier (t/\varepsilon) with uniformly bounded transformed derivatives across both sides. At the (q=1) endpoint the spectral collision also prevents a uniformly conditioned diagonalization. This is a direct mismatch with the fixed-periodic-phase/regularity hypotheses, not a claim that the derivative-free iteration itself fails for systems that satisfy those hypotheses.

Segmenting the path could expose a narrower use: on a smooth subhorizon interval (q\ge1+\delta), define a mode-dependent phase from the oscillatory eigenvalue and try a time-dependent extension of the two-scale method. That would require a new proof for nonuniform phase speed, bounds on derivatives of the phase and transformed field independent of (\delta), exact matching at the selected boundary, and a separate solution on the evanescent side. Those constants deteriorate as (\delta\to0), so it does not remove the turning layer. The true reheating kink would require another exact restart and compatible preparation. Low-frequency present-horizon modes may have little or no substantial oscillatory interval, where phase averaging has little to average.

The derivative-free feature is also poorly matched to this linear two-dimensional kernel. Fast already evaluates a closed-form (2\times2) frozen-frequency transfer with a bounded phase increment and uses Bessel/WKB machinery in the appropriate regimes. A two-scale scheme would add a phase grid and Fourier/interpolation work per mode (or require a new shared transform across frequencies); the micro–macro variant would recompute a transform each step. The paper’s (O(N_\theta\log N_\theta)) cost does not establish a speed advantage over the current batched Numba path. Any oscillatory-segment experiment would need measure total cost, including phase preparation, boundary matching, full sampled history, spectral assembly, and the self-consistent map.

## Difference from neighboring methods already assessed

- **Earlier periodic UA / averaging schemes:** this article keeps the two-scale and micro–macro architecture but makes the near-identity averaging transform systematic and derivative-free by replacing its directional Jacobian action with finite differences. That is a substantial implementation simplification for nonlinear vector fields with difficult Jacobians. It does not change the fixed periodic carrier requirement or introduce a turning-point treatment. This is the distinction from Chartier–Lemou–Méhats–Vilmart, “A new class of uniformly accurate numerical schemes for highly oscillatory evolution equations,” already read in [the unified framework screen](unified_fast_framework_screen_20261106.md).
- **Gautschi/trigonometric or exponential integrators:** those methods exploit an explicitly identified linear oscillatory operator and integrate its oscillatory flow (or its φ-function weights) directly, then approximate the remainder. The present paper instead averages a periodic vector field in an added phase variable and handles nonlinear vector-field evaluations through a phase Fourier grid. Fast’s two-state transfer is already close to a specialized linear propagator. The paper supplies no new scalar (\omega(N,k)) turning-point formula or phase-integral error bound to replace that transfer.
- **ARDC/Riccati–Chebyshev:** ARDC converts the damped scalar mode equation to a log-derivative/Riccati problem and switches between oscillatory and nonoscillatory numerical charts using local defect control; it is designed to face a turning point through regime switching. This paper does not evolve a Riccati variable or estimate a local transfer defect. The two methods are algorithmically distinct; the derivative-free UA paper does not validate ARDC nor duplicate its turn-switching mechanism. The existing source-level ARDC audit is [riccati_source_audit_20261105.md](riccati_source_audit_20261105.md).

## Decision and smallest useful falsifier

Do not adopt this method as a unified Fast replacement on the present evidence. The full-chart periodic-phase hypotheses fail at (q=1) and on the superhorizon real-spectrum side, and exact reheating adds a nonsmooth coefficient boundary. The main derivative-free advantage is not compelling for Fast’s already linear, two-state mode kernel.

If a future, restricted subhorizon experiment is proposed, first require a written decomposition on (q\ge1+\delta) that identifies the periodic carrier and proves uniform bounds for the phase-transformed vector field and its derivatives. Then apply its preparation and one AB-TSF step to the existing two cases with (q_\mathrm{today}=2,4), comparing full mode histories and assembled channels against a high-accuracy direct solve at identical output nodes. Vary (\delta\) toward the turn: growth of the transform/derivative bounds, phase-node count, error, or runtime before reaching the full-frequency map is a concrete rejection gate. This is a proposed falsifier only; no prototype or numerical experiment was run for this paper.

## References checked

1. Chartier, P., Lemou, M., Méhats, F., & Zhao, X. (2022). [Derivative-free high-order uniformly accurate schemes for highly-oscillatory systems](https://doi.org/10.1093/imanum/drab014). *IMA Journal of Numerical Analysis*, 42(2), 1623–1644. [Official OUP record](https://academic.oup.com/imajna/article-abstract/42/2/1623/6260882); [full author-manuscript PDF inspected](https://scispace.com/pdf/derivative-free-high-order-uniformly-accurate-schemes-for-14utse9nfe.pdf); [HAL record](https://hal.inria.fr/hal-03141156).
2. Chartier, P., Lemou, M., Méhats, F., & Vilmart, G. (2020). [A new class of uniformly accurate numerical schemes for highly oscillatory evolution equations](https://doi.org/10.1007/s10208-019-09413-3). *Foundations of Computational Mathematics*, 20, 1–33.
3. Project context: [tensor-mode equation and event treatment](../physics.md), [unified Fast framework source screen](unified_fast_framework_screen_20261106.md), and [ARDC source audit](riccati_source_audit_20261105.md).
