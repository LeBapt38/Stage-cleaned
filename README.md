# M2 intership at LAPTh

### Table of content

- [Froissart-Gribov Projection](#froissart-gribov-projection)
    - [Project Goal (Froissart-Gribov)](#project-goal)
    - [Files (Froissart-Gribov)](#files)
        - [Code_report.nb](#code_reportnb)
        - [Test_building_FG.nb](#test_building_fgnb)
        - [FG_disc_analytic.wl](#fg_disc_analyticwl)
        - [Regge_traj.nb](#regge_trajnb)
- [The Regge residues of the Veneziano amplitude](#the-regge-residues-of-the-veneziano-amplitude)
    - [Project goal (Veneziano)](#project-goal-1)
        - [1. Stokes Lines and Correspondence with Roots](#1-stokes-lines-and-correspondence-with-roots)
        - [2. Migrating Saddle and Subleading Orders in k](#2-migrating-saddle-and-subleading-orders-in-k)
        - [3. Asymptotic Formulas and Fine-Structure Matching](#3-asymptotic-formulas-and-fine-structure-matching)
        - [4. Analytical Matching and Topological Structures](#4-analytical-matching-and-topological-structures)
        - [5. Future Directions (Open Research Horizons)](#5-future-directions-open-research-horizons)
    - [Files (Veneziano)](#files-1)
        - [Toymodel_Bessel.nb](#toymodel_besselnb)
        - [Pol_Veneziano.nb](#pol_venezianonb)


# Froissart-Gribov Projection

## Project Goal

This part of the project focuses on the Froissart-Gribov projection of the ansatz from acking S-matrix bounds across dimensions. Details on the physics can be found in the report. Important formulas are reminded below. 

The goal is to turn the crossed-channel representation of the amplitude into a numerical J-plane object that can be scanned and later used for the Regge analysis.

The starting point is the crossing-symmetric rho ansatz,

$$
\rho(s,\sigma)=\frac{\sqrt{\sigma-4}-\sqrt{4-s}}{\sqrt{\sigma-4}+\sqrt{4-s}},
$$

with the corresponding discontinuity

$$
\mathrm{disc}\,\rho(s,\sigma)=\frac{2\sqrt{s-4}\sqrt{\sigma-4}}{s+\sigma-8}.
$$

The amplitude is built as a symmetric sum over the three Mandelstam variables,

$$
\mathcal M(s,t,u)=\frac16N_d\sum_{i,j,n,m} \alpha_{ij}^{nm}\left[\,\mathcal 
\rho(s,\sigma_i)^n\rho(,\sigma_j)^m
\;+\mathrm{crossed \ sym}\right]+
\;\text{threshold terms},
$$

with $s+t+u=4$ and with extra leading and subleading threshold pieces used to control the behavior near $s=4$.

The Froissart-Gribov projection is evaluated with the partial-wave kernel $P_J^{(d)}$ and the discontinuity across the relevant cut,

$$
\mathcal Q_J(s)=\int dz\, (z^2-1)^{\frac{d-4}2}Q_J^{(d)}(z)\,\mathrm{Disc_t}\,\mathcal M\big(s,t(s,z),u(s,z)\big),
$$

where the change of variables used in the notebooks is

$$
t(s,z)=\frac{4-s}{2}(1-z),\qquad u(s,z)=\frac{4-s}{2}(1+z).
$$

In practice, the work moves from the report version to a more robust analytic implementation of the discontinuity, then to a server-ready script that computes a J-plane grid for the final analysis.

## Files

### Code_report.nb

Notebook used to generate the figures that appear in the report.

Main differences with recent code:
- `iQInfinity12[start, J, d]`: analytic tail integral used for the large-z asymptotic correction.
- `tail[amplitudeT, s, J, d, start]`: estimates the asymptotic tail by sampling the t-channel amplitude at large `z` and matching it to the analytic tail.
- `exprJmain[J, s]`, `exprJtail[J, s]`, `exprJ[J, s]`: report-era evaluation of the projected quantity.


How to use it:

- Use it as the notebook that reproduces the report plots.
- Old implementation so less precise.

### Test_building_FG.nb

Notebook used to test the different pieces of the ansatz and build the last, more robust version of the projection. This version uses the analytic discontinuity and improves the report implementation.

Main entry points:
- `PolQ[J, d, z]`: Froissart-Gribov kernel used in the projection integral. `J` is the spin, `d` the spacetime dimension, and `z` the integration variable.
- `PolQDer[J, d, z]`: first derivative of `PolQ` with respect to `z` when a higher-order kernel expansion is needed.
- `z0t[s]`, `z0u[s]`: location of the t- and u-channel branch points.
- `tt[s, z]`, `uu[s, z]`: maps from the variable `z` back to Mandelstam `t` and `u`.
- `Pgen[J, d, z]`, `IPt[J, s]`: comparison with the standard partial-wave projection based on the Gegenbauer kernel.
- `MyNIntegrateFaster[expr, {x, a, b}, p]`: numerical integration helper with controlled precision. `expr` is the integrand, `{x, a, b}` the integration domain, and `p` the target precision.
- `HalfCompactify[s0, \[Phi]]`, `CompactIntegrand[x, var, s0]`: compactification helpers used to map semi-infinite integration domains to a finite interval.
- `discZ[expr, s, z, eps]`: finite-difference discontinuity check used for validation.
- `exprJmain[J, s]`, `exprJthresh[J, s]`, `exprJ[J, s]`: the tested projection pipeline, split into the regular integral and the threshold branch point contribution.
- `Pgen[J, d, z]`, `IPt[J, s]`: standard projection used for comparisons against the Froissart-Gribov result.

How to use it:

- Use it as the validation notebook for the ansatz pieces, discontinuities, and projection formulas.
- The notebook is mainly for cross-checks between the analytic Froissart-Gribov result and the standard partial-wave projection.
- It is the right place to compare the original report version against the improved analytic discontinuity version.

### FG_disc_analytic.wl

Final analytic implementation of the Froissart-Gribov projection. This is the end product for the expression and contains the code that computes the J-plane grid. It is meant to run on a server (just change first line). The current implementation works only for d even...

Main entry points:

- `\rho[s, \sigma]`, `disc\rho[s, \sigma]`: conformal map and its discontinuity.
- `\rhoansatz[s, t, u, n, m, \sigmai, \sigmaj]`: crossing-symmetric ansatz built from rho blocks. `n` and `m` are the rho powers, while `\sigmai` and `\sigmaj` are the grid points in the sigma plane.
- `discT\rhoansatz[...]`, `discU\rhoansatz[...]`: analytic discontinuities of the ansatz in the t and u channels. Works only for `n,m=1`.
- `Mraw[...]`, `discTMraw[...]`, `discUMraw[...]`: raw ansatz tables before threshold terms are added.
- `M[d, s, t, u, \sigmag, pg]`, `discTM[...]`, `discUM[...]`: full amplitude and discontinuity tables, including the threshold pieces for M. For the discontinuities, only the regular term is taken, the rest will be added manualy in `exprJkeyhole[J, s, n]`, `exprJthresh[J, s]` already computed.
- `myM[d, Nmax]`, `myDiscTM[d, Nmax]`, `myDiscUM[d, Nmax]`: convenience wrappers that build the amplitude and discontinuities using a chosen dimension `d` and a grid size `Nmax`.
- `exprJmain1[J, s]`, `exprJmain2[J, s]`: regular Froissart-Gribov integrals for the t- and u-channel cuts of the regular part.
- `exprJkeyhole[J, s, n]`, `exprJthresh[J, s]`: handling of the threshold singular pieces.
- `exprJ[J, s]`: full analytic projection.
- `Pgen[J, d, z]`, `IPt[J, s]`: comparison with the standard projection kernel.
- `adaptiveGrid[exprJNum, s, largeRes, smallRes, xRange, yRange, coarseIn]`: adaptive meshing routine used to compute the J-plane grid efficiently.

How to use it on a server:

1. Put the code and its dependency files in a local folder, for example `Exchange`.
2. Copy the folder to the server:

```bash
scp -r Exchange guillem@lapthui.in2p3.fr:~
```

3. Connect to the server and move to the job folder:

```bash
ssh guillem@lapthui.in2p3.fr
cd ~/Exchange
ls
```

4. Check that `wolframscript` is available:

```bash
wolframscript -version
```

If it is not known, add the Wolfram installation to `~/.bashrc` and reconnect so that `wolframscript` is available in the login shell.

5. Launch the computation. The script takes three arguments in order: the coefficient-amplitude input file, the dimension `d`, and the output `.mx` file.

```bash
nohup wolframscript -file FG_disc_analytic.wl coeffAmpPath d outPath \
> log.txt 2>&1 &
```


6. Check that it is running:

```bash
tail -f log.txt
```

Ctrl+C only exits `tail`; it does not stop the computation. You can also leave the server without stopping the job.

7. Retrieve the output from the laptop once the run is finished:

```bash
scp guillem@lapthui.in2p3.fr:'~/Exchange/*.mx' .
```

### Regge_traj.nb

Notebook used to build the Regge trajectories.

Main entry points:

- `minPtsFromGraph[data_]`: scans a sampled grid and extracts candidate minima to initialize the secant method. `data` is the list of `{x, y, z}` points.
- `secantOrder2[myFunc, start, nmax, tol]`: second-order secant solver used to refine a root from an initial guess. `myFunc` is the function to solve, `start` is the initial point, `nmax` is the maximum number of iterations, and `tol` is the target tolerance.
- `trackZero[myFunc, start, smin, smax, ds]`: follows a zero of `myFunc` by repeatedly calling `secantOrder2` while scanning `s` from `smin` to `smax` in steps of `ds`. `start` is the initial root used for the first step.
- `cleanTraj[traj_]`: removes duplicate or failed points from a trajectory list before plotting or fitting.

How to use it:

- This notebook is the trajectory-building step that comes after the Froissart-Gribov grid is available.
- It still requires manual input to choose the starting points for the trajectories.

# The Regge residues of the Veneziano amplitude

## Project goal

This part of the project studies the Regge residues of the Veneziano amplitude as a complex-analytic object, with the aim of understanding how the zeros of the residue polynomials organize in the complex plane and what this organization reveals about asymptotic dynamics.

In the report-level perspective, the objective is to characterize the root locus of the residue polynomials $p_k(s)$, to connect this locus with asymptotic/saddle-point control at large $k$, and to identify robust geometric structures that persist across orders. The practical workflow starts from numerical root tracking and visualization (domain coloring, phase diagnostics), then moves to asymptotic interpretation, and finally to analytical matching formulas.

The recent work extends this baseline substantially and gives a global interpretation in terms of Stokes geometry, migrating saddles, and multi-sheet topology.

### 1. Stokes Lines and Correspondence with Roots

The discrete roots of $p_k(s)$ are controlled by the Stokes phenomenon of the large-$k$ asymptotic expansion.

- Saddle interferences: roots appear where two competing saddle contributions have equal magnitude and opposite phase, yielding destructive interference:

$$
\beta_k(s) \sim c_A e^{k S_A(z)} + c_B e^{k S_B(z)} = 0
\quad \Longrightarrow \quad
\operatorname{Re}(S_A)=\operatorname{Re}(S_B).
$$

- C-shell caustic: in the thermodynamic limit $k\to\infty$, discrete roots condense onto anti-Stokes lines. The observed C-shell (cardioid-like envelope) is interpreted as the geometric caustic of this phase-transition boundary. Inside the shell, roots remain dense and real; outside, they move into complex branches as one saddle becomes exponentially dominant.

### 2. Migrating Saddle and Subleading Orders in $k$

A purely leading-order analysis in $k$ is insufficient to recover the C-shell boundary. The decisive mechanism appears at subleading order.

- $1/k$ mechanism: subleading corrections inject an effective migrating saddle into the complex $t$-plane action.
- Saddle interaction: this migrating saddle is not static. As $z$ changes, it travels across the $t$-plane, breaks translation symmetry, and collides/interacts with a fixed periodic tower of saddles.

This interaction is the key ingredient needed to reproduce the global root geometry rather than only local asymptotic behavior.

### 3. Asymptotic Formulas and Fine-Structure Matching

Restricting to the complex $t$-plane, $t=x+i\chi$, one obtains a structured transcendental picture.

- Approxiamte position of the migrating tower:

$$t_{\text{mig}} \approx -\frac{z \sqrt{k}}{\sqrt{1 + z + \frac{1}{k}}}$$

- Tower asymptotics: tower saddles are asymptotically separated by $2\pi$ along the imaginary direction, $\chi\approx 2\pi n$. On high sheets, their real part drifts logarithmically:

$$
x\approx 2\ln\left|2\pi n\left(1+\frac{1}{z}\right)\right|.
$$

- Fine-structure matching equation: matching tower phase to the dominant background gives a transcendental relation for root layers/rays:

$$
\arg(z)-\frac{1}{2}\arg\left(1+z+\frac{1}{k}\right)
=-\frac{\pi}{2}-\frac{\ln\left|2\pi n\left(1+\frac{1}{z}\right)\right|}{\pi n}.
$$

This equation tracks the successive transitions of the migrating root across discrete branches.

### 4. Analytical Matching and Topological Structures

Numerical domain-coloring maps and stabilized phase-field plots support a genuinely multi-sheeted analytic structure.

- Multi-sheeted Riemann surface: logarithmic and hyperbolic branch cuts partition the $t$-plane into an infinite stack of sheets.
- Moving topological defect: the migrating saddle behaves as a moving defect crossing branch cuts as $z$ evolves, passing from the principal sheet to higher sheets and reorganizing the relevant Stokes graph.

### 5. Future Directions (Open Research Horizons)

These points are prospective research directions and not yet final results.

- Lefschetz-thimble bifurcation analysis: decompose the physical contour as
$
\mathcal{C}=\sum_{\sigma} n_{\sigma}\,\mathcal{J}_{\sigma}
$
and test whether C-shell crossing induces jumps in intersection numbers $n_\sigma$, making formerly dominant saddles topologically inactive.

- Geometric uniformization: construct a uniformizing map that unfolds the multi-sheeted surface into a single periodic domain, to avoid explicit saddle tracking across branch cuts.


- Physical interpretation (UV completion and resurgence): the migrating saddle could be viewed as an effective localized (particle-like) sector, while the periodic tower encodes non-local string-like excitations. A resurgence analysis of the $1/k$ sector could clarify how global tower data non-perturbatively restores high-energy consistency.

## Files

### Toymodel_Bessel.nb

Notebook containing two toy models where the Stokes description is developed analytically. It is used as the controlled setting to derive and test the saddle-interference picture before applying it to the full Veneziano residue problem.

How to use it:

- Use it as the analytic benchmark for the Stokes-line mechanism.
- Compare the exact and asymptotic behavior of the toy models with the qualitative structures observed later in the Veneziano analysis.

### Pol_Veneziano.nb

Main notebook for the Veneziano residue analysis. It explores the saddle structure in detail, then attempts a Stokes analysis of the corresponding phase geometry. Lastly, it builds the curve obtained from the transcendental expression of the argument.

Current status:

- The saddle-structure exploration is implemented and can be used directly.
- The Stokes-analysis plots are present but require careful graph-cleaning modifications to get fully robust visuals.
- The transcendental-argument curve construction is included as the final comparison layer.


