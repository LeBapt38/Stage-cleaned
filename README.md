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
