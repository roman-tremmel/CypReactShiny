# CypReact Shiny

A small **R Shiny** web app around [CypReact](https://bitbucket.org/Leon_Ti/cypreact)
(Tian *et al.*, *J. Chem. Inf. Model.* 2018). You only enter **SMILES** — a single
molecule or a whole batch — and get back, for each of the nine major human
cytochrome-P450 (CYP450) isoforms, whether the molecule is predicted to be a
**reactant** (metabolised) or a **non-reactant**.

Isoforms: `1A2, 2A6, 2B6, 2C8, 2C9, 2C19, 2D6, 2E1, 3A4`.

```
R = predicted reactant (metabolised by that CYP)
N = non-reactant
```

## Features

- **Single molecule** input with one-click examples (caffeine, aspirin, ...).
- **Batch** input: paste many SMILES (one per line / comma-separated) or upload
  a `.csv` / `.txt` / `.smi` file.
- Pick which CYP isoforms to run.
- Colour-coded results table and **CSV download**.
- Robust handling of edge cases (see *How it works*).

## Requirements

The host that serves the app needs:

1. **R** with the packages `shiny` and `DT`:
   ```r
   install.packages(c("shiny", "DT"))
   ```
2. **Java** — a JRE **8 or newer** must be on `PATH` (CypReact is a Java tool).
   Check with `java -version`. On Debian/Ubuntu: `apt-get install default-jre`.

> Python/R/Shiny on the server is **not enough on its own** — without Java the
> predictions cannot run.

That's it: the CypReact tool itself (`cypreact.jar` + model files) is vendored
in this repository under [`cypreact_bundle/`](cypreact_bundle/), so there is
nothing else to download.

## Run locally

From the repository root:

```r
shiny::runApp(".")
```

or

```bash
R -e 'shiny::runApp(".", port = 3838, host = "0.0.0.0")'
```

Then open http://localhost:3838.

## Deploy (Shiny Server / Posit Connect)

Push/copy the whole repository (keep the directory layout) to the server's apps
directory, e.g. `/srv/shiny-server/cypreact/`. Shiny Server picks up `app.R`
automatically. Make sure `java` is available to the user that runs the Shiny
process.

If you keep `cypreact.jar` somewhere else, point the app at it with an
environment variable instead of using the vendored copy:

```bash
export CYPREACT_BUNDLE=/path/to/CypReactBundle   # folder with cypreact.jar + supportfiles/
export JAVA=/usr/bin/java                          # optional, if java is not on PATH
```

## Layout

```
app.R                 # the Shiny app (UI + server)
R/cypreact.R          # wrapper that drives cypreact.jar and aligns results
cypreact_bundle/      # vendored CypReact tool
  cypreact.jar
  supportfiles/       # per-isoform model + support files
```

## How it works

`R/cypreact.R` shells out to `cypreact.jar`. The CLI has two awkward
behaviours that the wrapper hides so the result rows always line up 1:1 with
your input, in order:

1. **One bad SMILES aborts the whole batch.** The wrapper reads the offending
   SMILES from the Java error, marks it as `parse error`, drops it and re-runs
   — so a single typo no longer kills the run.
2. **Out-of-scope molecules are silently dropped** (mixtures, very large
   molecules, sugars/lipids, ...) and the remaining rows renumbered. The
   wrapper detects these and re-inserts them as `out of scope` (CypReact treats
   them as non-reactants).

## Credit & license

CypReact is by Siyang Tian *et al.* and distributed under the MIT License
(see `cypreact_bundle` origin: https://bitbucket.org/Leon_Ti/cypreact).
Please cite the CypReact paper when you use these predictions. This app is a
convenience front-end and does not change the underlying model.

> Predictions are *in-silico* estimates and are not a substitute for
> experimental data.
