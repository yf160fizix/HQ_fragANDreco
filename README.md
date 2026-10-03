# FragReco

## Overview

FragReco is the current implementation of the `HQ-fragANDreco` heavy-quark
(HQ) hadronization framework for vacuum and the quark--gluon plasma (QGP). It
provides vacuum fragmentation and a hybrid QGP treatment that combines
fragmentation with recombination (coalescence, implemented here as iSRM). The
current production setup is charm focused and includes charm-hadron feeddown.

In the fragmentation path, heavy quarks are converted into hadrons with an
explicit Peterson momentum-fragmentation function. Hadron-species probabilities
are sampled separately from thermal/statistical charm-chemistry weights. The
original `HQ-fragANDreco` README described this path more broadly as
"HQET-inspired fragmentation"; the associated model reference is
[Charm-hadron production in pp and AA collisions](https://arxiv.org/abs/2002.00392),
but this implementation should not be read as including every ingredient of
that work.

In the hybrid QGP path, recombination is attempted first and heavy quarks that
do not recombine follow the fragmentation path. The corresponding model and
hadron-chemistry context are described in
[Charmed hadron chemistry in relativistic heavy-ion collisions](https://arxiv.org/abs/1911.00456).

## Build

FragReco requires CMake 3.16 or newer, a C++17 compiler, and Boost
Program_options.

```bash
cmake -S . -B build -DISRM_STRICT_WARNINGS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

The chemistry calculation uses `std::cyl_bessel_k` when available and a
Boost.Math fallback otherwise.

## Basic usage

Run from the package root so the default table paths resolve correctly:

```bash
./build/hq_hadronizer_cpp --help
```

For a charm production run in mode 3, specify the complete physics
configuration and all input tables explicitly:

```bash
./build/hq_hadronizer_cpp \
  input.oscar \
  output.oscar \
  --mode 3 \
  --io-format urqmd \
  --hq 4 \
  --recomb-table data/recomb_c_raw_M020_B0267.dat \
  --wigner-table data/max_wigner_c_M020_B0267.dat \
  --meson-table data/charm_meson_states.csv \
  --baryon-table data/charm_baryon_states.csv \
  --epsilon-meson 0.01 \
  --epsilon-baryon 0.03 \
  --gamma-s 0.7 \
  --gamma-hb 1.0 \
  --tchem 0.160 \
  --omega-m 0.20 \
  --omega-b 0.267
```

For fragmentation-only mode 1, use the same command and change only `--mode 3`
to `--mode 1`. The Recomb/Wigner and omega options may remain in the command;
they are accepted but are not used in mode 1.

These values currently match the charm production defaults, but spelling them
out records the complete physics configuration.

`--seed` is optional. If omitted, a random seed is generated automatically and
printed in the run summary. Specify `--seed <value>` when exact reproducibility
is required.

## Hadronization modes

- `1 = Frag`: Peterson fragmentation and charm feeddown; Recomb/Wigner tables
  are not loaded.
- `2 = Recomb`: currently a pass-through mode. HQs are written unchanged, and
  neither Frag chemistry nor Recomb/Wigner tables are loaded.
- `3 = Frag + Recomb`: the normal production mode. Recombination is attempted
  first and unsuccessful candidates follow the fragmentation path.

## I/O formats

`--io-format` explicitly selects the input and output schema; formats are not
auto-detected. The default is `urqmd`.

```text
standard: ... ipx ipy ipz wt
urqmd:    ... ipx ipy ipz iE wt 0
analysis: ... ipx ipy ipz wt origin
```

`iE` and the final zero are UrQMD compatibility fields and appear only in the
UrQMD format. Standard output is unchanged by analysis metadata.

To write the per-particle hadronization origin, use the same command with
`--io-format analysis`.

The final analysis column is the hadronization origin:

```text
origin = 0  Unknown
origin = 1  Frag
origin = 2  Recomb
```

## Main physics parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `Tchem` | 0.160 GeV | charm chemistry/hadronization temperature |
| `gamma_s` | 0.7 | strange-hadron chemistry factor |
| `gamma_HB` | 1.0 | primary heavy-baryon fragmentation weight factor |
| `eps_M` | 0.01 | Peterson epsilon for charm mesons |
| `eps_B` | 0.03 | Peterson epsilon for charm baryons |
| `omega_M` | 0.20 GeV | charm-meson oscillator scale |
| `omega_B` | 0.267 GeV | charm-baryon oscillator scale |
| `m_c` | 1.8 GeV | charm-quark mass |
| `m_b` | 5.2 GeV | bottom-quark mass |
| `m_q` | 0.30 GeV | light-quark mass |
| `m_s` | 0.40 GeV | strange-quark mass |
| `m_g` | 0.30 GeV | effective gluon mass |

The corresponding CLI options are shown by `--help`. Constituent masses are
configuration defaults rather than separate CLI options.

## Default charm inputs

The production defaults are:

```text
data/recomb_c_raw_M020_B0267.dat
data/max_wigner_c_M020_B0267.dat
data/charm_meson_states.csv
data/charm_baryon_states.csv
```

The Recomb and Wigner tables correspond to `omega_M = 0.20 GeV` and
`omega_B = 0.267 GeV`. If either omega is changed in mode 3, both
`--recomb-table` and `--wigner-table` must be supplied explicitly. FragReco
does not infer table parameters from filenames.

For example, to use the available alternate omega/table pair, replace all four
matching arguments in a complete mode-3 command:

```bash
  --omega-m 0.21 \
  --omega-b 0.259 \
  --recomb-table data/recomb_c_raw_M021_B0259.dat \
  --wigner-table data/max_wigner_c_M021_B0259.dat
```

## Hadron origin metadata

Each particle carries a hadronization-origin tag internally:

```text
0  Unknown
1  Frag
2  Recomb
```

The tag is written only by the `analysis` format. Charm feeddown preserves the
Frag or Recomb origin of the primary hadron.

## Validation and tests

The CTest suite covers Lorentz transformations, chemistry and Peterson
parameters, feeddown and Recomb behavior, mode handling, the three I/O schemas,
and origin metadata.

The fixed-input, fixed-seed production regression has the expected output
SHA-256:

```text
a29a85e3460e968ce587b27eba87c630a527998ef9d7011558345505162d05c0
```

This regression checks that the validated default simulation output remains
unchanged for the fixed input and RNG seed.

## Current limitations

- The `Xi_c`/`Omega_c` Recomb channel groups corresponding to channels 19--24
  are not implemented. If one of these channels is selected, the HQ is omitted
  from the output and counted as an unsupported-channel drop.
- Mode 2 is accepted as the Recomb-only mode, but its current top-level path
  passes HQs through unchanged rather than forcing recombination.
- Recombination uses externally tabulated Wigner maxima as rejection envelopes.
  The runtime code assumes that these envelopes bound the evaluated Wigner
  weights; it does not independently derive or validate the maxima.
- Bottom support is incomplete. The fragmentation chemistry is charm-only, and
  the bottom Recomb channel/table coverage is not a complete production model.
- Mixed charm/bottom input is not a supported production workflow.
