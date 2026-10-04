# FragReco

## Overview

Implementation of the heavy quark hadronization framework for the vacuum and the quark-gluon plasma.
It provides vacuum fragmentation and a hybrid QGP treatment that combines
fragmentation with recombination (coalescence, implemented here as iSRM). 

In the fragmentation path, heavy quarks are converted into hadrons with an
explicit Peterson momentum-fragmentation function. Hadron-species probabilities
are sampled separately from a statistical hadronization model. The associated reference is
[Charm-hadron production in pp and AA collisions](https://arxiv.org/abs/2002.00392).

In the hybrid QGP path, recombination is attempted first, and heavy quarks that
do not recombine follow the fragmentation path. The corresponding model and
hadron-chemistry context are described in
[Charmed hadron chemistry in relativistic heavy-ion collisions](https://arxiv.org/abs/1911.00456).

## Build

FragReco requires CMake 3.16 or newer, a C++17 compiler, and Boost
Program_options.

```bash
mkdir build && cmake ..
make 
```

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

For fragmentation-only mode 1, use the same command and change `--mode 3`
to `--mode 1`, `--tchem 0.160' to `--tchem 0.170 '. The Recomb/Wigner and omega options may remain in the command;
they are accepted but are not used in mode 1.

These values currently match the charm production defaults.



## Hadronization modes

- `1 = Frag`: Peterson fragmentation and charm feeddown.
- `2 = Recomb`: currently a pass-through mode. HQs are written unchanged.
- `3 = Frag + Recomb`: the normal production mode. Recombination is attempted
  first, and unsuccessful candidates follow the fragmentation path.

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
  --omega-m 0.24 \
  --omega-b 0.24 \
  --recomb-table data/recomb_c_raw.dat \
  --wigner-table data/max_wigner_c.dat
```

## Hadron origin flag

Each particle carries a hadronization-origin tag internally:

```text
0  Unknown
1  Frag
2  Recomb
```

The tag is written only by the `analysis` format. Charm feeddown preserves the
Frag or Recomb origin of the primary hadron.


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
