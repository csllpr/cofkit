# External Tools

External binaries are optional. They are needed only for Zeo++, LAMMPS, EQeq, gRASPA/RASPA2, and hybrid MD/MC workflows.

The CLI can read tool paths from the nearest `.env` file, unless the variable is already set in the shell.

## Zeo++

Use Zeo++ `0.3` from the upstream tarball:

- download: `http://www.zeoplusplus.org/zeo++-0.3.tar.gz`
- expected binary: `network`

Typical build:

```bash
curl -LO http://www.zeoplusplus.org/zeo++-0.3.tar.gz
gunzip zeo++-0.3.tar.gz
tar xvf zeo++-0.3.tar
cd zeo++-0.3/voro++/src
make
cd ../..
make
export COFKIT_ZEOPP_PATH="$PWD/network"
```

## LAMMPS

Install LAMMPS from `conda-forge`:

```bash
conda create -n cofkit-lammps -c conda-forge lammps
conda activate cofkit-lammps
export COFKIT_LMP_PATH="$(command -v lmp || command -v lmp_mpi)"
```

`cofkit` accepts either `lmp` or `lmp_mpi` as long as `COFKIT_LMP_PATH` points to a working executable.

## EQeq

Use the CIF-capable fork:

- source: `https://github.com/csllpr/EQeq`
- expected binary: `eqeq`

Typical build:

```bash
git clone https://github.com/csllpr/EQeq.git
cd EQeq
g++ main.cpp -O3 -o eqeq
export COFKIT_EQEQ_PATH="$PWD/eqeq"
```

That fork reads `data/ionizationdata.dat` and `data/chargecenters.dat` relative to the executable by default. `EQEQ_IONIZATION_DATA_PATH` and `EQEQ_CHARGE_CENTERS_PATH` remain available if you need to override them.

## gRASPA

Use one of these repositories:

- preferred: `https://github.com/csllpr/gRASPA`
- fallback upstream: `https://github.com/snurr-group/gRASPA`
- expected binary: `nvc_main.x`

The preferred fork carries local performance tweaks for higher-throughput GPU scheduling. gRASPA is source-first, so the exact build depends on your NVIDIA HPC SDK / CUDA installation. The checked-in `NVC_COMPILE` script in both repositories is the reference starting point.

Typical flow:

```bash
git clone https://github.com/csllpr/gRASPA.git
cd gRASPA
# Edit NVC_COMPILE if your NVIDIA HPC SDK / CUDA paths differ.
bash NVC_COMPILE
export COFKIT_GRASPA_PATH="$PWD/nvc_main.x"
```

If your local build places the binary somewhere else, such as `src_clean/nvc_main.x`, point `COFKIT_GRASPA_PATH` there instead.

## RASPA2

Use upstream RASPA2 when you need CPU-only Monte Carlo runs:

- source: `https://github.com/iRASPA/RASPA2`
- expected binary: `simulate`

Build RASPA2 with its upstream instructions, then point `cofkit` at the compiled simulator:

```bash
export COFKIT_RASPA2_PATH=/path/to/simulate
```

The existing `graspa-*` command names are retained for compatibility. Select RASPA2 with `--backend raspa2`; the default backend remains `graspa`.
