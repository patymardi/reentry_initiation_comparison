# Reentry Initiation Comparison

This repository contains the code to compare atrial reentry initiation using different pacing protocols in patient-specific atrial computer models.

This project is a collaboration between the University of Bordeaux, Karlsruhe Institute of Technology (KIT), University of Washington, Queen Mary University of London, Medical University of Graz, and Pontificia Universidad Católica de Chile.


## Data

The repository includes two patient-specific biatrial bilayer meshes derived from MRI, prepared for simulations with **openCARP**. 

## Repository structure

```
.
├── src/
│   ├── pacing_protocols/
│   │   ├── cross_field.py
│   │   └── other_protocol.py
│   ├── filament/
│   │   └── run_igbfilament.py
│   └── requirements.txt
│
├── data/
│   ├── meshes/
│   │   ├── patient12/
│   │   └── patient51/
│   ├── parameters.par
│   ├── element_tag.csv
│   └── link.txt
│
│
└── README.md
```

## Description

### `src/`

Source code used for simulation and post-processing.

- `pacing_protocols/` — implementation of the different pacing protocols used in the study.
- `filament/` — scripts for filament extraction and analysis.
- `requirements.txt` — Python dependencies.

### `data/`

Each patient folder contains:

- Bilayer meshes (`.pts`, `.elem`, `.lon`)
- Stimulation points
- Fibrosis distributions (`.regele`)

Additional files:

- `parameters.par` — ionic model scaling factors.


## Installation

Create a virtual environment `.venv` and install the required Python packages:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r src/requirements.txt
```

Once the virtual environment `.venv` is active, install carputils from the local source:

```bash
python -m pip install --upgrade pip setuptools
python -m pip install -e /path/to/openCARP/external/carputils
```
Replace `/path/to/openCARP/external/carputils` with the path to your local carputils directory.

Verify the installation:
```bash
python --version
python -c "import carputils; print(carputils.__file__)"
python -m pip show carputils
python -m pip check
```

## Dependencies

- Python 3.13
- openCARP
- carputils

## Contributing

Contributions are welcome :)
