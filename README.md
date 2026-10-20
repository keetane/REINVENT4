# Molecular Generator — REINVENT4 Web Application

This repository contains a local web application built on [REINVENT4](https://github.com/MolecularAI/REINVENT4). It provides one workspace for molecule generation, Transfer Learning, QSAR modeling, Reinforcement Learning, job history, and result analysis.

The application is called **Molecular Generator**. It is a local-only service: the browser UI runs on `http://127.0.0.1:5173/`, the FastAPI backend runs on `http://127.0.0.1:8000/`, and jobs continue in an independent worker after submission.

## What Molecular Generator provides

- **Molecule Generation**: sample from Standard Prior models or from registered Transfer Learning / Reinforcement Learning models. Reinvent does not require a Parent; Mol2Mol requires one, while LibInvent and LinkInvent use attachment-point SMILES.
- **Transfer Learning**: validate and standardize SMILES data, split Training and Validation data, and create a specialized Reinvent or Mol2Mol generator.
- **QSAR**: train a single model or compare algorithms such as Chemprop, QSARtuna, MLP, Random Forest, LightGBM, XGBoost, and CatBoost. Registered models can predict a drawn structure and can be selected as an RL Scoring Function.
- **Reinforcement Learning**: update an Agent toward a combined Scoring Function built from QSAR predictions, molecular properties, and other constraints.
- **Jobs & Results**: inspect progress, logs, training curves, generated structures, properties, model-comparison scores, and downloadable artifacts.

The Tutorial tab in the application explains the workflow in the same order as the tabs: environment setup, molecule generation, Transfer Learning, QSAR, Reinforcement Learning, and open-source components.

## Installation and setup

### Requirements

- Python 3.10 or newer
- Node.js and npm (Node.js 24 or newer is recommended for the frontend)
- Git
- On Ubuntu with CUDA, a compatible NVIDIA driver and CUDA runtime

The setup scripts create isolated environments and install the required packages:

| Environment | Contents |
| --- | --- |
| `.runtimes/reinvent` | REINVENT4, PyTorch, RDKit, scikit-learn, LightGBM, XGBoost, CatBoost |
| `webapp/backend/.venv` | FastAPI backend and web application dependencies |
| `.runtimes/chemprop-v2` | Chemprop v2.3.1 |
| `.runtimes/qsartuna` | QSARtuna and its pinned dependencies |
| `webapp/frontend/node_modules` | Frontend dependencies from `package-lock.json` |

### macOS Apple Silicon

```bash
git clone https://github.com/keetane/REINVENT4.git
cd REINVENT4
bash scripts/setup_macos.sh
python3 scripts/check_environment.py
bash scripts/start_webapp.sh
```

### Ubuntu CPU

```bash
git clone https://github.com/keetane/REINVENT4.git
cd REINVENT4
bash scripts/setup_ubuntu.sh
python3 scripts/check_environment.py
bash scripts/start_webapp.sh
```

### Ubuntu NVIDIA CUDA

Set the PyTorch processor type before running the Ubuntu setup. Use the CUDA version appropriate for the installed driver, for example:

```bash
REINVENT_PROCESSOR_TYPE=cu124 bash scripts/setup_ubuntu.sh
```

The setup scripts use `-d none` for the REINVENT installation. OpenEye ROCS and other license-dependent packages are therefore **not installed or used**. Structure and similarity calculations use open-source components such as RDKit.

The scripts do not install Python, Node.js, npm, system compilers, or GPU drivers. Those prerequisites must already be available on the machine.

### Verify the installation

```bash
python3 scripts/check_environment.py
curl http://127.0.0.1:8000/api/health
```

The launcher starts the API, worker, and frontend. To stop them, use the process controls provided by the launcher or stop the PIDs recorded in `webapp/.run/`.

## Model and data workflow

1. Prepare a generation model, TL model, or QSAR model.
2. Select the model and device in the relevant tab.
3. Start the asynchronous job.
4. Inspect progress and logs in **Jobs & History**.
5. Inspect structures, properties, distributions, correlations, and model scores in **Results**.

QSAR files can be uploaded using names such as `yyyymmdd_QSAR_RandomForest_comment.joblib`. Generation model uploads use `yyyymmdd_TL_Mol2MolHigh_comment.model` or `yyyymmdd_RL_Reinvent_comment.chkpt`. The filename supplies the model type, purpose, algorithm, and display comment automatically; duplicate names receive a serial suffix.

Generated CSV downloads use the form `yyyymmdd_MolGen_<ModelType>_<comment>.csv`, for example `20261019_MolGen_Mol2MolHigh_aspirin.csv`.

## Storage locations

By default, artifacts are stored under the repository root:

```text
models/                  TL/RL generation models
models/qsar/             QSAR model files
jobs/<job-id>/input/     uploaded and standardized inputs
jobs/<job-id>/output/    generated CSV/SDF, predictions, evaluations
jobs/<job-id>/logs/      worker and runtime logs
jobs/<job-id>/config/    job configuration snapshots
jobs.db                  persistent job and model registry
```

The locations can be overridden with `REINVENT_WEB_REPOSITORY_ROOT`, `REINVENT_WEB_JOBS_DIR`, `REINVENT_WEB_MODELS_DIR`, and `REINVENT_WEB_DATABASE_PATH`.

## Open-source software

The application uses REINVENT4, RDKit, PyTorch, scikit-learn, Chemprop, QSARtuna, LightGBM, XGBoost, CatBoost, FastAPI, and related open-source packages. See the Tutorial tab for versions, licenses, source links, and documentation links. OpenEye ROCS is not used.

## Upstream REINVENT4 documentation

The sections below preserve the upstream REINVENT4 background, command-line usage, prior-model references, scoring-plugin notes, and development information. They describe the underlying REINVENT4 project; for this repository's web application, use the setup and launch commands above.

REINVENT 4
==========


Description
-----------

REINVENT is a molecular design tool for de novo design, scaffold hopping,
R-group replacement, linker design, molecule optimization, and other small
molecule design tasks.  REINVENT uses a Reinforcement Learning
(RL) algorithm to generate optimized molecules compliant with a user defined
property profile defined as a multi-component score.  Transfer Learning (TL)
can be used to create or pre-train a model that generates molecules closer
to a set of input molecules. 

A paper describing the software has been published as Open Access in the
Journal of Cheminformatics:
[Reinvent 4: Modern AI–driven generative molecule design](https://link.springer.com/article/10.1186/s13321-024-00812-5?utm_source=rct_congratemailt&utm_medium=email&utm_campaign=oa_20240221&utm_content=10.1186/s13321-024-00812-5).
See AUTHORS.md for references to previous papers.


Requirements
------------

REINVENT is being developed on Linux and supports both GPU and CPU.  The Linux
version is fully validated.  REINVENT on Windows and MacOSX supports
both GPU and CPU but is only partially tested on these platforms and therefore
support is limited.

The code is written in Python 3 (>= 3.10).  The list of
dependencies can be found in the repository (see also Installation below).

A GPU is not strictly necessary but strongly recommended for performance
reasons especially for transfer learning and model training.  Reinforcement
learning (RL) requires the computation of scores where most scoring
components run on the CPU.  Thus, a GPU is less important for RL (depending
on how much time is spent on the CPU).

Note that if no GPU is installed in your computer the code will run on the
CPU automatically.  REINVENT [supports](https://pytorch.org/get-started/locally/) NVIDIA GPUs and also some AMD GPUs.
For most design tasks a memory of about 8 GiB for both CPU main memory and
GPU memory is sufficient.


Upstream REINVENT4 CLI installation (optional)
----------------------------------------------

The web application setup above is the supported installation path for this repository. If you also need to run the underlying REINVENT4 command-line tool directly, follow the upstream installation guidance and select the PyTorch processor type appropriate for the machine:

```shell
python install.py --help
python install.py cu124  # or rocm6.2.4, cpu, mac, etc.
reinvent --help
```

The web application itself uses isolated runtimes created by `scripts/setup_macos.sh` or `scripts/setup_ubuntu.sh`. Do not install OpenEye ROCS for the web application; it is not used by this project.

Basic Usage
-----------

REINVENT is a command line tool and works principally as follows
```shell
reinvent -l sampling.log sampling.toml
```

This writes logging information to the file `sampling.log`.  If you wish to write
this to the screen, leave out the `-l sampling.log` part. `sampling.toml` is the
configuration file.  The main user format is [TOML](https://toml.io/en/) as it tends to be more
use friendly.  JSON can be used too, add `-f json`, but a specialised editor is
recommended as the format is very sensitive to minor changes.

Sample configuration files for all run modes are
located in `configs/toml` in the repository and file paths in these files would need to be
adjusted to your local installation.  In particular, ready made prior models are
located in `priors` and you would choose a model and the
appropriate run mode depending on the research problem you are trying to address.
There is additional information in `configs/toml` in several `*.md` files with
instructions on how to configure the TOML file.  Internal priors can be referenced with a
dot notation (see `reinvent/prior_registry.py`).


Prior models
------------

All public prior models can be found on [Zenodo](https://doi.org/10.5281/zenodo.15641296).


Tutorials / `Jupyter` notebooks
-------------------------------

Basic instructions can be found in the comments in the config examples in `configs/toml`.

Notebooks are provided in the `notebooks/` directory.  Please note that we
provide the notebooks in jupytext "light script" format.  To work with the light
scripts you will need to install jupytext.  A few other packages will come in handy too.

```shell
pip install jupytext mols2grid seaborn
```

The Python files in `notebooks/` can then be converted to a notebook e.g.

```shell
jupytext -o Reinvent_demo.ipynb Reinvent_demo.py
```


Updating dependencies
---------------------

Update the lock files with [pip-tools](https://pypi.org/project/pip-tools/) (please, do not edit the files manually):
```shell
pip-compile --extra-index-url=https://download.pytorch.org/whl/cu121 --extra-index-url=https://pypi.anaconda.org/OpenEye/simple --resolver=backtracking pyproject.toml
```
To update a single package, use `pip-compile --upgrade-package somepackage`
(see the documentation for pip-tools).


Scoring Plugins
---------------

The scoring subsystem uses a simple plugin mechanism (Python
[native namespace packages](https://packaging.python.org/en/latest/guides/packaging-namespace-packages/#native-namespace-packages)).  If you
wish to write your own plugin, follow the instructions below.
There is no need to touch any of the REINVENT code. The public
repository contains a [contrib](https://github.com/MolecularAI/REINVENT4/tree/main/contrib/reinvent_plugins/components) directory with some useful examples.

1. Create `/top/dir/somewhere/reinvent\_plugins/components` where `/top/dir/somewhere` is a convenient location for you.
1. Do **not** place a `__init__.py` in either `reinvent_plugins` or `components` as this would break the mechanism.  It is fine to create normal packages within `components` as long as you import those correctly.
1. Place a file whose name starts with `comp_*` into `reinvent_plugins/components`.   Files with different names will be ignored i.e. not imported. The directory will be searched recursively so structure your code as needed but directory/package names must be unique.
1. Tag the scoring component class(es) in that file with the @add\_tag decorator.  More than one component class can be added to the same *comp\_* file. See existing code.
1. Tag at most one dataclass for parameters in the same file, see existing code.  This is optional.
1. Set or add `/top/dir/somewhere` to the `PYTHONPATH` environment variable or use any other mechanism to extend `sys.path`.
1. The scoring component should now automatically be picked up by REINVENT.


Unit and Integration Tests 
--------------------------

This is primarily for developers and admins/users who wish to ensure that the
installation works.  The information here is not relevant to the practical use
of REINVENT.  Please refer to _Basic Usage_ for instructions on how to use the 
`reinvent` command.

The REINVENT project uses the `pytest` framework for its tests.  Before you run
them you first have to create a configuration file for the tests.

In the project directory, create a `config.json` file in the `configs/` directory.
You can use the example config `example.config.json` as a base.  Make sure that
you set `MAIN_TEST_PATH` to a non-existent directory.  That is where temporary
files will be written during the tests.  If it is set to an existing directory,
that directory will be removed once the tests have finished.

Some tests require a proprietary OpenEye license.  You have to set up a few
things to make the tests read your license.  The simple way is to just set the
`OE_LICENSE` environment variable to the path of the file containing the
license.  

Once you have a configuration and your license can be read, you can run the tests.

```
$ pytest tests --json /path/to/config.json --device cuda
```
