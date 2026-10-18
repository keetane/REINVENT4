# Direct third-party dependencies

This inventory records direct application and test dependencies pinned by this project. Transitive
packages and their license files are installed through the lock files and should be regenerated and
reviewed before redistribution.

| Dependency | Pinned version | License |
| --- | ---: | --- |
| FastAPI | 0.143.2 | MIT |
| Pydantic | 2.11.3 | MIT |
| Uvicorn | 0.38.0 | BSD-3-Clause |
| RDKit | 2024.9.6 | BSD-3-Clause |
| HTTPX | 0.28.1 | BSD-3-Clause |
| pytest | 8.4.1 | MIT |
| React | 18.2.0 | MIT |
| React DOM | 18.2.0 | MIT |
| i18next | 26.4.7 | MIT |
| react-i18next | 17.0.21 | MIT |
| Ketcher Core | 3.18.0 | Apache-2.0 |
| Ketcher React | 3.18.0 | Apache-2.0 |
| Ketcher Standalone | 3.18.0 | Apache-2.0 |
| events | 3.3.0 | MIT |
| process | 0.11.10 | MIT |
| Vite | 7.3.7 | MIT |
| Vitest | 4.1.11 | MIT |
| Playwright test | 1.55.1 | Apache-2.0 |
| TypeScript | 5.9.3 | Apache-2.0 |
| REINVENT4 | 4.5.11 | Apache-2.0 |
| scikit-learn | 1.6.1 | BSD-3-Clause |
| joblib | 1.4.2 | BSD-3-Clause |
| NumPy (REINVENT runtime) | 1.x installed runtime | BSD-3-Clause |
| PyTorch | 2.6.0 | BSD-3-Clause |
| Chemprop (isolated runtime) | 2.3.1 | MIT |
| QSARtuna (isolated runtime) | 4.0.1, commit `664f73e3a9f68aaf788ee6a92a4889778ddeaf84` | Apache-2.0 declared in `pyproject.toml`; see note below |
| LightGBM | 4.7.0 | MIT |
| XGBoost | 3.2.0 | Apache-2.0 |
| CatBoost | 1.2.10 | Apache-2.0 |

## Commercial use and redistribution notes

- Apache-2.0 projects permit commercial use and hosted services. Redistribution requires a copy of
  the license, preservation of applicable copyright/attribution notices, prominent notices on
  modified files, and preservation of an upstream NOTICE file when one exists. Patent and trademark
  clauses remain applicable.
- BSD-3-Clause projects permit commercial use and hosted services. Source and binary redistribution
  must preserve the copyright notice, conditions, and disclaimer. Project or contributor names may
  not be used to endorse a product without permission.
- MIT projects permit commercial use and hosted services. Redistribution must preserve the
  copyright and permission notice in copies or substantial portions.
- SaaS use without distributing the libraries generally does not itself trigger source-distribution
  obligations under these permissive licenses. Container, desktop bundle, wheel, source, or compiled
  frontend distribution does. This inventory is operational guidance, not legal advice.

REINVENT4's repository license is Apache-2.0; RDKit and scikit-learn declare BSD-3-Clause; and
Chemprop declares MIT. The exact license texts from the versions actually shipped must be included
in a distributable `licenses/` bundle rather than relying only on this summary.

QSARtuna 4.0.1 is installed from the official repository at the pinned commit shown above and is
used for AutoML training and inference. Its `pyproject.toml` declares Apache-2.0, but the reviewed
commit does not contain a standalone LICENSE file. Before redistribution, obtain and preserve the
complete license materials for the pinned source and audit its transitive dependencies. This caveat
does not change the application's runtime capability; it is a packaging/legal review item.

OpenEye toolkits and ROCS are optional proprietary dependencies upstream of REINVENT4. They are not
installed, invoked, or required by this Web app. Do not select the upstream `openeye` optional extra
or ship OpenEye model/query assets without a separate commercial-license review.

Chemprop v2 and QSARtuna are installed into isolated runtimes. This isolation prevents their newer
PyTorch/RDKit dependency set from changing the pinned REINVENT4 runtime, but it does not remove any
license or attribution duties when the environments are redistributed.

Ketcher is used locally and is not connected to an external structure service. Its Apache-2.0
license is declared by each installed Ketcher package. Preserve required license and notice text
when distributing a compiled frontend or dependency bundle.
