# QSAR and RL operating guide

## Supported runtime

The verified local runtime is REINVENT4 4.5.11, RDKit 2024.09.6, scikit-learn 1.6.1,
LightGBM 4.7.0, XGBoost 3.2.0, CatBoost 1.2.10, joblib 1.4.2, and PyTorch 2.6.0.
Chemprop 2.3.1 and QSARtuna 4.0.1 run in isolated Python 3.12 environments to avoid dependency
conflicts. The capability API reports actual installed versions and disables a missing runtime.

## QSAR workflow

1. Upload a CSV and select the SMILES and target columns.
2. Choose regression/classification, unit conversion, and duplicate-conflict handling.
3. Review valid/excluded counts, preview rows, and download the exclusion report when needed.
4. Select Single Model or Model Comparison, then choose algorithms, features, split,
   cross-validation folds, and AD threshold. Chemprop uses its molecular graph and QSARtuna runs
   AutoML, so neither requires a fingerprint choice.
5. Start the asynchronous job and monitor it under `05 進捗と履歴`.
6. Inspect Test metrics and download the model, cleaned dataset, evaluation JSON, and prediction CSV.

Only the Training partition fits the feature transformer, imputer, scaler, and estimator. Test data
is not used for model selection. Scaffold split keeps a Bemis–Murcko scaffold in only one partition;
datasets that cannot form two scaffold groups must use Random Split. Classification random split is
stratified where scikit-learn can satisfy the requested partition sizes.

The initial AD method is nearest-neighbour Morgan Tanimoto (radius 2, 2048 bits). It is a domain
heuristic, not a prediction-confidence guarantee.

## RL workflow

1. Select a standard Prior, TL model, or prior RL checkpoint as the initial Agent.
2. Select one or more registered QSAR models and optional RDKit property components.
3. Test selected QSAR models and RDKit properties against example SMILES before starting. The
   preview reports raw predictions, AD, transformed score, weight, and the same weighted Total
   Score aggregation used by the job configuration.
4. Configure DAP parameters and the diversity filter under Advanced Settings.
5. Start the server-side RL job. The worker writes an immutable TOML snapshot, log, summary CSV,
   TensorBoard data, and checkpoint.
6. On successful completion the checkpoint appears in the generation model selector. Its summary
   CSV can be opened by the existing structure-analysis view.

The `QsarJoblib` REINVENT scoring plugin accepts only model paths resolved by the server from the
QSAR registry. Chemprop v2 and QSARtuna use a fixed server-owned `ExternalProcess` predictor and
runtime. The API never accepts an executable command or arbitrary model path. joblib files must be
treated like executable content; downloaded files must not be uploaded and deserialized.

## Limits and unverified environments

- Chemprop currently trains one selected target per model. Chemprop multitask output heads and
  predictive uncertainty are not exposed yet.
- A complete staged-learning run using a Chemprop or QSARtuna scorer was not executed in the final
  local verification; both saved-model predictors and their ExternalProcess adapters were exercised.
- CPU and Apple Silicon MPS are exposed for REINVENT/TL/RL; scikit-learn QSAR runs on CPU.
- The implementation was exercised locally on Apple Silicon. The provided clean setup scripts were
  not run on a fresh Mac, Ubuntu CPU, or Ubuntu CUDA host in this workspace.
- Checkpoint resume is not advertised. REINVENT writes a checkpoint on completion/interruption, but
  this UI does not claim resume support.
- Resource serialization uses the existing single worker queue. Run only one worker process when a
  single GPU/MPS device is shared.
