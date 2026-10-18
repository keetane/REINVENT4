# QSAR / Reinforcement Learning change history

## 2026-10-18

- Added LightGBM, XGBoost, CatBoost, Kernel Ridge, Gaussian Process, Ridge, Elastic Net,
  Logistic Regression, and MLP with regression/classification capability filtering.
- Added isolated Chemprop v2 and QSARtuna runtimes, real training/prediction workers, saved-model
  reload, and REINVENT ExternalProcess scoring integration.
- Added Model Comparison with one shared split, training-only cross-validation and grid search,
  partial-failure handling, recommendation, comparison plots, and selected-model registration.
- Added multiple weighted QSAR scoring targets, maximize/minimize/range/target transforms,
  applicability-domain behavior, and a server-side scoring preview that includes RDKit properties
  in the displayed Total Score.
- Added Japanese/English switching with persisted language, translated forms, tables, charts,
  model descriptions, Tutorial, setup guide, job status, QSAR, and RL screens.
- Added macOS/Ubuntu setup scripts, isolated dependency checks, OSS source/license links, and an
  explicit note that OpenEye ROCS is not used.
- Fixed runtime discovery so virtual-environment Python symlinks keep their own `sys.prefix`.

## 2026-10-17

- Added persistent QSAR model registry and job types for QSAR training and Reinforcement Learning.
- Added CSV inspection, RDKit validation/canonicalization, concentration-to-pActivity conversion,
  duplicate conflict policies, rejected-row export, and cleaned dataset snapshots.
- Added Random Forest, SVM, Gradient Boosting, and Extra Trees for regression/classification with
  Pipeline-based preprocessing, random/scaffold splits, cross-validation, Test metrics, predictions,
  and Morgan-Tanimoto applicability-domain output.
- Added capability reporting that disables an optional runtime with a concrete environment reason.
- Added trusted batched QSAR prediction and the REINVENT4 `QsarJoblib` scoring plugin.
- Added Single-stage DAP RL, official transforms/aggregation, diversity filter, responder progress,
  checkpoint/summary/config/log persistence, RL model registration, and generation reuse.
- Added tabs `03 QSAR`, `04 強化学習`, renumbered job/results tabs, and extended common job status.
- Added automated backend/frontend coverage plus real scikit-learn and one-step REINVENT4 RL smoke
  verification.
