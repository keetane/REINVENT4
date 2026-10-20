interface TutorialItem {
  title: string;
  body: string;
}

interface ModelChoice extends TutorialItem {
  key: "standard" | "transfer_learning" | "reinforcement_learning";
}

export interface TutorialContent {
  title: string;
  intro: string;
  environmentTitle: string;
  environmentIntro: string;
  generationTitle: string;
  generationIntro: string;
  modelChoiceTitle: string;
  modelChoices: ModelChoice[];
  generationSteps: TutorialItem[];
  sampleTitle: string;
  noParent: string;
  transferTitle: string;
  transferIntro: string;
  transferSteps: TutorialItem[];
  transferResultsTitle: string;
  transferResults: TutorialItem[];
  transferParametersTitle: string;
  transferParameters: TutorialItem[];
  openTraining: string;
  qsarTitle: string;
  qsarIntro: string;
  qsarWorkflowTitle: string;
  qsarWorkflowSteps: string[];
  qsarUses: TutorialItem[];
  qsarScoresTitle: string;
  qsarScores: TutorialItem[];
  rlTitle: string;
  rlIntro: string;
  rlSteps: TutorialItem[];
  rlWorkflowTitle: string;
  rlWorkflowSteps: string[];
  labels: {
    modelSource: string;
    standardSample: string;
    use: string;
    copy: string;
    copied: string;
    mac: string;
    ubuntu: string;
    verify: string;
    source: string;
    docs: string;
  };
  generators: { key: string; name: string; example: string }[];
  installation: { title: string; mac: string[]; ubuntu: string[]; verify: string[] };
  ossTitle: string;
  noRocs: string;
}

export const tutorialContent: Record<"ja" | "en", TutorialContent> = {
  ja: {
    title: "Molecular Generatorの使い方",
    intro: "環境構築から分子生成、Transfer Learning、QSAR、強化学習まで、各タブの役割と結果の見方を順に説明します。",
    environmentTitle: "環境構築",
    environmentIntro: "アプリ本体と計算環境を準備します。Chemprop v2とQSARtunaは依存関係の競合を避けるため、独立した実行環境にインストールされます。",
    generationTitle: "分子生成",
    generationIntro: "モデルが学習している化学空間から分子をサンプリングします。使用するモデルの由来によって、生成される分子の傾向が変わります。",
    modelChoiceTitle: "モデルの種類",
    modelChoices: [
      { key: "standard", title: "Standardモデル", body: "すぐに分子生成へ使用できる標準Priorです。Reinvent、LibInvent、LinkInvent、Mol2Molなど各モデル固有の説明は、分子生成タブのModel選択欄に表示されます。" },
      { key: "transfer_learning", title: "Transfer Learningモデル", body: "SMILESデータから構造パターンを追加学習した生成モデルです。学習データに多い骨格・部分構造・化学空間を反映した分子が生成されます。" },
      { key: "reinforcement_learning", title: "Reinforcement Learningモデル", body: "Scoring Functionが高くなるよう最適化した生成モデルです。QSAR活性や物性など、学習時に設定した評価基準を反映した分子が生成されます。" }
    ],
    generationSteps: [
      { title: "Model", body: "Standard、TL、RLから目的に合うモデルを選びます。選択すると、Parentの要否とモデルの説明が同じ欄に表示されます。" },
      { title: "Parent", body: "ReinventはParent不要です。Mol2Mol、LibInvent、LinkInventなどは、モデルに合ったParent SMILESまたは接続点を入力します。" },
      { title: "Sampling", body: "生成数、コメント、デバイスを設定してジョブを開始します。完了後は結果解析タブで構造と物性を確認できます。" }
    ],
    sampleTitle: "Standardモデルの入力サンプル",
    noParent: "ReinventはParent SMILESを入力せずに使用します。",
    transferTitle: "Transfer Learning",
    transferIntro: "既知化合物のSMILESから構造情報を学習し、目的の化学空間に特化した生成モデルを作成します。活性値は必須ではなく、SMILESの構造パターンそのものを学びます。",
    transferSteps: [
      { title: "学習データ", body: "SMILES/SMI/TXTまたはCSVを読み込みます。CSVではSMILES列を選択でき、不正構造と重複は標準化時に除外されます。" },
      { title: "Base Model", body: "Reinventは分子全体の生成分布を学びます。Mol2Molは入力分子から別の分子へ変換するパターンを学びます。" },
      { title: "Training / Validation", body: "データを学習用と検証用に分けます。Validationは未学習データでの再現性を確認し、過学習を見つけるために使います。" }
    ],
    transferResultsTitle: "学習結果の見方",
    transferResults: [
      { title: "Lossの収束", body: "Training LossとValidation Lossがともに低下し、変化が小さくなれば学習は収束傾向です。絶対値だけでなく推移を確認します。" },
      { title: "過学習", body: "Training Lossだけが下がり続け、Validation Lossが上昇し始めた場合は過学習の可能性があります。上昇前後のEpochを目安にします。" },
      { title: "生成確認", body: "完成したモデルで分子を生成し、妥当な構造の割合、重複、目的の骨格や部分構造が反映されているかを確認します。" }
    ],
    transferParametersTitle: "パラメータ調節",
    transferParameters: [
      { title: "Epoch", body: "不足すると学習が浅く、多すぎると過学習しやすくなります。Validation Lossが悪化する前の範囲へ調整します。" },
      { title: "Learning rate", body: "Lossが不安定なら下げ、ほとんど改善しない場合は少し上げます。一度に大きく変更せず段階的に調節します。" },
      { title: "Batch size", body: "大きくすると計算は安定しやすい一方、メモリ消費が増えます。メモリ不足では値を下げます。" },
      { title: "Tanimoto範囲", body: "Mol2Molで学習する構造変換の近さを設定します。下限を下げるほど大きな変換を含み、上限を上げるほど近い変換を含みます。" }
    ],
    openTraining: "Transfer Learningを開く",
    qsarTitle: "QSAR",
    qsarIntro: "SMILESと活性値・物性値の対応から、未知構造の値を予測するモデルを作成します。作成したモデルは単一分子の予測と、強化学習のScoring Functionの両方に使用できます。",
    qsarWorkflowTitle: "QSARモデル作成の流れ",
    qsarWorkflowSteps: ["SMILES＋目的変数", "標準化・分割", "学習・比較", "Test評価", "モデル登録"],
    qsarUses: [
      { title: "Single Model", body: "アルゴリズム、特徴量、分割方法、交差検証条件を指定して1つのモデルを学習します。" },
      { title: "Model Comparison", body: "同じデータ分割で複数アルゴリズムを比較します。推奨モデルは選択した主要指標を基準に決まります。" },
      { title: "活性予測", body: "登録済みモデルの「新規構造を予測」から、描画または入力したSMILESの予測値と適用領域を確認できます。" },
      { title: "RL Scoring", body: "推論可能なQSARモデルをRLのScoring Functionに追加し、予測活性を考慮した生成モデルを学習できます。" }
    ],
    qsarScoresTitle: "評価スコアの見方",
    qsarScores: [
      { title: "R²", body: "回帰モデルの説明力です。1に近いほど良く、0以下は単純な平均予測より悪い可能性を示します。" },
      { title: "RMSE / MAE", body: "予測誤差の大きさです。小さいほど良好です。RMSEは大きな誤差をより強く評価し、MAEは平均的なずれを表します。" },
      { title: "ROC-AUC / F1", body: "分類モデルの指標です。ROC-AUCは順位付け性能、F1は適合率と再現率のバランスを示し、どちらも1に近いほど良好です。" },
      { title: "Train / Validation / Test", body: "Trainだけが良い場合は過学習を疑います。モデル比較はValidationで選び、最後にTestで未知データへの性能を確認します。" },
      { title: "Applicability Domain", body: "学習化合物との類似度から予測の信頼範囲を示します。Outの場合は数値を参考値として慎重に扱います。" }
    ],
    rlTitle: "強化学習",
    rlIntro: "生成分子のScoring Functionを改善する方向へAgentを更新する学習です。設定した評価項目の総合スコアが高い構造を生成しやすいモデルを作成します。",
    rlSteps: [
      { title: "初期Agent", body: "StandardモデルまたはTransfer Learningモデルを出発点にします。TLモデルを使うと、対象化学空間を保ちながら最適化を始められます。" },
      { title: "Scoring Function", body: "QSAR予測、MW、TPSA、QEDなどを組み合わせ、重みと変換方法を設定します。目的と制約を一つの総合スコアにまとめます。" },
      { title: "QSARモデル", body: "QSARタブで作成またはアップロードして登録したモデルを選択します。活性値などの予測結果を報酬へ組み込みます。" },
      { title: "結果確認", body: "Scoreの改善だけでなく、妥当な構造の割合、重複、化学的多様性、各Scoring項目の推移も確認します。" }
    ],
    rlWorkflowTitle: "QSARを利用した強化学習",
    rlWorkflowSteps: ["生成モデル", "分子生成", "QSAR・物性評価", "総合Score", "Agent更新", "最適化モデル"],
    labels: { modelSource: "モデルの由来", standardSample: "Standardモデル", use: "分子生成で試す", copy: "コピー", copied: "コピー済み", mac: "macOS", ubuntu: "Ubuntu", verify: "動作確認", source: "Source", docs: "Docs" },
    generators: [
      { key: "reinvent", name: "Reinvent", example: "" },
      { key: "libinvent", name: "LibInvent", example: "N#CC[C@H](C1CCCC1)[*:1]" },
      { key: "linkinvent", name: "LinkInvent", example: "Oc1cncc(*)c1|*c1ccoc1" },
      { key: "mol2mol_high", name: "Mol2Mol", example: "CC(=O)OC1=CC=CC=C1C(=O)O" }
    ],
    installation: {
      title: "対応: macOS Apple Silicon / Ubuntu CPU・CUDA（利用可能な場合）",
      mac: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_macos.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      ubuntu: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_ubuntu.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      verify: ["python3 --version", "python3 scripts/check_environment.py", "curl http://127.0.0.1:8000/api/health", "ブラウザ: http://127.0.0.1:5173/"]
    },
    ossTitle: "使用している主なOSS",
    noRocs: "OpenEye ROCSは使用していません。構造・類似度計算はRDKit等のOSSで実装しています。"
  },
  en: {
    title: "Using Molecular Generator",
    intro: "Follow the workflow from environment setup through molecule generation, Transfer Learning, QSAR, and Reinforcement Learning, with guidance for interpreting results.",
    environmentTitle: "Environment setup",
    environmentIntro: "Prepare the application and its compute runtimes. Chemprop v2 and QSARtuna are installed in isolated environments to avoid dependency conflicts.",
    generationTitle: "Molecule Generation",
    generationIntro: "Sample molecules from the chemical space learned by a model. The origin of the selected model determines the tendencies of the generated structures.",
    modelChoiceTitle: "Model source",
    modelChoices: [
      { key: "standard", title: "Standard model", body: "A ready-to-use Prior. Model-specific guidance for Reinvent, LibInvent, LinkInvent, and Mol2Mol is shown directly in the Model selector on the Molecule Generation tab." },
      { key: "transfer_learning", title: "Transfer Learning model", body: "A generator fine-tuned on SMILES data. Generated molecules reflect the scaffolds, substructures, and chemical space represented in the training set." },
      { key: "reinforcement_learning", title: "Reinforcement Learning model", body: "A generator optimized for a Scoring Function. Generated molecules reflect the QSAR activity, properties, and other objectives used during training." }
    ],
    generationSteps: [
      { title: "Model", body: "Choose a Standard, TL, or RL model. Its description and Parent requirement are displayed in the same selector area." },
      { title: "Parent", body: "Reinvent needs no Parent. Mol2Mol, LibInvent, and LinkInvent require the appropriate Parent SMILES or attachment points." },
      { title: "Sampling", body: "Set the number of molecules, comment, and device, then start the job. Inspect structures and properties in Results when it completes." }
    ],
    sampleTitle: "Standard model input samples",
    noParent: "Reinvent is used without a Parent SMILES.",
    transferTitle: "Transfer Learning",
    transferIntro: "Learn structural information from known SMILES and create a generator specialized for a target chemical space. Activity values are not required; the model learns patterns in the structures themselves.",
    transferSteps: [
      { title: "Training data", body: "Load SMILES/SMI/TXT or CSV. Choose the SMILES column for CSV data; invalid structures and duplicates are removed during standardization." },
      { title: "Base Model", body: "Reinvent learns a distribution over complete molecules. Mol2Mol learns transformations from input molecules to related outputs." },
      { title: "Training / Validation", body: "Split the data into training and validation sets. Validation measures generalization to unseen data and helps detect overfitting." }
    ],
    transferResultsTitle: "Reading training results",
    transferResults: [
      { title: "Loss convergence", body: "Training and Validation Loss should both decrease and gradually level off. Interpret the trend rather than the absolute value alone." },
      { title: "Overfitting", body: "If Training Loss continues to fall while Validation Loss starts rising, the model may be overfitting. Use the turn in Validation Loss as a guide." },
      { title: "Generation check", body: "Generate molecules and inspect validity, duplicates, and whether the intended scaffolds or substructures are represented." }
    ],
    transferParametersTitle: "Tuning parameters",
    transferParameters: [
      { title: "Epoch", body: "Too few epochs underfit; too many increase overfitting risk. Prefer the range before Validation Loss deteriorates." },
      { title: "Learning rate", body: "Lower it when Loss is unstable; raise it slightly when learning barely progresses. Adjust gradually." },
      { title: "Batch size", body: "Larger batches can be more stable but use more memory. Reduce the value when memory is insufficient." },
      { title: "Tanimoto range", body: "For Mol2Mol, this controls transformation similarity. A lower minimum includes larger changes; a higher maximum includes closer pairs." }
    ],
    openTraining: "Open Transfer Learning",
    qsarTitle: "QSAR",
    qsarIntro: "Train a model that predicts an activity or property from a structure. A registered QSAR model can predict newly drawn molecules and serve as a Reinforcement Learning Scoring Function.",
    qsarWorkflowTitle: "QSAR modeling workflow",
    qsarWorkflowSteps: ["SMILES + target", "Standardize and split", "Train or compare", "Test evaluation", "Register model"],
    qsarUses: [
      { title: "Single Model", body: "Train one model with a selected algorithm, features, split strategy, and cross-validation settings." },
      { title: "Model Comparison", body: "Compare algorithms on the same split. The recommended model is selected by the chosen primary metric." },
      { title: "Activity prediction", body: "Use Predict new structure on a registered model to inspect the prediction and applicability domain for a drawn or entered SMILES." },
      { title: "RL Scoring", body: "Add a runnable QSAR model to an RL Scoring Function to train a generator that considers predicted activity." }
    ],
    qsarScoresTitle: "Interpreting evaluation scores",
    qsarScores: [
      { title: "R²", body: "Regression explanatory power. Values closer to 1 are better; a value at or below 0 can be worse than predicting the mean." },
      { title: "RMSE / MAE", body: "Prediction error magnitudes; lower is better. RMSE penalizes large errors more strongly, while MAE represents the average absolute deviation." },
      { title: "ROC-AUC / F1", body: "Classification metrics. ROC-AUC measures ranking performance and F1 balances precision and recall; both improve toward 1." },
      { title: "Train / Validation / Test", body: "Strong Train results alone can indicate overfitting. Select models on Validation, then use Test once for final unseen-data performance." },
      { title: "Applicability Domain", body: "Estimates whether a prediction is supported by similar training compounds. Treat Out-of-domain predictions cautiously." }
    ],
    rlTitle: "Reinforcement Learning",
    rlIntro: "Update an Agent to improve its Scoring Function. The resulting model becomes more likely to generate structures with high combined scores for the selected objectives.",
    rlSteps: [
      { title: "Initial Agent", body: "Start from a Standard or Transfer Learning model. A TL model helps optimization remain near the target chemical space." },
      { title: "Scoring Function", body: "Combine QSAR predictions with MW, TPSA, QED, and other properties, then configure weights and transformations." },
      { title: "QSAR model", body: "Choose a model trained or uploaded and registered on the QSAR tab. Its activity or property prediction becomes part of the reward." },
      { title: "Review results", body: "Check Score improvement together with validity, duplicates, chemical diversity, and the trend for each scoring component." }
    ],
    rlWorkflowTitle: "Reinforcement Learning with QSAR",
    rlWorkflowSteps: ["Generator", "Generate molecules", "QSAR and properties", "Combined Score", "Update Agent", "Optimized model"],
    labels: { modelSource: "Model source", standardSample: "Standard model", use: "Try in generation", copy: "Copy", copied: "Copied", mac: "macOS", ubuntu: "Ubuntu", verify: "Verification", source: "Source", docs: "Docs" },
    generators: [
      { key: "reinvent", name: "Reinvent", example: "" },
      { key: "libinvent", name: "LibInvent", example: "N#CC[C@H](C1CCCC1)[*:1]" },
      { key: "linkinvent", name: "LinkInvent", example: "Oc1cncc(*)c1|*c1ccoc1" },
      { key: "mol2mol_high", name: "Mol2Mol", example: "CC(=O)OC1=CC=CC=C1C(=O)O" }
    ],
    installation: {
      title: "Supported: Apple Silicon macOS / Ubuntu CPU and CUDA where available",
      mac: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_macos.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      ubuntu: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_ubuntu.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      verify: ["python3 --version", "python3 scripts/check_environment.py", "curl http://127.0.0.1:8000/api/health", "Browser: http://127.0.0.1:5173/"]
    },
    ossTitle: "Primary open-source software",
    noRocs: "OpenEye ROCS is not used. Structure and similarity calculations use open-source software such as RDKit."
  }
};
