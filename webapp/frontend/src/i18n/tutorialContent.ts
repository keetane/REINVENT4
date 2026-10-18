export interface TutorialContent {
  title: string; intro: string; chooseTitle: string; step2Title: string; qsarTitle: string; workflowTitle: string;
  sampling: string; tl: string; rl: string; tlrl: string; qsar: string; installationTitle: string; installationIntro: string;
  decisionTree: string[]; noParent: string; openTraining: string; workflowSteps: string[];
  modelingSteps: { title: string; body: string }[];
  labels: { purpose: string; research: string; input: string; parent: string; output: string; use: string; mac: string; ubuntu: string; verify: string; source: string; docs: string };
  generators: { key: string; name: string; purpose: string; research: string; input: string; parent: string; output: string; example: string }[];
  installation: { title: string; mac: string[]; ubuntu: string[]; verify: string[] };
  ossTitle: string; noRocs: string;
}

export const tutorialContent: Record<"ja" | "en", TutorialContent> = {
  ja: {
    title: "分子生成をはじめる", intro: "目的からGeneratorと学習・最適化方法を選び、QSARをRLの予測スコアとして組み合わせます。",
    chooseTitle: "STEP 1：どのような分子を生成したいか？", step2Title: "STEP 2：モデルを学習・最適化する必要があるか？",
    qsarTitle: "QSARの位置付け", workflowTitle: "実測データから分子生成まで",
    decisionTree: ["Parentなし → Reinvent", "固定Scaffold + R-group → LibInvent", "2 Fragments + Linker → LinkInvent", "Parent類似体 → Mol2Mol"],
    labels: { purpose: "機能", research: "適した研究", input: "入力", parent: "Parent", output: "出力", use: "このモデルを使う", mac: "macOS", ubuntu: "Ubuntu", verify: "動作確認", source: "Source", docs: "Docs" },
    noParent: "Parent SMILESを入力せずに利用できます。", openTraining: "Transfer Learningを開く",
    workflowSteps: ["実測データ", "QSAR学習", "予測モデル", "RL Scoring Function", "最適化モデル", "分子生成"],
    modelingSteps: [
      { title: "Single Model", body: "アルゴリズム、特徴量、Scaffold/Random分割、交差検証条件を指定して1モデルを学習します。" },
      { title: "Model Comparison", body: "共通分割とTraining内部CVで比較し、選択指標で推奨した後、Testで1回だけ評価します。" },
      { title: "Chemprop v2 / QSARtuna", body: "独立したPython環境で実行します。Chempropは分子グラフを使い、QSARtunaはAutoMLを行います。" },
      { title: "RL scoring", body: "登録済みで実行可能なQSAR目的変数だけを表示し、MW、TPSA、QED等の物性スコアと組み合わせます。" }
    ],
    sampling: "Sampling：既存モデルを変更せず、そのまま分子を生成します。",
    tl: "Transfer Learning：既知化合物群のSMILESから構造的特徴を学び、目的の化学空間に特化した生成モデルを作ります。",
    rl: "Reinforcement Learning：生成分子をScoring Functionで評価し、高スコア分子の生成確率が高くなるようAgentを更新します。",
    tlrl: "TL → RL：TLで目的の化学空間を学習した後、QSAR活性や物性で最適化します。TLとQSARは独立して作成でき、RLで組み合わせられます。",
    qsar: "QSARは分子生成モデルではなく、構造から活性・物性を予測するモデルです。登録済みで推論可能なQSARモデルだけがRLのScoring Function候補になります。",
    generators: [
      { key: "reinvent", name: "Reinvent", purpose: "Parentを指定せず新規分子を生成", research: "広い化学空間の探索、de novo design", input: "生成件数・モデル・デバイス", parent: "不要", output: "Priorが学習した分布に基づく多様な分子", example: "REINVENT.priorで100分子をSampling" },
      { key: "libinvent", name: "LibInvent", purpose: "固定Scaffoldの接続点からR-groupを生成", research: "母核を保持した側鎖最適化、SAR展開", input: "接続点 * を持つScaffold SMILES", parent: "必要（Scaffold）", output: "Scaffoldを維持し側鎖が異なる分子", example: "N#CC[C@H](C1CCCC1)[*:1]" },
      { key: "linkinvent", name: "LinkInvent", purpose: "2つの固定Fragment間にLinkerを生成", research: "Fragment linking、PROTAC linker探索", input: "各Fragmentに * を1個置き | で区切る", parent: "必要（2 Fragment）", output: "両Fragmentを維持しLinkerが異なる分子", example: "Oc1cncc(*)c1|*c1ccoc1" },
      { key: "mol2mol_high", name: "Mol2Mol", purpose: "Parentに類似した分子を生成", research: "類縁体探索、Matched Molecular Pair、Scaffold hopping", input: "単一分子のParent SMILES", parent: "必須", output: "指定Similarity範囲の構造変換分子", example: "CC(=O)OC1=CC=CC=C1C(=O)O" }
    ],
    installationTitle: "環境構築 / Installation & Setup", installationIntro: "クリーン環境では付属スクリプトを実行します。Chemprop v2とQSARtunaは依存衝突を避ける独立環境に入ります。",
    installation: {
      title: "対応: macOS Apple Silicon / Ubuntu CPU・CUDA（利用可能な場合）",
      mac: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_macos.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      ubuntu: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_ubuntu.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      verify: ["python3 --version", "python3 scripts/check_environment.py", "curl http://127.0.0.1:8000/api/health", "ブラウザ: http://127.0.0.1:5173/"]
    },
    ossTitle: "使用している主なOSS", noRocs: "OpenEye ROCSは使用していません。構造・類似度計算はRDKit等のOSSで実装しています。"
  },
  en: {
    title: "Choosing generation, training and optimization workflows", intro: "Choose a generator and training method from your objective, then combine QSAR predictions with RL scoring.",
    chooseTitle: "STEP 1: What molecules do you want to generate?", step2Title: "STEP 2: Do you need to train or optimize a model?",
    qsarTitle: "Where QSAR fits", workflowTitle: "From measured data to generated molecules",
    decisionTree: ["No parent → Reinvent", "Fixed scaffold + R-group → LibInvent", "2 fragments + linker → LinkInvent", "Parent analogues → Mol2Mol"],
    labels: { purpose: "Function", research: "Best for", input: "Input", parent: "Parent", output: "Output", use: "Use this model", mac: "macOS", ubuntu: "Ubuntu", verify: "Verification", source: "Source", docs: "Docs" },
    noParent: "No Parent SMILES is required.", openTraining: "Open Transfer Learning",
    workflowSteps: ["Measured data", "QSAR training", "Predictor", "RL Scoring Function", "Optimized Agent", "Molecule generation"],
    modelingSteps: [
      { title: "Single Model", body: "Choose one algorithm, features, Scaffold/Random split and cross-validation settings." },
      { title: "Model Comparison", body: "Apply one shared split, compare training-only CV, recommend by the selected metric, then evaluate once on Test." },
      { title: "Chemprop v2 / QSARtuna", body: "Run in isolated Python environments; Chemprop uses molecular graphs and QSARtuna performs AutoML." },
      { title: "RL scoring", body: "Only registered runnable QSAR targets appear. Combine multiple targets with MW, TPSA, QED and other property scores." }
    ],
    sampling: "Sampling: generate molecules directly from an existing model without updating it.",
    tl: "Transfer Learning: learn structural patterns from known SMILES and create a generator specialized for a target chemical space.",
    rl: "Reinforcement Learning: score generated molecules and update the Agent to increase the probability of high-scoring molecules.",
    tlrl: "TL → RL: learn the target chemical space with TL, then optimize activity or properties with RL. TL and QSAR can be built independently and combined in RL.",
    qsar: "QSAR is a property/activity predictor, not a molecular generator. Only registered and runnable QSAR models become RL scoring candidates.",
    generators: [
      { key: "reinvent", name: "Reinvent", purpose: "Generate new molecules without a parent", research: "Broad chemical-space exploration and de novo design", input: "Count, model and device", parent: "Not required", output: "Diverse molecules sampled from the Prior", example: "Sample 100 molecules from REINVENT.prior" },
      { key: "libinvent", name: "LibInvent", purpose: "Generate R-groups from attachment points on a fixed scaffold", research: "Side-chain optimization and SAR with a retained core", input: "Scaffold SMILES containing attachment point *", parent: "Required (scaffold)", output: "Molecules retaining the scaffold with varied side chains", example: "N#CC[C@H](C1CCCC1)[*:1]" },
      { key: "linkinvent", name: "LinkInvent", purpose: "Generate a linker between two fixed fragments", research: "Fragment linking and PROTAC linker exploration", input: "One * per fragment, separated by |", parent: "Required (two fragments)", output: "Molecules retaining both fragments with varied linkers", example: "Oc1cncc(*)c1|*c1ccoc1" },
      { key: "mol2mol_high", name: "Mol2Mol", purpose: "Generate molecules similar to a parent", research: "Analogue search, matched pairs and scaffold hopping", input: "One parent SMILES", parent: "Required", output: "Transformations within the selected similarity range", example: "CC(=O)OC1=CC=CC=C1C(=O)O" }
    ],
    installationTitle: "Installation & Setup", installationIntro: "Run the included setup script from a clean environment. Chemprop v2 and QSARtuna are installed in isolated runtimes to avoid dependency conflicts.",
    installation: {
      title: "Supported: Apple Silicon macOS / Ubuntu CPU and CUDA where available",
      mac: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_macos.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      ubuntu: ["git clone https://github.com/MolecularAI/REINVENT4.git", "cd REINVENT4", "bash scripts/setup_ubuntu.sh", "python3 scripts/check_environment.py", "bash scripts/start_webapp.sh"],
      verify: ["python3 --version", "python3 scripts/check_environment.py", "curl http://127.0.0.1:8000/api/health", "Browser: http://127.0.0.1:5173/"]
    },
    ossTitle: "Primary open-source software", noRocs: "OpenEye ROCS is not used. Structure and similarity calculations use open-source software such as RDKit."
  }
};
