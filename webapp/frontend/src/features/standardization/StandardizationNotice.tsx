interface Props { version?: string; }

export function StandardizationNotice({ version = "1.0.0" }: Props) {
  return <details className="standardization-notice">
    <summary>共通SMILES標準化条件 · v{version}</summary>
    <ul>
      <li>RDKit MolStandardizeでCleanup / Normalizeを実行</li>
      <li>重原子数が最大のフラグメントを採用（同数時も決定的に選択）</li>
      <li>中和可能な電荷を中和し、永久電荷は保持</li>
      <li>標準化後の分子量が200以下、無効SMILES、CXSMILES区切りを除外</li>
      <li>立体情報を保持したcanonical SMILESを学習・予測・RL評価に使用</li>
    </ul>
    <p>生成元SMILESは変更せず、元SMILESと標準化SMILESを監査用CSVに保存します。</p>
  </details>;
}
