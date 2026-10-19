export type GenerationUploadMetadata = {
  kind: "generation";
  source: "transfer_learning" | "reinforcement_learning";
  modelType: "reinvent" | "mol2mol";
  baseModelKey: string;
  comment: string;
};

export type QsarUploadMetadata = {
  kind: "qsar";
  algorithm: string;
  comment: string;
};

export type UploadMetadata = GenerationUploadMetadata | QsarUploadMetadata;

const generationKeys: Record<string, Omit<GenerationUploadMetadata, "kind" | "source" | "comment">> = {
  reinvent: { modelType: "reinvent", baseModelKey: "reinvent" },
  mol2mol: { modelType: "mol2mol", baseModelKey: "mol2mol_scaffold_generic" },
  mol2molhigh: { modelType: "mol2mol", baseModelKey: "mol2mol_high" },
  mol2molmedium: { modelType: "mol2mol", baseModelKey: "mol2mol_medium" },
  mol2molmmp: { modelType: "mol2mol", baseModelKey: "mol2mol_mmp" },
  mol2molscaffold: { modelType: "mol2mol", baseModelKey: "mol2mol_scaffold" },
  mol2molscaffoldgeneric: { modelType: "mol2mol", baseModelKey: "mol2mol_scaffold_generic" },
  mol2molsimilarity: { modelType: "mol2mol", baseModelKey: "mol2mol_similarity" }
};

const qsarAlgorithms = new Set([
  "randomforest", "svm", "gradientboosting", "extratrees", "lightgbm", "xgboost", "catboost",
  "kernelridge", "gaussianprocess", "ridge", "elasticnet", "logisticregression", "mlp", "chempropv2", "qsartuna"
]);

function token(value: string): string {
  return value.toLowerCase().replace(/[^a-z0-9]/g, "");
}

function algorithmKey(value: string): string {
  const aliases: Record<string, string> = {
    randomforest: "random_forest", gradientboosting: "gradient_boosting", extratrees: "extra_trees",
    kernelridge: "kernel_ridge", gaussianprocess: "gaussian_process", elasticnet: "elastic_net",
    logisticregression: "logistic_regression", chempropv2: "chemprop_v2"
  };
  return aliases[token(value)] ?? token(value);
}

export function inferModelFilename(filename: string): UploadMetadata | null {
  const stem = filename.replace(/\.[^.]+$/, "");
  const parts = stem.split("_");
  if (parts.length < 4 || !/^\d{8}$/.test(parts[0] ?? "")) return null;
  const kind = parts[1]?.toUpperCase();
  const comment = parts.slice(3).join("_");
  if (!comment) return null;
  if (kind === "QSAR" && qsarAlgorithms.has(token(parts[2] ?? ""))) {
    return { kind: "qsar", algorithm: algorithmKey(parts[2] ?? ""), comment };
  }
  if ((kind === "TL" || kind === "RL") && generationKeys[token(parts[2] ?? "")]) {
    return {
      kind: "generation",
      source: kind === "RL" ? "reinforcement_learning" : "transfer_learning",
      ...generationKeys[token(parts[2] ?? "")],
      comment
    };
  }
  return null;
}
