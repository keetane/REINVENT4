export interface HealthResponse {
  status: "ok";
  database: "ok";
  version: string;
}

export type JobStatus =
  | "queued"
  | "running"
  | "completed"
  | "failed"
  | "cancelling"
  | "cancelled"
  | "interrupted";

export interface Job {
  id: string;
  status: JobStatus;
  created_at: string;
  queued_at: string;
  started_at: string | null;
  finished_at: string | null;
  updated_at: string;
  comment: string;
  parent_smiles: string;
  parent_smiles_canonical: string | null;
  progress: {
    current: number | null;
    total: number | null;
    percent: number | null;
    message: string | null;
  };
  exit_code: number | null;
  error_summary: string | null;
  job_type: "generation" | "transfer_learning" | "qsar_training" | "qsar_comparison" | "reinforcement_learning";
  model_id: string | null;
  metrics: {
    training_loss?: number | null;
    validation_loss?: number | null;
    training?: Record<string, number | null>;
    test?: Record<string, number | null>;
    cross_validation?: Record<string, number | boolean | string | null>;
    mean_total_score?: number | null;
    fraction_valid_smiles?: number | null;
    components?: Record<string, number | null>;
    history?: Record<string, unknown>[];
  };
}

export interface QsarCapability {
  key: string;
  label: string;
  available: boolean;
  reason: string | null;
  version: string | null;
  task_types: ("regression" | "classification")[];
}

export interface QsarCapabilities {
  algorithms: QsarCapability[];
  devices: QsarCapability[];
}

export interface QsarDataInspection {
  columns: string[];
  total_rows: number;
  valid_count: number;
  invalid_smiles_count: number;
  missing_target_count: number;
  invalid_target_count: number;
  duplicate_count: number;
  conflicting_duplicate_count: number;
  standardized_count: number;
  fragment_extracted_count: number;
  neutralized_count: number;
  charged_structure_retained_count: number;
  low_molecular_weight_excluded_count: number;
  final_accepted_count: number;
  standardization_version: string;
  preview: Record<string, string>[];
  rejected_rows: { row_number: number; smiles: string; target: string; reason: string; standardized_smiles: string | null }[];
}

export interface QsarModel {
  id: string;
  name: string;
  algorithm: string;
  task_type: "regression" | "classification";
  target_name: string;
  target_unit: string;
  created_at: string;
  training_data_count: number;
  test_data_count: number;
  metrics: Record<string, unknown>;
  metadata: Record<string, unknown>;
  ad_threshold: number;
  in_use: boolean;
}

export interface QsarDataPayload {
  file_name: string;
  content: string;
  smiles_column: string | null;
  target_column: string | null;
  task_type: "regression" | "classification";
  target_transform: "none" | "p_activity";
  concentration_unit: "M" | "mM" | "uM" | "μM" | "nM" | "pM";
  activity_type: "IC50" | "EC50" | "Ki" | "Kd" | "custom";
  duplicate_strategy: "exclude_conflicts" | "keep_first" | "mean";
}

export interface ModelInfo {
  id: string;
  name: string;
  model_type: string;
  source: "standard" | "transfer_learning" | "reinforcement_learning";
  description: string;
  created_at: string | null;
  training_data_count: number | null;
  validation_data_count: number | null;
  training_conditions: Record<string, unknown>;
  base_model_key: string | null;
  deletable: boolean;
}

export interface TrainingDataInspection {
  columns: string[];
  total_rows: number;
  valid_count: number;
  invalid_count: number;
  duplicate_count: number;
  standardized_count: number;
  fragment_extracted_count: number;
  neutralized_count: number;
  charged_structure_retained_count: number;
  low_molecular_weight_excluded_count: number;
  model_compatibility_error_count: number;
  standardization_version: string;
}

export interface JobListResponse {
  items: Job[];
  total: number;
  limit: number;
  offset: number;
}

export interface ParentValidationResponse {
  input_smiles: string;
  canonical_smiles: string;
  molblock: string;
  atom_count: number;
  molecular_formula: string;
}

export interface PubChemLookupResponse {
  cid: number;
  molecular_name: string;
  isomeric_smiles: string;
}

export interface JobLogResponse {
  job_id: string;
  lines: string[];
  truncated: boolean;
}

export interface DescriptorSchema {
  name: string;
  source: string;
  aliases: string[];
  minimum: number | null;
  maximum: number | null;
  missing_count: number;
}

export interface ResultSchemaResponse {
  job_id: string;
  columns: string[];
  descriptors: DescriptorSchema[];
  total_rows: number;
  valid_rows: number;
  invalid_smiles: number;
  warnings: string[];
  rdkit_version: string;
}

export interface ResultMolecule {
  row_number: number;
  smiles: string;
  canonical_smiles: string;
  structure_svg: string;
  values: Record<string, number | null>;
}

export interface MoleculeListResponse {
  job_id: string;
  items: ResultMolecule[];
  total: number;
  limit: number;
  offset: number;
  invalid_smiles: number;
}

export interface NumericSummary {
  count: number;
  missing_count: number;
  minimum: number | null;
  maximum: number | null;
  mean: number | null;
  median: number | null;
  standard_deviation: number | null;
  q1: number | null;
  q3: number | null;
}

export interface StatisticsResponse {
  job_id: string;
  total_rows: number;
  filtered_rows: number;
  invalid_smiles: number;
  summaries: Record<string, NumericSummary>;
  histograms: Record<string, { start: number; end: number; count: number }[]>;
  scatter_x: string | null;
  scatter_y: string | null;
  scatter_sampled: boolean;
  scatter: { row_number: number; x: number; y: number }[];
}

export type ResultQuery = Record<string, string | number | null | undefined>;

interface ErrorEnvelope {
  error?: { code?: string; message?: string };
}

export class ApiError extends Error {
  constructor(
    message: string,
    public readonly status: number,
    public readonly code: string
  ) {
    super(message);
  }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(path, {
    ...init,
    headers: {
      Accept: "application/json",
      ...(init?.body ? { "Content-Type": "application/json" } : {}),
      ...init?.headers
    }
  });
  if (!response.ok) {
    let body: ErrorEnvelope = {};
    try {
      body = (await response.json()) as ErrorEnvelope;
    } catch {
      // Keep the fallback uniform instead of exposing raw server text.
    }
    throw new ApiError(
      body.error?.message ?? "The server could not process the request.",
      response.status,
      body.error?.code ?? "request_failed"
    );
  }
  return (await response.json()) as T;
}

export function getHealth(signal?: AbortSignal): Promise<HealthResponse> {
  return request<HealthResponse>("/api/health", { signal });
}

export function validateParent(
  smiles: string,
  modelKey: string
): Promise<ParentValidationResponse> {
  return request<ParentValidationResponse>("/api/parents/validate", {
    method: "POST",
    body: JSON.stringify({ smiles, model_key: modelKey })
  });
}

export function lookupPubChem(name: string): Promise<PubChemLookupResponse> {
  const parameters = new URLSearchParams({ name });
  return request<PubChemLookupResponse>(`/api/parents/pubchem?${parameters.toString()}`);
}

export function createJob(payload: {
  parent_smiles: string;
  comment: string;
  run_conditions: Record<string, unknown>;
}): Promise<Job> {
  return request<Job>("/api/jobs", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function listModels(signal?: AbortSignal): Promise<{ items: ModelInfo[] }> {
  return request<{ items: ModelInfo[] }>("/api/models", { signal });
}

export function deleteModel(modelId: string): Promise<{ id: string; deleted: boolean }> {
  return request<{ id: string; deleted: boolean }>(`/api/models/${encodeURIComponent(modelId)}`, {
    method: "DELETE"
  });
}

export function modelArtifactUrl(modelId: string, kind: "model"): string {
  return `/api/models/${encodeURIComponent(modelId)}/artifacts/${kind}`;
}

export function uploadModel(file: File, metadata: {
  name: string;
  modelType: "reinvent" | "mol2mol";
  baseModelKey: string;
  source: "transfer_learning" | "reinforcement_learning";
}): Promise<ModelInfo> {
  const parameters = new URLSearchParams({
    name: metadata.name,
    model_type: metadata.modelType,
    base_model_key: metadata.baseModelKey,
    source: metadata.source
  });
  return request<ModelInfo>(`/api/models/upload?${parameters.toString()}`, {
    method: "POST",
    body: file,
    headers: {
      "Content-Type": "application/octet-stream",
      "X-Model-Filename": file.name
    }
  });
}

export function inspectTrainingData(payload: {
  file_name: string;
  content: string;
  smiles_column?: string | null;
  base_model?: "reinvent" | "mol2mol";
}): Promise<TrainingDataInspection> {
  return request<TrainingDataInspection>("/api/training-data/inspect", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function createTrainingJob(payload: {
  file_name: string;
  content: string;
  smiles_column?: string | null;
  base_model: "reinvent" | "mol2mol";
  comment: string;
  train_fraction: number;
  num_epochs: number;
  batch_size: number;
  learning_rate: number;
  tanimoto_lower: number;
  tanimoto_upper: number;
  device: "cpu" | "cuda" | "mps";
}): Promise<Job> {
  return request<Job>("/api/training-jobs", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function getQsarCapabilities(): Promise<QsarCapabilities> {
  return request<QsarCapabilities>("/api/qsar/capabilities");
}

export function inspectQsarData(payload: QsarDataPayload): Promise<QsarDataInspection> {
  return request<QsarDataInspection>("/api/qsar-data/inspect", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function createQsarJob(payload: QsarDataPayload & Record<string, unknown>): Promise<Job> {
  return request<Job>("/api/qsar-jobs", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export interface QsarComparisonResult {
  primary_metric: string;
  recommended_algorithm: string;
  selection_basis: string;
  successful_count: number;
  failed_count: number;
  results: {
    algorithm: string;
    status: "completed" | "failed";
    duration_seconds: number;
    metrics?: { training?: Record<string, number | null>; test?: Record<string, number | null>; cross_validation?: Record<string, number | boolean | string | null | number[]>; performance_gap?: Record<string, number>; applicability_domain_fraction?: number; inference_molecules_per_second?: number };
    plots?: Record<string, number[][]>;
    error?: string;
  }[];
}

export function createQsarComparisonJob(payload: QsarDataPayload & Record<string, unknown>): Promise<Job> {
  return request<Job>("/api/qsar-comparison-jobs", { method: "POST", body: JSON.stringify(payload) });
}

export function getQsarComparison(jobId: string): Promise<QsarComparisonResult> {
  return request<QsarComparisonResult>(`/api/qsar-comparisons/${encodeURIComponent(jobId)}`);
}

export function registerQsarComparisonModel(jobId: string, algorithm: string, comment: string): Promise<QsarModel> {
  return request<QsarModel>(`/api/qsar-comparisons/${encodeURIComponent(jobId)}/register`, {
    method: "POST",
    body: JSON.stringify({ algorithm, comment })
  });
}

export function listQsarModels(): Promise<{ items: QsarModel[] }> {
  return request<{ items: QsarModel[] }>("/api/qsar-models");
}

export function deleteQsarModel(modelId: string): Promise<{ id: string; deleted: boolean }> {
  return request<{ id: string; deleted: boolean }>(`/api/qsar-models/${encodeURIComponent(modelId)}`, {
    method: "DELETE"
  });
}

export function qsarArtifactUrl(modelId: string, kind: "model" | "dataset" | "predictions" | "evaluation"): string {
  return `/api/qsar-models/${encodeURIComponent(modelId)}/artifacts/${kind}`;
}

export function getQsarEvaluation(modelId: string): Promise<{ metrics: Record<string, unknown>; metadata: Record<string, unknown>; plots?: Record<string, number[][]> }> {
  return request<{ metrics: Record<string, unknown>; metadata: Record<string, unknown>; plots?: Record<string, number[][]> }>(qsarArtifactUrl(modelId, "evaluation"));
}

export interface RlQsarComponent {
  model_id: string;
  weight: number;
  transform: "sigmoid" | "reverse_sigmoid" | "double_sigmoid" | "none";
  low: number;
  high: number;
  k: number;
  class_index: number;
  ad_behavior: "penalize" | "allow";
}

export interface RlPropertyComponent {
  component: string;
  weight: number;
  transform: "sigmoid" | "reverse_sigmoid" | "double_sigmoid" | "none";
  low: number;
  high: number;
  k: number;
  reference_smiles: string[];
  smarts: string[];
}

export interface RlPreviewRow {
  smiles: string;
  generated_smiles: string;
  standardized_smiles: string | null;
  molecular_weight: number | null;
  formal_charge: number | null;
  exclusion_reason: string | null;
  valid: boolean;
  components: { kind: "qsar" | "property"; name: string; prediction: number | null; in_ad: boolean | null; score: number; weight: number }[];
  total_score: number;
}

export function createRlJob(payload: Record<string, unknown>): Promise<Job> {
  return request<Job>("/api/rl-jobs", { method: "POST", body: JSON.stringify(payload) });
}

export function getRlCapabilities(): Promise<{ devices: { key: "cpu" | "cuda" | "mps"; label: string; available: boolean }[] }> {
  return request<{ devices: { key: "cpu" | "cuda" | "mps"; label: string; available: boolean }[] }>("/api/rl/capabilities");
}

export interface QsarPredictionRow {
  smiles: string;
  standardized_smiles: string | null;
  molecular_weight: number | null;
  formal_charge: number | null;
  valid: boolean;
  prediction: number | null;
  scoring_prediction: number;
  nearest_tanimoto: number | null;
  in_ad: boolean;
  error: string | null;
}

export function predictQsar(modelId: string, smiles: string[]): Promise<{ model_id: string; predictions: QsarPredictionRow[] }> {
  return request<{ model_id: string; predictions: QsarPredictionRow[] }>("/api/qsar/predict", {
    method: "POST",
    body: JSON.stringify({ model_id: modelId, smiles })
  });
}

export function previewRlScoring(payload: {
  smiles: string[];
  qsar_components: RlQsarComponent[];
  property_components: RlPropertyComponent[];
  aggregation: "geometric_mean" | "arithmetic_mean";
}): Promise<{ rows: RlPreviewRow[] }> {
  return request<{ rows: RlPreviewRow[] }>("/api/rl/preview", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function listJobs(signal?: AbortSignal): Promise<JobListResponse> {
  return request<JobListResponse>("/api/jobs?limit=100&offset=0", { signal });
}

export function getJob(jobId: string, signal?: AbortSignal): Promise<Job> {
  return request<Job>(`/api/jobs/${encodeURIComponent(jobId)}`, { signal });
}

export function jobStandardizationArtifactUrl(jobId: string, kind: "accepted" | "rejected"): string {
  return `/api/jobs/${encodeURIComponent(jobId)}/standardization/${kind}`;
}

export function getJobLogs(jobId: string, signal?: AbortSignal): Promise<JobLogResponse> {
  return request<JobLogResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/logs?tail=200`,
    { signal }
  );
}

export function cancelJob(jobId: string): Promise<Job> {
  return request<Job>(`/api/jobs/${encodeURIComponent(jobId)}/cancel`, {
    method: "POST"
  });
}

export function resultQueryString(query: ResultQuery): string {
  const parameters = new URLSearchParams();
  Object.entries(query).forEach(([key, value]) => {
    if (value !== null && value !== undefined && value !== "") {
      parameters.set(key, String(value));
    }
  });
  const encoded = parameters.toString();
  return encoded ? `?${encoded}` : "";
}

export function getResultSchema(jobId: string, signal?: AbortSignal): Promise<ResultSchemaResponse> {
  return request<ResultSchemaResponse>(`/api/jobs/${encodeURIComponent(jobId)}/results/schema`, { signal });
}

export function getResultMolecules(
  jobId: string,
  query: ResultQuery,
  signal?: AbortSignal
): Promise<MoleculeListResponse> {
  return request<MoleculeListResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/results/molecules${resultQueryString(query)}`,
    { signal }
  );
}

export function getResultStatistics(
  jobId: string,
  query: ResultQuery,
  signal?: AbortSignal
): Promise<StatisticsResponse> {
  return request<StatisticsResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/results/statistics${resultQueryString(query)}`,
    { signal }
  );
}

export function resultDownloadUrl(jobId: string, format: "csv" | "sdf", query: ResultQuery): string {
  return `/api/jobs/${encodeURIComponent(jobId)}/results/download.${format}${resultQueryString(query)}`;
}
