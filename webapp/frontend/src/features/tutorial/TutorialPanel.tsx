import { useState, type ReactNode } from "react";
import { useTranslation } from "react-i18next";

import { tutorialContent } from "../../i18n/tutorialContent";

interface Props {
  onUseSample: (modelKey: string, smiles: string) => void;
  onOpenTraining: () => void;
}

const OSS = [
  ["REINVENT4", "4.5.11", "Apache-2.0", "https://github.com/MolecularAI/REINVENT4", "https://github.com/MolecularAI/REINVENT4/tree/main/docs"],
  ["RDKit", "2024.9.6", "BSD-3-Clause", "https://github.com/rdkit/rdkit", "https://www.rdkit.org/docs/"],
  ["scikit-learn", "1.6.1", "BSD-3-Clause", "https://github.com/scikit-learn/scikit-learn", "https://scikit-learn.org/stable/"],
  ["Chemprop", "2.3.1", "MIT", "https://github.com/chemprop/chemprop", "https://chemprop.readthedocs.io/"],
  ["QSARtuna", "4.0.1", "Apache-2.0 (pyproject)", "https://github.com/MolecularAI/QSARtuna", "https://github.com/MolecularAI/QSARtuna#readme"],
  ["LightGBM", "4.7.0", "MIT", "https://github.com/microsoft/LightGBM", "https://lightgbm.readthedocs.io/"],
  ["XGBoost", "3.2.0", "Apache-2.0", "https://github.com/dmlc/xgboost", "https://xgboost.readthedocs.io/"],
  ["CatBoost", "1.2.10", "Apache-2.0", "https://github.com/catboost/catboost", "https://catboost.ai/docs/"],
  ["PyTorch", "2.6.0 / 2.14.1", "BSD-3-Clause", "https://github.com/pytorch/pytorch", "https://pytorch.org/docs/stable/"],
];

type TutorialItem = { title: string; body: string };

export function TutorialPanel({ onUseSample, onOpenTraining }: Props) {
  const { i18n } = useTranslation();
  const language = i18n.language.startsWith("en") ? "en" : "ja";
  const content = tutorialContent[language];
  const [selectedSource, setSelectedSource] = useState<"standard" | "transfer_learning" | "reinforcement_learning">("standard");
  const [selectedGenerator, setSelectedGenerator] = useState("reinvent");
  const modelChoice = content.modelChoices.find((item) => item.key === selectedSource) ?? content.modelChoices[0];
  const generator = content.generators.find((item) => item.key === selectedGenerator) ?? content.generators[0];
  const labels = content.labels;

  return <div className="tutorial-panel">
    <div className="section-heading tutorial-hero">
      <div><p className="section-kicker">TUTORIAL</p><h2>{content.title}</h2><p>{content.intro}</p></div>
      <span className="tutorial-chip">LOCAL WORKSPACE</span>
    </div>

    <section className="tutorial-section installation-section">
      <SectionTitle number="00" title={content.environmentTitle} description={content.environmentIntro} />
      <p className="tutorial-note">{content.installation.title}</p>
      <div className="installation-grid">
        <CommandList title={labels.mac} commands={content.installation.mac} copyLabel={labels.copy} copiedLabel={labels.copied} />
        <CommandList title={labels.ubuntu} commands={content.installation.ubuntu} copyLabel={labels.copy} copiedLabel={labels.copied} />
        <CommandList title={labels.verify} commands={content.installation.verify} copyLabel={labels.copy} copiedLabel={labels.copied} />
      </div>
    </section>

    <section className="tutorial-section">
      <SectionTitle number="01" title={content.generationTitle} description={content.generationIntro} />
      <div className="tutorial-model-choice">
        <div className="tutorial-model-picker field-group">
          <label htmlFor="tutorial-model-source">{labels.modelSource}</label>
          <select id="tutorial-model-source" value={selectedSource} onChange={(event) => setSelectedSource(event.target.value as typeof selectedSource)}>
            {content.modelChoices.map((item) => <option key={item.key} value={item.key}>{item.title}</option>)}
          </select>
        </div>
        <article className="tutorial-model-card">
          <div className="tutorial-model-copy"><h4>{modelChoice.title}</h4><p>{modelChoice.body}</p></div>
        </article>
      </div>
      <InfoGrid items={content.generationSteps} />
      <div className="tutorial-sample-block">
        <div className="tutorial-sample-heading">
          <strong>{content.sampleTitle}</strong>
          <select aria-label={labels.standardSample} value={selectedGenerator} onChange={(event) => setSelectedGenerator(event.target.value)}>
            {content.generators.map((item) => <option key={item.key} value={item.key}>{item.name}</option>)}
          </select>
        </div>
        <div className="tutorial-example">
          {generator.example ? <code>{generator.example}</code> : <span className="tutorial-no-parent">{content.noParent}</span>}
          <button className="button secondary" type="button" onClick={() => onUseSample(generator.key, generator.example)}>{labels.use}</button>
        </div>
      </div>
    </section>

    <section className="tutorial-section">
      <SectionTitle number="02" title={content.transferTitle} description={content.transferIntro} />
      <InfoGrid items={content.transferSteps} />
      <Subsection title={content.transferResultsTitle}><InfoGrid items={content.transferResults} /></Subsection>
      <Subsection title={content.transferParametersTitle}><InfoGrid items={content.transferParameters} /></Subsection>
      <div className="tutorial-training-action"><p>{content.transferIntro}</p><button className="button secondary" type="button" onClick={onOpenTraining}>{content.openTraining}</button></div>
    </section>

    <section className="tutorial-section">
      <SectionTitle number="03" title={content.qsarTitle} description={content.qsarIntro} />
      <Workflow title={content.qsarWorkflowTitle} steps={content.qsarWorkflowSteps} />
      <InfoGrid items={content.qsarUses} />
      <Subsection title={content.qsarScoresTitle}><InfoGrid items={content.qsarScores} /></Subsection>
    </section>

    <section className="tutorial-section">
      <SectionTitle number="04" title={content.rlTitle} description={content.rlIntro} />
      <Workflow title={content.rlWorkflowTitle} steps={content.rlWorkflowSteps} />
      <InfoGrid items={content.rlSteps} />
    </section>

    <section className="tutorial-section oss-section">
      <SectionTitle number="05" title={content.ossTitle} />
      <div className="qsar-preview"><table><thead><tr><th>OSS</th><th>Version</th><th>License</th><th>{labels.source}</th><th>{labels.docs}</th></tr></thead><tbody>{OSS.map((item) => <tr key={item[0]}><td>{item[0]}</td><td>{item[1]}</td><td>{item[2]}</td><td><a href={item[3]} target="_blank" rel="noreferrer">GitHub ↗</a></td><td><a href={item[4]} target="_blank" rel="noreferrer">Docs ↗</a></td></tr>)}</tbody></table></div>
      <p className="oss-exclusion-note"><strong>{content.noRocs}</strong></p>
    </section>
  </div>;
}

function SectionTitle({ number, title, description }: { number: string; title: string; description?: string }) {
  return <div className="tutorial-section-title"><span>{number}</span><div><h3>{title}</h3>{description && <p>{description}</p>}</div></div>;
}

function InfoGrid({ items }: { items: TutorialItem[] }) {
  return <div className="tutorial-info-grid">{items.map((item) => <article key={item.title}><strong>{item.title}</strong><span>{item.body}</span></article>)}</div>;
}

function Subsection({ title, children }: { title: string; children: ReactNode }) {
  return <section className="tutorial-subsection"><h4>{title}</h4>{children}</section>;
}

function Workflow({ title, steps }: { title: string; steps: string[] }) {
  return <div className="tutorial-workflow"><h4>{title}</h4><div className="workflow-flow" aria-label={title}>{steps.map((item, index) => <div key={item}><span>{item}</span>{index < steps.length - 1 && <b>→</b>}</div>)}</div></div>;
}

function CommandList({ title, commands, copyLabel, copiedLabel }: { title: string; commands: string[]; copyLabel: string; copiedLabel: string }) {
  const [copied, setCopied] = useState(false);
  const commandText = commands.join("\n");

  const copy = async () => {
    try {
      await navigator.clipboard.writeText(commandText);
      setCopied(true);
      window.setTimeout(() => setCopied(false), 1600);
    } catch {
      setCopied(false);
    }
  };

  return <article><div className="installation-command-heading"><h4>{title}</h4><button className="button quiet installation-copy-button" type="button" onClick={() => void copy()}>{copied ? copiedLabel : copyLabel}</button></div><pre>{commandText}</pre></article>;
}
