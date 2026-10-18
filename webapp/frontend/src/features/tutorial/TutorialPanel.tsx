import { useState } from "react";
import { useTranslation } from "react-i18next";

import { tutorialContent } from "../../i18n/tutorialContent";

interface Props { onUseSample: (modelKey: string, smiles: string) => void; onOpenTraining: () => void; }

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

export function TutorialPanel({ onUseSample, onOpenTraining }: Props) {
  const { i18n } = useTranslation();
  const language = i18n.language.startsWith("en") ? "en" : "ja";
  const content = tutorialContent[language];
  const [selected, setSelected] = useState("reinvent");
  const generator = content.generators.find((item) => item.key === selected) ?? content.generators[0];
  const labels = content.labels;
  return <div className="tutorial-panel">
    <div className="section-heading tutorial-hero"><div><p className="section-kicker">00 · TUTORIAL</p><h2>{content.title}</h2><p>{content.intro}</p></div><span className="tutorial-chip">LOCAL WORKSPACE</span></div>
    <section className="tutorial-section"><div className="tutorial-section-title"><span>01</span><div><h3>{content.chooseTitle}</h3></div></div>
      <div className="decision-tree">{content.decisionTree.map((item) => <span key={item}>{item}</span>)}</div>
      <div className="tutorial-model-picker field-group"><label htmlFor="tutorial-model">Model</label><select id="tutorial-model" value={selected} onChange={(event) => setSelected(event.target.value)}>{content.generators.map((item) => <option key={item.key} value={item.key}>{item.name}</option>)}</select></div>
      <article className="tutorial-model-card"><div className="tutorial-model-copy"><span className="example-label">{generator.name}</span><p><strong>{labels.purpose}</strong>{generator.purpose}</p><p><strong>{labels.research}</strong>{generator.research}</p><p><strong>{labels.input}</strong>{generator.input}</p><p><strong>{labels.parent}</strong>{generator.parent}</p><p><strong>{labels.output}</strong>{generator.output}</p></div><div className="tutorial-example">{generator.key === "reinvent" ? <span className="tutorial-no-parent">{content.noParent}</span> : <code>{generator.example}</code>}<button className="button secondary" type="button" onClick={() => onUseSample(generator.key, generator.key === "reinvent" ? "" : generator.example)}>{labels.use}</button></div></article>
    </section>
    <section className="tutorial-section"><div className="tutorial-section-title"><span>02</span><div><h3>Transfer Learning</h3><p>{content.step2Title}</p></div></div><ol className="tutorial-training-steps"><li><strong>Sampling</strong><span>{content.sampling}</span></li><li><strong>Transfer Learning</strong><span>{content.tl}</span></li><li><strong>Reinforcement Learning</strong><span>{content.rl}</span></li><li><strong>TL → RL</strong><span>{content.tlrl}</span></li></ol><button className="button secondary" type="button" onClick={onOpenTraining}>{content.openTraining}</button></section>
    <section className="tutorial-section"><div className="tutorial-section-title"><span>03</span><div><h3>{content.qsarTitle}</h3><p>{content.qsar}</p></div></div><div className="workflow-flow" aria-label={content.workflowTitle}>{content.workflowSteps.map((item, index) => <div key={item}><span>{item}</span>{index < content.workflowSteps.length - 1 && <b>→</b>}</div>)}</div></section>
    <section className="tutorial-section"><div className="tutorial-section-title"><span>04</span><div><h3>QSAR / Model Comparison / RL</h3></div></div><ol className="tutorial-training-steps">{content.modelingSteps.map((item) => <li key={item.title}><strong>{item.title}</strong><span>{item.body}</span></li>)}</ol></section>
    <section className="tutorial-section oss-section"><div className="tutorial-section-title"><span>05</span><div><h3>{content.ossTitle}</h3></div></div><div className="qsar-preview"><table><thead><tr><th>OSS</th><th>Version</th><th>License</th><th>{labels.source}</th><th>{labels.docs}</th></tr></thead><tbody>{OSS.map((item) => <tr key={item[0]}><td>{item[0]}</td><td>{item[1]}</td><td>{item[2]}</td><td><a href={item[3]} target="_blank" rel="noreferrer">GitHub ↗</a></td><td><a href={item[4]} target="_blank" rel="noreferrer">Docs ↗</a></td></tr>)}</tbody></table></div><p className="oss-exclusion-note"><strong>{content.noRocs}</strong></p></section>
    <section className="tutorial-section installation-section"><div className="tutorial-section-title"><span>06</span><div><h3>{content.installationTitle}</h3><p>{content.installationIntro}</p></div></div><p className="tutorial-note">{content.installation.title}</p><div className="installation-grid"><CommandList title={labels.mac} commands={content.installation.mac} /><CommandList title={labels.ubuntu} commands={content.installation.ubuntu} /><CommandList title={labels.verify} commands={content.installation.verify} /></div></section>
  </div>;
}

function CommandList({ title, commands }: { title: string; commands: string[] }) {
  return <article><h4>{title}</h4><pre>{commands.join("\n")}</pre></article>;
}
