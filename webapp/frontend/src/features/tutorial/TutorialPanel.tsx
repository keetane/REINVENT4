interface Props {
  onUseSample: (modelKey: string, smiles: string) => void;
}

const LIBINVENT_SAMPLE = "N#CC[C@H](C1CCCC1)[*:1]";
const LINKINVENT_SAMPLE = "Oc1cncc(*)c1|*c1ccoc1";

export function TutorialPanel({ onUseSample }: Props) {
  return (
    <div className="tutorial-panel">
      <div className="section-heading tutorial-hero">
        <div>
          <p className="section-kicker">00 · TUTORIAL</p>
          <h2>分子生成をはじめる</h2>
          <p>モデルごとの親分子の書き方から、生成結果の比較・保存までを短く確認できます。</p>
        </div>
        <span className="tutorial-chip">LOCAL WORKSPACE</span>
      </div>

      <section className="tutorial-section">
        <div className="tutorial-section-title"><span>01</span><div><h3>基本の流れ</h3><p>入力から解析まで4ステップです。</p></div></div>
        <ol className="tutorial-steps">
          <li><strong>モデルを選ぶ</strong><span>変換方法に合わせてMol2Mol、LibInvent、LinkInventなどを選択します。</span></li>
          <li><strong>親分子を入力</strong><span>SMILESを貼り付けるか、Ketcherで構造を描画します。</span></li>
          <li><strong>生成を実行</strong><span>生成件数とCPU・GPU・MPSを選び、ジョブを開始します。</span></li>
          <li><strong>結果を解析</strong><span>構造、物性、類似度を絞り込み、CSVまたはSDFで保存します。</span></li>
        </ol>
      </section>

      <section className="tutorial-section">
        <div className="tutorial-section-title"><span>02</span><div><h3>LibInvent</h3><p>scaffoldを固定し、接続点から側鎖を生成します。</p></div></div>
        <div className="tutorial-example">
          <div><span className="example-label">SAMPLE SMILES</span><code>{LIBINVENT_SAMPLE}</code></div>
          <button className="button secondary" type="button" onClick={() => onUseSample("libinvent", LIBINVENT_SAMPLE)}>この例を使う</button>
        </div>
        <p className="tutorial-note"><code>*</code> が側鎖の接続点です。各接続点はscaffoldの原子1個に結合させます。</p>
      </section>

      <section className="tutorial-section">
        <div className="tutorial-section-title"><span>03</span><div><h3>LinkInvent</h3><p>2つのwarheadを固定し、その間をつなぐlinkerを生成します。</p></div></div>
        <div className="tutorial-example">
          <div><span className="example-label">SAMPLE SMILES</span><code>{LINKINVENT_SAMPLE}</code></div>
          <button className="button secondary" type="button" onClick={() => onUseSample("linkinvent", LINKINVENT_SAMPLE)}>この例を使う</button>
        </div>
        <p className="tutorial-note">各warheadに接続点 <code>*</code> を1個置き、2つのSMILESを <code>|</code> で区切ります。</p>
      </section>

      <section className="tutorial-section">
        <div className="tutorial-section-title"><span>04</span><div><h3>解析結果の見方</h3><p>医薬品らしさ、合成容易性、構造多様性をまとめて比較できます。</p></div></div>
        <div className="tutorial-property-grid">
          <div><strong>基本物性</strong><span>MW、LogP、TPSA、HBD、HBA、RotB</span></div>
          <div><strong>品質指標</strong><span>QED、SA Score、Lipinski違反数</span></div>
          <div><strong>構造指標</strong><span>環数、芳香環数、Fsp³、重原子数、形式電荷、MolMR</span></div>
          <div><strong>類似度</strong><span>ParentSimilarity、ScaffoldSimilarity</span></div>
        </div>
      </section>

      <section className="tutorial-section oss-section">
        <div className="tutorial-section-title"><span>05</span><div><h3>使用している主なOSS</h3><p>本アプリは次のオープンソースソフトウェアを利用しています。</p></div></div>
        <div className="oss-grid">
          <div><strong>REINVENT4</strong><span>分子生成 · Apache-2.0</span></div>
          <div><strong>RDKit</strong><span>化学構造・物性計算 · BSD-3-Clause</span></div>
          <div><strong>Ketcher</strong><span>構造描画エディタ · Apache-2.0</span></div>
          <div><strong>FastAPI</strong><span>APIサーバー · MIT</span></div>
          <div><strong>React</strong><span>ユーザーインターフェース · MIT</span></div>
          <div><strong>Vite</strong><span>フロントエンド開発環境 · MIT</span></div>
          <div><strong>SQLite</strong><span>ジョブデータベース · Public Domain</span></div>
        </div>
        <p className="tutorial-note">配布時は同梱される各OSSのライセンス本文とバージョン別の表示を確認してください。</p>
      </section>
    </div>
  );
}
