# REINVENTer: REINVENT4 統合GUI

`REINVENTer` は、創薬支援ツール [REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5) のためのStreamlit製インタラクティブGUIです。
このアプリケーション一つで、**分子生成**、**転移学習**、**結果の分析**という創薬化学研究における一連のワークフローをシームレスに実行できます。

![メインページ](screenshot/Screenshot%202025-06-26%20at%205.48.55.png)

## 主な機能

REINVENTerは、サイドバーからアクセスできる複数のページで構成されています。

1.  **メインページ (`reinventer.py`)**:
    -   REINVENT4が生成した分子ライブラリ（CSVファイル）を読み込み、インタラクティブに分析します。
    -   分子量(MW)、LogP、TPSAなど、様々な物理化学的特性に基づいて、リアルタイムでのフィルタリングやソートが可能です。
    -   フィルタリング後のデータはCSVとしてダウンロードできます。

2.  **Auto R-group Decomposition (`pages/auto_R_group.py`)**:
    -   親分子（骨格）を入力すると、Recapアルゴリズムで自動的にR-group分解を行います。
    -   分解して得られたフラグメントを元に、以下の3つのモードで新規分子を生成（サンプリング）します。
        -   **Mol2Mol**: 親分子に類似した化合物を生成します。
        -   **LibInvent**: 指定したフラグメント（子分子）からライブラリを構築します。
        -   **LinkInvent**: 複数のフラグメント（弾頭）を連結して新規分子を設計します。
    -   ![Auto R-group Decompositionページ](screenshot/Screenshot%202025-06-26%20at%205.49.04.png)

3.  **Transfer Learning (`pages/Transfer_Learning.py`)**:
    -   既存の化合物情報（ターゲットタンパク質など）を用いて、事前学習済みモデルを追加学習（転移学習）させることができます。
    -   **ChEMBL連携**: UniProt IDを入力するだけで、ChEMBLデータベースから関連化合物を自動で取得し、学習データとして利用できます。
    -   **CSVアップロード**: 独自のSMILESリスト（CSVファイル）をアップロードして学習させることも可能です。
    -   学習後に生成された新しいモデルを使って、特定の化学空間に特化した分子をサンプリングできます。
    -   ![Transfer Learningページ](screenshot/Screenshot%202025-06-26%20at%205.49.15.png)

## インストールと起動

### 必要なライブラリ

本アプリケーションを実行するには、以下のライブラリが必要です。

```bash
pip install streamlit pandas rdkit pubchempy chembl_webresource_client
```

また、`reinvent` コマンドが実行可能な環境（例: `conda` 環境）が設定されている必要があります。

### 起動方法

ターミナルで以下のコマンドを実行して、Streamlitアプリケーションを起動します。

```bash
streamlit run reinventer.py
```

## 使い方

### 0. 作業ディレクトリの設定（全ページ共通）

最初に、サイドバーで作業ディレクトリを設定します。

1.  アプリが起動すると、サイドバーに作業ディレクトリが表示されます。デフォルトは `(カレントディレクトリ)/1_wd/all_users` です。
2.  テキストボックスに任意のユーザー名（またはプロジェクト名）を入力すると、その名前のディレクトリが `1_wd` の下に作成され、結果やモデルが整理されます。

![サイドバー](screenshot/Screenshot%202025-06-26%20at%205.48.10.png)

### 1. 結果の分析 (メインページ)

1.  REINVENT4で生成された結果のCSVファイルを、対応する作業ディレクトリ内の `results` フォルダに配置します。
2.  サイドバーの `Results CSV File Selection` から分析したいCSVファイルを選択します。
3.  メイン画面に分子構造が表示されたら、サイドバーの `Sort by` や `Filtering Parameters` を使って、分子をインタラクティブに絞り込みます。
4.  `Download CSV` ボタンで、フィルタリング後のデータを保存できます。

### 2. R-groupベースの分子生成 (Auto R-group Decompositionページ)

1.  サイドバーから `auto_R_group` ページに移動します。
2.  `Parent Molecule Drawer` で、骨格となる親分子を **化合物名で検索** するか、**SMILESを直接入力** します。
3.  `Children by Recap Decomposition` に、親分子が自動で分解されたフラグメントが表示されます。
4.  `ReInvent Sampling` セクションで、実行したいサンプリングモード（Mol2Mol, LibInvent, LinkInvent）を選択します。
    -   **Mol2Mol**: `Model File` を選択し、`Mol2Mol Sampling from Parent` ボタンをクリックします。
    -   **LibInvent**: `Select Child Molecule` で使用するフラグメントを選択し、`LibInvent Sampling` ボタンをクリックします。
    -   **LinkInvent**: `Select Warhead` で連結したいフラグメントを複数選択し、`LinkInvent Sampling` ボタンをクリックします。
5.  実行が完了すると、`results` ディレクトリに結果のCSVファイルが生成されます。メインページに戻って結果を確認しましょう。

### 3. 転移学習 (Transfer Learningページ)

1.  サイドバーから `Transfer_Learning` ページに移動します。
2.  `Transfer learning from ChEMBL or CSV` セクションで、学習データのソースを選択します。
    -   **from ChEMBL**: `Uniprot ID` を入力し、表示されたリストからターゲットを選択して `Transfer learning from ChEMBL` ボタンをクリックします。
    -   **from csv**: `Upload CSV file` からSMILESを含むCSVファイルをアップロードし、`Transfer learning from csv` ボタンをクリックします。
3.  学習が完了すると、`model` ディレクトリに新しいモデルファイル（`.model`）が保存されます。
4.  `Sampling from Transfer Learning Model` セクションで、今作成したモデルを選択し、`Sample from TL Model` ボタンをクリックすると、そのモデルから新しい分子が生成されます。

## ディレクトリ構造

本ツールは、指定された作業ディレクトリ（例: `1_wd/my_project`）の中に、以下のディレクトリを自動的に作成・使用します。

-   `results/`: 分子生成（サンプリング）の結果（CSVファイル）が格納されます。
-   `input/`: REINVENT4の実行に必要な入力ファイル（SMILESファイルなど）が格納されます。
-   `model/`: 転移学習で生成されたモデルファイルが保存されます。
-   `toml/`: REINVENT4の実行設定ファイル（.toml）が一時的に作成されます。
-   `results/log/`: 実行時のログファイルが保存されます。

## 注意事項

-   `reinvent` コマンドへのパスが通っている必要があります。conda環境などでREINVENT4をインストールした場合は、その環境を有効化してから`streamlit run`を実行してください。
-   macOSで`mps`デバイスを使用すると、転移学習時にエラーが発生することがあります。その場合、自動的にデバイスが`cpu`に変更されます。
-   ChEMBLからのデータ取得にはインターネット接続が必要です。