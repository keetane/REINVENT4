import os
import subprocess
import pandas as pd
from rdkit.Chem import PandasTools, Draw
from rdkit import Chem
import streamlit as st
import datetime # この行だけを残す

# Get the current timestamp
now = datetime.datetime.now() # このまま
# Format the timestamp as a string
time = now.strftime("%Y%m%d_%H%M")

st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon=":pill:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Set up directories relative to the script's location
wd = st.session_state.wd  # Use the working directory from session state
os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
reinvent_dir = os.getcwd()  # Adjusted REINVENT4 directory path
input_dir = os.path.join(wd, 'input')  # Input directory path
os.makedirs(input_dir, exist_ok=True)  # Create input directory if it doesn't exist
model_dir = os.path.join(wd, 'model')  # Model directory path
os.makedirs(model_dir, exist_ok=True)  # Create model directory if it doesn't exist
results_dir = os.path.join(wd, 'results')  # Results directory path
os.makedirs(results_dir, exist_ok=True)  # Create results directory if it doesn't exist
toml_dir = os.path.join(wd, 'toml')  # TOML directory path
os.makedirs(toml_dir, exist_ok=True)  # Create TOML directory if it doesn't exist
sampling_log = os.path.join(results_dir, 'log')  # Sampling log directory path
os.makedirs(sampling_log, exist_ok=True)  # Create sampling log directory if it doesn't exist

# Set up file paths for priors
priors_dir = os.path.join(reinvent_dir, "priors")
reinvent = os.path.join(priors_dir, "reinvent.prior")
lib = os.path.join(priors_dir, "libinvent.prior")
link = os.path.join(priors_dir, "linkinvent.prior")
mol2mol_high = os.path.join(priors_dir, "mol2mol_high_similarity.prior")
mol2mol_med = os.path.join(priors_dir, "mol2mol_medium_similarity.prior")
mol2mol_mmp = os.path.join(priors_dir, "mol2mol_mmp.prior")
mol2mol_scaffold_generic = os.path.join(priors_dir, "mol2mol_scaffold_generic.prior")
mol2mol_scaffold = os.path.join(priors_dir, "mol2mol_scaffold.prior")
mol2mol_similarity = os.path.join(priors_dir, "mol2mol_similarity.prior")
pubchem = os.path.join(priors_dir, "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior")

# # ここにサイドバー設定のコードを直接記述
# st.sidebar.title("working directory")

# # st.session_state.wd の初期化を確実に行う
# # Streamlitのマルチページアプリでは、メインのファイルでsession_stateを初期化し、
# # 各ページでそれを使うのが一般的です。
# # もし 'wd' がまだセッションステートにない場合は、デフォルト値を設定します。
# if 'wd' not in st.session_state:
#     st.session_state.wd = os.path.expanduser("~/Documents/infomatics/reinvent_data") # デフォルト値を設定

# st.sidebar.text(st.session_state.wd)

# # Set up directories relative to the script's location
# wd = st.session_state.wd  # Use the working directory from session state
# wd = os.path.abspath(os.path.expanduser(wd))
# os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
# # script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the script's directory
# home_dir = os.path.expanduser("~")  # Get the user's home directory
# # REINVENT4ディレクトリのパスを修正。絶対パスにする場合はスラッシュを最初につける
# # 相対パス（home_dirからの相対）にするならスラッシュは不要
# reinvent_dir = os.path.join(home_dir, "Documents", "apps", "REINVENT4")
# input_dir = os.path.join(wd, 'input')  # Input directory path
# os.makedirs(input_dir, exist_ok=True)  # Create input directory if it doesn't exist
# model_dir = os.path.join(wd, 'model')  # Model directory path
# os.makedirs(model_dir, exist_ok=True)  # Create model directory if it doesn't exist
# results_dir = os.path.join(wd, 'results')  # Results directory path
# os.makedirs(results_dir, exist_ok=True)  # Create results directory if it doesn't exist
# toml_dir = os.path.join(wd, 'toml')  # TOML directory path
# os.makedirs(toml_dir, exist_ok=True)  # Create TOML directory if it doesn't exist
# sampling_log = os.path.join(results_dir, 'log')  # Sampling log directory path
# os.makedirs(sampling_log, exist_ok=True)  # Create sampling log directory if it doesn't exist

# # Set up file paths for priors
# priors_dir = os.path.join(reinvent_dir, "priors")
# reinvent = os.path.join(priors_dir, "reinvent.prior")
# lib = os.path.join(priors_dir, "libinvent.prior")
# link = os.path.join(priors_dir, "linkinvent.prior")
# mol2mol_high = os.path.join(priors_dir, "mol2mol_high_similarity.prior")
# mol2mol_med = os.path.join(priors_dir, "mol2mol_medium_similarity.prior")
# mol2mol_mmp = os.path.join(priors_dir, "mol2mol_mmp.prior")
# mol2mol_scaffold_generic = os.path.join(priors_dir, "mol2mol_scaffold_generic.prior")
# mol2mol_scaffold = os.path.join(priors_dir, "mol2mol_scaffold.prior")
# mol2mol_similarity = os.path.join(priors_dir, "mol2mol_similarity.prior")
# pubchem = os.path.join(priors_dir, "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior")

# Sidebar for input parameters
st.sidebar.header("Sampling Parameters")
num_mols = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device = st.sidebar.selectbox("Device", ["cpu", "cuda", "mps"])  # "cpu"を先頭に

unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)
overwrite = st.sidebar.checkbox("Overwrite", value=True)

# Allow users to select a model file from priors_dir
model_files = {
    "High Similarity": mol2mol_high,
    "Medium Similarity": mol2mol_med,
    "MMP": mol2mol_mmp,
    "Scaffold Generic": mol2mol_scaffold_generic,
    "Scaffold": mol2mol_scaffold,
    "Similarity": mol2mol_similarity,
    "PubChem": pubchem,
    "Reinvent": reinvent,
}
# model_file = [] #st.sidebar.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
selected_model_file = model_files['Reinvent'] # Reinventをデフォルト値として設定

st.sidebar.header('mol2mol options')
# other options for mol2mol sampling
sample_stategies = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", 'multinomial'])
temperatures = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)

# parameters for TL in sidebar
st.sidebar.header("Transfer Learning Parameters")
# Epochs
num_epochs = st.sidebar.slider("Number of Epochs", min_value=1, max_value=100, value=3)
save_every_n_epochs = st.sidebar.slider("Save Every N Epochs", min_value=1, max_value=num_epochs, value=3)

# Batch sizes
batch_size = st.sidebar.number_input("Batch Size", min_value=1, max_value=1024, value=50)
sample_batch_size = st.sidebar.number_input("Sample Batch Size", min_value=1, max_value=2048, value=100)

# Number of reference molecules
num_refs = st.sidebar.number_input("Number of Reference Molecules", min_value=0, max_value=2000, value=100)

# Similarity settings
similarity_types = ["tanimoto", "cosine", "dice", "euclidean"]
pairs_type = st.sidebar.selectbox("Similarity Type", options=similarity_types, index=0)

pairs_upper_threshold = st.sidebar.slider("Upper Similarity Threshold", min_value=0.0, max_value=1.0, value=1.0)
pairs_lower_threshold = st.sidebar.slider("Lower Similarity Threshold", min_value=0.0, max_value=1.0, value=0.7)

pairs_min_cardinality = st.sidebar.number_input("Min Cardinality", min_value=1, max_value=1000, value=1)
pairs_max_cardinality = st.sidebar.number_input("Max Cardinality", min_value=pairs_min_cardinality, max_value=2000, value=199)

# 構成データを辞書形式で定義
config = {
    "num_epochs": num_epochs,
    "save_every_n_epochs": save_every_n_epochs,
    "batch_size": batch_size,
    "sample_batch_size": sample_batch_size,
    "num_refs": num_refs,
    "pairs": {
        "type": pairs_type,
        "upper_threshold": pairs_upper_threshold,
        "lower_threshold": pairs_lower_threshold,
        "min_cardinality": pairs_min_cardinality,
        "max_cardinality": pairs_max_cardinality
    }
}

# --- 関数定義をここに移動 ---

# run_reinvent 関数は他の関数から呼ばれるので、一番上に定義する
def run_reinvent(toml_path):
    log_file = os.path.join(sampling_log, "sampling.log")
    # reinvent実行ファイルのパスを正しく指定。MiniforgeやConda環境の場合、binディレクトリにあることが多い
    # `subprocess.run` を使ってエラーコードをチェックする方が堅牢
    try:
        result = subprocess.run(
            ["reinvent", "-l", log_file, toml_path],
            check=True,  # 0以外の終了コードでCalledProcessErrorを発生させる
            capture_output=True, # 標準出力と標準エラー出力をキャプチャ
            text=True # テキストモードでキャプチャ
        )
        st.success(f"REINVENT4 executed successfully for {os.path.basename(toml_path)}.")
        if result.stdout:
            st.code(result.stdout)
        if result.stderr:
            st.warning(f"REINVENT4 stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        st.error(f"REINVENT4 execution failed for {os.path.basename(toml_path)} with error code {e.returncode}.")
        st.error(f"STDOUT: {e.stdout}")
        st.error(f"STDERR: {e.stderr}")
    except FileNotFoundError:
        st.error(f"Error: reinvent executable not found at {home_dir}/miniforge3/envs/r4/bin/reinvent. Please check the path.")
    except Exception as e:
        st.error(f"An unexpected error occurred during REINVENT4 execution: {e}")


# generate toml file for transfer learning
def TL(TL_input, TL_val, TL_method): # TL_methodを引数として受け取るように修正
    TL_input_filename = TL_input.split("/")[-1].replace(".smi", "")
    TL_method_basename = os.path.basename(TL_method).replace(".prior", "").replace(".model", "")
    output_model_file = f"{model_dir}/TL_{TL_input_filename}_{TL_method_basename}.model"
    # deviceを"mps"選択時は"cpu"に強制する（PyTorch capturableエラー対策）
    safe_device = device
    if device == "mps":
        st.warning("PyTorchのcapturable=Trueエラー回避のため、デバイスを'cpu'に変更します。")
        safe_device = "cpu"
    TL_toml_content=f'''
run_type = "transfer_learning"
device = "{safe_device}"  # set torch device e.g. "cpu". For macOS, use "mps"
tb_logdir = "{sampling_log}/tb_TL"  # name of the TensorBoard logging directory
json_out_config = "{sampling_log}/json_transfer_learning.json"  # write this TOML to JSON

[parameters]
num_epochs = {num_epochs}  # number of steps to run
save_every_n_epochs = {save_every_n_epochs}  # save checkpoint model file very N steps
batch_size = {batch_size}  # batch size for training
num_refs = {num_refs}  # number of reference molecules randomly chosen for similarity
                    # set this to zero for large datasets (>200 molecules)!
sample_batch_size = {sample_batch_size}  # number of sampled molecules to compute sample loss
# Uncomment one of the comment blocks below.  Each generator needs a model
# file and possibly a SMILES file with seed structures.

input_model_file = "{TL_method}"
smiles_file = "{TL_input}"  # read 1st column
validation_smiles_file = "{TL_val}"  # read 1st column
output_model_file = '{output_model_file}'  # sampled SMILES and NLL in CSV format

# Define the type of similarity and its parameters
pairs.type = '{pairs_type}'  # e.g. "tanimoto", "cosine", "dice", "euclidean"
pairs.upper_threshold = {pairs_upper_threshold}  # upper threshold for similarity
pairs.lower_threshold = {pairs_lower_threshold}  # lower threshold for similarity
pairs.min_cardinality = {pairs_min_cardinality}  # minimum number of similar molecules
pairs.max_cardinality = {pairs_max_cardinality}  # maximum number of similar molecules

'''

    toml_path = f"{toml_dir}/TL.toml"

    with open(toml_path, "w") as f:
        f.write(TL_toml_content) # 変数名を修正
    return toml_path

# Generate TOML file for sampling from TL model (これはTL_tomlとは異なる機能のようです)
def TL_toml_sampling(method, num_smiles, smiles_file, device, show, sample_strategy, temperature): # 引数を明示的に受け取る
    # ファイル名の衝突を避けるため、datetimeを使う場合はここで定義し直す
    current_time_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    filename = os.path.basename(method).replace(".model", "").replace(".prior", "").replace(" ", "_")
        
    if not overwrite: # overwriteがFalseの場合のみタイムスタンプを追加
        filename = f"{filename}_{current_time_str}"
    
    output_file = os.path.join(results_dir, f"{filename}.csv")
    toml_content = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "{sampling_log}/_TL_sampling.json"

[parameters]
model_file = "{method}"
output_file = "{output_file}"
num_smiles = {num_smiles}
"""
    if unique_molecules: # True/Falseの比較は不要
        toml_content += f'''
        unique_molecules = true
        '''
    if randomize_smiles:
        toml_content += f'''
        randomize_smiles = true
        '''

    if sample_strategy == 'beamsearch':
        toml_content += f'''
        sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
        '''
    elif sample_strategy == 'multinomial':
        toml_content += f'''
        sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
        temperature = {temperature} # temperature in multinomial sampling
        '''

    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml_content)
    return toml_path, output_file

# Generate TOML file (これは汎用的なサンプリング用)
def generate_toml(method, num_smiles, smiles_file, device, show, sample_strategy, temperature): # 引数を明示的に受け取る
    current_time_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    if method == link:
        filename = 'LinkInvent'
    elif method == lib:
        filename = 'LibInvent'
    else:
        filename = os.path.basename(method).replace(".model", "").replace(".prior", "").replace(" ", "_")
        
    if not overwrite: # overwriteがFalseの場合のみタイムスタンプを追加
        filename = f"{filename}_{current_time_str}"
    
    output_file = os.path.join(results_dir, f"{filename}.csv")
    toml_content = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "{sampling_log}/_sampling.json"

[parameters]
model_file = "{method}"
output_file = "{output_file}"
num_smiles = {num_smiles}
"""
    if unique_molecules:
        toml_content += f'''
        unique_molecules = true
        '''
    if randomize_smiles:
        toml_content += f'''
        randomize_smiles = true
        '''
    if smiles_file is not None:
        toml_content += f'''
        smiles_file = "{smiles_file}"  # 1 compound per line
        '''
    
    if sample_strategy == 'beamsearch':
        toml_content += f'''
        sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
        '''
    elif sample_strategy == 'multinomial':
        toml_content += f'''
        sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
        temperature = {temperature} # temperature in multinomial sampling
        '''

    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml_content)
    return toml_path, output_file

# --- 関数定義の終わり ---

# Transfer learning from ChEMBL or CSV
from chembl_webresource_client.new_client import new_client
import requests
from rdkit.Chem import Descriptors

st.header("Transfer learning from ChEMBL or CSV")
TL_model_file_key = "tl_model_file_select"
TL_model_file = st.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x, key=TL_model_file_key)
selected_TL_model_file = model_files[TL_model_file]


col1, col2 = st.columns([3, 1])
with col1:
    st.markdown('#### from ChEMBL')
    uniprot_id = st.text_input("Enter Uniprot ID", value="P10415") #BCL2

    # Always fetch ChEMBL targets for the current Uniprot ID (not only on button press)
    try:
        target_records = list(new_client.target.get(
            target_components__accession=uniprot_id
        ).only(
            'target_chembl_id',
            'organism',
            'pref_name',
            'target_type'
        ))
        df = pd.DataFrame.from_records(target_records)
    except Exception as e:
        st.error(f"Error fetching ChEMBL data for {uniprot_id}: {e}")
        df = pd.DataFrame()

    if not df.empty and 'target_chembl_id' in df.columns:
        st.dataframe(df, use_container_width=True)
        selected_target = st.selectbox(
            "Select ChEMBL target",
            options=df['target_chembl_id'].tolist()
        )
    else:
        st.warning("No valid ChEMBL targets found for the given Uniprot ID.")
        selected_target = None

    if st.button('Transfer learning from ChEMBL'):
        if selected_target is not None:
            # 1. 活性のある化合物のデータ取得
            activities = new_client.activity.filter(
                target_chembl_id__in=selected_target,
                pchembl_value__gte=5,
                assay_type='B',
            ).only([
                'molecule_chembl_id',
                'molecule_pref_name',
                'target_pref_name',
                'parent_molecule_chembl_id',
                'pchembl_value',
                'canonical_smiles',
                'document_chembl_id',
                'document_journal',
            ])
            df = pd.DataFrame.from_records(activities)
            # st.dataframe(df, use_container_width=True)
            # 2. 主成分（分子量最大の構成要素）を抽出する関数
            def extract_main_component(smiles):
                try:
                    fragments = smiles.split('.')
                    max_mol = None
                    max_weight = -1
                    for frag in fragments:
                        mol = Chem.MolFromSmiles(frag)
                        if mol is None:
                            continue
                        mw = Descriptors.MolWt(mol)
                        if mw > max_weight:
                            max_mol = mol
                            max_weight = mw
                    if max_mol:
                        return Chem.MolToSmiles(max_mol, isomericSmiles=False)
                    else:
                        return None
                except Exception as e:
                    print(f"Error splitting SMILES: {smiles} -> {e}")
                    return None

            # 3. 分子量最大部分のSMILES抽出 → 標準化処理
            df["smiles"] = df["canonical_smiles"].apply(extract_main_component)

            # 4. 無効なエントリを削除
            df_clean = df.dropna(subset=["smiles"])

            # 5. 学習用とvalidation用のデータセットを分割
            # ここでは、全体の80%を学習用、20%を検証用に分割
            df_clean = df_clean.sample(frac=1, random_state=42)  # シャッフル
            train_size = int(0.8 * len(df_clean))
            train_df = df_clean.iloc[:train_size]
            val_df = df_clean.iloc[train_size:]
            # st.dataframe(train_df, use_container_width=True)
            # 6. トレーニングデータセットのSMILESを保存
            input_file = input_dir + f'/{selected_target}_tl.smi'
            train_df[['smiles']].to_csv(input_file, sep='\t', index=False, header=False)
            # 7. バリデーションデータセットのSMILESを保存
            val_file = input_dir + f'/{selected_target}_val.smi'
            val_df[['smiles']].to_csv(val_file, sep='\t', index=False, header=False)
            # 8. csvファイルを保存
            df_clean.to_csv(input_dir + f'/{selected_target}.csv', sep='\t', index=False)
            df_clean[['smiles']].to_csv(input_dir + f'/{selected_target}.smi', sep='\t', index=False)
            st.text(f"Saved {len(df_clean)} ligands with largest fragment SMILES to input/{selected_target}.csv and input/{selected_target}.smi")

            # # 5. 保存
            # # ヘッダーなしで保存
            # df_clean[['smiles']].to_csv(input_dir + f'/{selected_target}.smi', sep='\t', index=False, header=False)

            # 6. TOMLファイルを生成
            toml_path = TL(
                input_file,  # 学習用SMILESファイル
                val_file,    # バリデーション用SMILESファイル
                TL_method=selected_TL_model_file # selected_TL_model_fileを渡す
            )

            # 7. REINVENT4を実行
            run_reinvent(toml_path)
            st.success("Transfer learning from ChEMBL completed!")
with col2:
    st.markdown('#### from csv')
    uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file, sep='\t')
            if 'smiles' in df.columns:
                if st.button('Transfer learning from csv'):
                    # Save uploaded file to input_dir and use its path for TL
                    csv_save_path = os.path.join(input_dir, uploaded_file.name)
                    with open(csv_save_path, "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    # Split into train/val for TL
                    df = df.sample(frac=1, random_state=42)  # Shuffle
                    train_size = int(0.8 * len(df))
                    train_df = df.iloc[:train_size]
                    val_df = df.iloc[train_size:]
                    smi_save_path = csv_save_path.rsplit('.', 1)[0] + "_tl.smi"
                    val_smi_save_path = csv_save_path.rsplit('.', 1)[0] + "_val.smi"
                    train_df[['smiles']].to_csv(smi_save_path, sep='\t', index=False, header=False)
                    val_df[['smiles']].to_csv(val_smi_save_path, sep='\t', index=False, header=False)
                    toml_path = TL(smi_save_path, val_smi_save_path, TL_method=selected_TL_model_file) # selected_TL_model_fileを渡す
                    run_reinvent(toml_path)
                    st.success("Transfer learning from CSV completed!")
            else:
                st.error("CSV file must contain a 'smiles' column.")
        except Exception as e:
            st.error(f"Error reading CSV file: {e}")

# Sampling from TL model
st.header("Sampling from Transfer Learning Model")
tl_model_files = [f for f in os.listdir(model_dir) if f.endswith(".model")]
TL_model = st.selectbox("Select TL Model", options=tl_model_files, key="tl_model_select")


if st.button("Sample from TL Model", key="sample_from_tl_model_btn"):
    # st.session_stateから値を取得する前に、キーが存在するか確認
    if "tl_model_select" in st.session_state:
        selected_tl_model = st.session_state.tl_model_select
        toml_path, output_file = TL_toml_sampling( # TL_toml -> TL_toml_sampling に変更
            method=os.path.join(model_dir, selected_tl_model),
            # smiles_file=f"{input_dir}/parent.smi", # コメントアウトされているので、そのまま
            num_smiles=num_mols,
            device=device,
            show=False, # 引数として追加
            sample_strategy=sample_stategies,
            temperature=temperatures
        )
        run_reinvent(toml_path)
        st.success(f"Sampling from TL model {selected_tl_model} completed!")
    else:
        st.warning("Please select a TL Model first.")