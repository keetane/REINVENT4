import os
import subprocess
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, Recap, Descriptors, rdMolDescriptors
import pubchempy as pcp
from datetime import datetime

# ─── Page & Session State Setup ──────────────────────────────────────────────
st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded",
)

if 'wd' not in st.session_state:
    st.session_state.wd = os.getcwd()
wd = st.session_state.wd

now = datetime.now()
time = now.strftime("%Y%m%d_%H%M")

# ─── Directory Setup ─────────────────────────────────────────────────────────
input_dir    = os.path.join(wd, 'input');    os.makedirs(input_dir, exist_ok=True)
results_dir  = os.path.join(wd, 'results');  os.makedirs(results_dir, exist_ok=True)
toml_dir     = os.path.join(wd, 'toml');     os.makedirs(toml_dir, exist_ok=True)
sampling_log = os.path.join(results_dir, 'log'); os.makedirs(sampling_log, exist_ok=True)

priors_dir = os.path.join(os.getcwd(), "priors")
priors = {
    # "LibInvent":    os.path.join(priors_dir, "libinvent.prior"),
    # "LinkInvent":   os.path.join(priors_dir, "linkinvent.prior"),
    "mol2mol_high": os.path.join(priors_dir, "mol2mol_high_similarity.prior"),
    "mol2mol_med":  os.path.join(priors_dir, "mol2mol_medium_similarity.prior"),
    "mol2mol_mmp":  os.path.join(priors_dir, "mol2mol_mmp.prior"),
    "scaffold_gen": os.path.join(priors_dir, "mol2mol_scaffold_generic.prior"),
    "scaffold":     os.path.join(priors_dir, "mol2mol_scaffold.prior"),
    "similarity":   os.path.join(priors_dir, "mol2mol_similarity.prior"),
    "PubChem":      os.path.join(priors_dir, "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior"),
    "Reinvent":     os.path.join(priors_dir, "reinvent.prior"),
}

# ─── Parent Molecule Drawer (2:1 Split) ────────────────────────────────────────
st.header("Parent Molecule Drawer")

# Create two columns, left (2/3) and right (1/3)
col1, col2 = st.columns([2, 1])

with col1:
    # --- PubChem から SMILES を取得 ---
    if 'compound_name' not in st.session_state:
        st.session_state.compound_name = "ruxolitinib"
    if 'smiles' not in st.session_state:
        st.session_state.smiles = "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"

    # 分子名入力
    compound_name = st.text_input(
        "Enter compound name",
        value=st.session_state.compound_name
    )
    # 取得ボタン
    if st.button("Fetch SMILES from PubChem"):
        try:
            new_smiles = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
            st.session_state.smiles = new_smiles
            st.session_state.compound_name = compound_name
            st.success(f"Fetched SMILES for {compound_name}")
        except Exception as e:
            st.error(f"Error fetching SMILES: {e}")

    # SMILES 入力／編集
    smiles = st.text_area(
        "SMILES",
        value=st.session_state.smiles,
        height=80
    )
    st.session_state.smiles = smiles

    # テキストファイルとして保存
    with open(os.path.join(input_dir, "parent.smi"), "w") as f:
        f.write(smiles + "\n")

with col2:
    # SMILES から構造を描画
    mol_parent = Chem.MolFromSmiles(st.session_state.smiles)
    if mol_parent:
        img = Draw.MolToImage(mol_parent, size=(300, 300))
        st.image(img, caption="Parent Molecule")
    else:
        st.warning("Invalid SMILES: 構造を生成できません")

# # ─── Parent Molecule Drawer ──────────────────────────────────────────────────
# st.header("Parent Molecule Drawer")

# if 'compound_name' not in st.session_state:
#     st.session_state.compound_name = "ruxolitinib"
# if 'smiles' not in st.session_state:
#     st.session_state.smiles = "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"

# compound_name = st.text_input("Enter compound name", value=st.session_state.compound_name)
# if st.button("Fetch SMILES from PubChem"):
#     try:
#         new = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
#         st.session_state.smiles = new
#         st.session_state.compound_name = compound_name
#         st.success(f"Fetched SMILES for {compound_name}")
#     except Exception as e:
#         st.error(f"Error fetching SMILES: {e}")

# smiles = st.text_area("SMILES", value=st.session_state.smiles, height=80)
# st.session_state.smiles = smiles

# mol_parent = Chem.MolFromSmiles(smiles)
# if mol_parent:
#     st.image(Draw.MolToImage(mol_parent, size=(300,300)), caption="Parent Molecule")

# with open(os.path.join(input_dir, "parent.smi"), "w") as f:
#     f.write(smiles + "\n")

# ─── Recap Decomposition ─────────────────────────────────────────────────────
st.header("Children by Recap Decomposition")
recap_tree = Recap.RecapDecompose(mol_parent)

warhead_mols, warhead_smis = [], []
child_mols,   child_smis   = [], []

def relabel_dummy_atoms(mol, map_num=1):
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(map_num)
    return rw.GetMol()

for node in recap_tree.children.values():
    m = node.mol
    warhead_mols.append(m)
    warhead_smis.append(Chem.MolToSmiles(m))
    child = relabel_dummy_atoms(m, map_num=1)
    child_mols.append(child)
    child_smis.append(Chem.MolToSmiles(child))

BB_namelist = [f"BB{i+1}" for i in range(len(warhead_mols))]
if warhead_mols:
    st.image(
        Draw.MolsToGridImage(
            warhead_mols, molsPerRow=3, subImgSize=(600,300), legends=BB_namelist
        ),
        caption="Recap Fragments"
    )

# ─── Sidebar: Sampling Parameters ─────────────────────────────────────────────
st.sidebar.header("Sampling Parameters")
num_mols         = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device           = st.sidebar.selectbox("Device", ["cpu", "cuda", "mps"])
unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)
overwrite        = st.sidebar.checkbox("Overwrite outputs", value=True)
comment = st.sidebar.text_input("File Comment")

st.sidebar.header("mol2mol Options")
sample_strategy  = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", "multinomial"])
temperature      = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)

# ─── REINVENT4 Runner ────────────────────────────────────────────────────────
def run_reinvent(method, smiles_file, num_smiles, device, sample_strategy=None, temperature=1.0):
    fname = os.path.basename(method).replace(".prior","")
    if not overwrite:
        fname += f"_{time}"
    if comment!= "":
        fname = str(comment) + "_" + fname
    output_csv = os.path.join(results_dir, f"{fname}.csv")

    toml = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "{sampling_log}/_sampling.json"

[parameters]
model_file = "{method}"
output_file = "{output_csv}"
num_smiles = {num_smiles}
"""
    if unique_molecules: toml += "unique_molecules = true\n"
    if randomize_smiles: toml += "randomize_smiles = true\n"
    if smiles_file:      toml += f'smiles_file = "{smiles_file}"\n'
    if sample_strategy == "beamsearch":
        toml += 'sample_strategy = "beamsearch"\n'
    elif sample_strategy == "multinomial":
        toml += 'sample_strategy = "multinomial"\n'
        toml += f"temperature = {temperature}\n"

    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml)

    log_file = os.path.join(sampling_log, "sampling.log")
    subprocess.call(["reinvent", "-l", log_file, toml_path])

    df = pd.read_csv(output_csv)
    df['ROMol'] = df['SMILES'].map(lambda s: Chem.MolFromSmiles(s))

    # ── Descriptor Calculation ─────────────────────────────
    df['MW']         = df['ROMol'].map(lambda m: Descriptors.MolWt(m)            if m else None)
    df['LogP']       = df['ROMol'].map(lambda m: Descriptors.MolLogP(m)         if m else None)
    df['HBD']        = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumHBD(m) if m else None)
    df['HBA']        = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumHBA(m) if m else None)
    df['TPSA']       = df['ROMol'].map(lambda m: rdMolDescriptors.CalcTPSA(m)   if m else None)
    df['CtRotBonds']     = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumRotatableBonds(m) if m else None)
    df['CtRings']    = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumRings(m)           if m else None)
    df['CtAmides']   = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumAmideBonds(m)     if m else None)
    df['CtAromaticRings']      = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumAromaticRings(m)  if m else None)
    df['CtHetAromaticRings']   = df['ROMol'].map(lambda m: rdMolDescriptors.CalcNumHeteroatoms(m)    if m else None)

    df.drop(columns=['ROMol'], inplace=True)
    df.to_csv(output_csv, index=False)

    return output_csv

# ─── Main Sampling UI ────────────────────────────────────────────────────────
st.header("ReInvent Sampling")
col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("Mol2Mol Sampling")
    model_key        = st.selectbox("Model", list(priors.keys()), index=0)
    selected_prior   = priors[model_key]
    if st.button("Run Mol2Mol", key="mol2mol"):
        out = run_reinvent(
            # method       = priors["mol2mol_high"],
            method       = selected_prior,
            smiles_file  = os.path.join(input_dir, "parent.smi"),
            num_smiles   = num_mols,
            device       = device,
            sample_strategy = sample_strategy,
            temperature     = temperature
        )
        st.success(f"Mol2Mol sampling done: {out}")

with col2:
    st.subheader("LibInvent Sampling")
    choice = st.multiselect("Select Child", BB_namelist, default=BB_namelist[:1], key="child_sel")
    if choice:
        idx = BB_namelist.index(choice[0])
        with open(os.path.join(input_dir, "child.smi"), "w") as f:
            f.write(child_smis[idx] + "\n")
    if st.button("Run LibInvent", key="libinv"):
        out = run_reinvent(
            method      = os.path.join(priors_dir, "libinvent.prior"),
            smiles_file = os.path.join(input_dir, "child.smi"),
            num_smiles  = num_mols,
            device      = device
        )
        st.success(f"LibInvent sampling done: {out}")

with col3:
    st.subheader("LinkInvent Sampling")
    choice = st.multiselect("Select Warheads", BB_namelist, default=BB_namelist[:2], key="warhead_sel")
    if choice:
        war_smis = [warhead_smis[BB_namelist.index(w)] for w in choice]
        with open(os.path.join(input_dir, "warheads.smi"), "w") as f:
            f.write("|".join(war_smis))
    if st.button("Run LinkInvent", key="linkinv"):
        out = run_reinvent(
            method      = os.path.join(priors_dir, "linkinvent.prior"),
            smiles_file = os.path.join(input_dir, "warheads.smi"),
            num_smiles  = num_mols,
            device      = device
        )
        st.success(f"LinkInvent sampling done: {out}")
