import glob, os, shutil, itertools, copy
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from tqdm import tqdm
import pandas as pd
from matplotlib import pyplot as plt
from . import logfile_process, FormatConverter, xtb_process, mol_manipulation, Tool
import seaborn as sns

def generate_combinations(reactant_file, result_file):
    """Analysis Reactants, Get React Atom ID

    Args:
        reactant_file (_type_): _description_
        result_file (_type_): _description_
    """  
    react_csv = pd.read_csv(reactant_file)
    BN_smiles = react_csv["B_N"].dropna().to_numpy()
    BN_index = react_csv['B_N_id'].dropna().to_numpy()
    S_smiles = react_csv["Thiol"].dropna().to_numpy()
    S_index = react_csv['Thiol_id'].dropna().to_numpy()
    sol_smiles = react_csv["solvent"].dropna().to_numpy()
    sol_index = react_csv['solvent_id'].dropna().to_numpy()
    ini_smiles = react_csv["ini"].dropna().to_numpy()
    ini_index = react_csv['ini_id'].dropna().to_numpy()
    ini_temp = react_csv['ini_t'].dropna().to_numpy()

    if os.path.exists(result_file):
        result_csv = pd.read_csv(result_file).to_dict()
    else:
        result_csv = {"Smiles":{}, "Index":{}, "Atomid":{}, "Type":{}}
    result_csv_id = len(result_csv["Index"])
    types_ = ["B", "S", "sol", "ini"]
    for all_smiles, all_ids, type_ in zip([BN_smiles, S_smiles, sol_smiles, ini_smiles], [BN_index, S_index, sol_index, ini_index], types_):
        for id_, [smiles, index] in enumerate(zip(all_smiles, all_ids)):
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            if smiles in result_csv["Smiles"]:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if type_ in ['B', "S"]:
                atom_idx = [each for each in mol.GetAtoms() if each.GetSymbol() == type_][0].GetIdx()
            elif type_ in ['ini']:
                atom_idx = ini_temp[id_]
            else:
                atom_idx = -1
            result_csv["Smiles"][result_csv_id] = smiles
            result_csv['Index'][result_csv_id] = int(index)
            result_csv['Atomid'][result_csv_id] = atom_idx
            result_csv['Type'][result_csv_id] = type_
            result_csv_id += 1
            if type_ == 'B':
                rwmol = Chem.RWMol(mol)
                new_atom_id = rwmol.AddAtom(Chem.Atom(17))
                rwmol.AddBond(new_atom_id, atom_idx, Chem.BondType.SINGLE)
                new_mol = rwmol.GetMol()
                new_smiles = Chem.MolToSmiles(new_mol)
                new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(new_smiles))
                result_csv["Smiles"][result_csv_id] = new_smiles
                result_csv['Index'][result_csv_id] = int(index) + 1000
                result_csv['Atomid'][result_csv_id] = new_atom_id
                result_csv['Type'][result_csv_id] = 'BNCl'
                result_csv_id += 1
    result_csv = pd.DataFrame(result_csv)
    result_csv.to_csv(result_file, index=False)

def B_N_Single_Xtb(root_file, result_file, mol_xtb_name = 'Mol_xtb', mol_name = "Mols"):
    """产生B自由基、配体单体的Xtb优化文件

    Args:
        result_file (str): 存储信息的文件路径
    """    
    result_file=pd.read_csv(result_file)
    mol_xtb_file = os.path.join(root_file, mol_xtb_name)
    old_mol_file = os.path.join(root_file, 'Mols')
    mol_file = os.path.join(root_file, mol_name)
    if not os.path.isdir(mol_file):
        os.mkdir(mol_file)
    if not os.path.isdir(mol_xtb_file):
        os.mkdir(mol_xtb_file)
    all_B_mols = []
    all_B_names = []
    all_Other_mols = []
    all_Other_names = []
    for ix, row in result_file.iterrows():
        mol_idx = row['Index']
        mol_name = f"{mol_idx:05}"
        type_ = row['Type']
        if os.path.exists(os.path.join(old_mol_file, f"{mol_name}.mol")):
            continue
        mol = mol_manipulation.smiles2mol(row['Smiles'], conf_num=1)
        mol_atom_idx = row["Atomid"]
        AllChem.UFFOptimizeMolecule(mol)
        if type_ == "B":
            all_B_mols.append(mol)
            all_B_names.append(mol_name)
        else:
            all_Other_mols.append(mol)
            all_Other_names.append(mol_name)
        Chem.MolToMolFile(mol, os.path.join(mol_file, f"{mol_name}.mol"))
    mol_xtb_file_ = os.path.join(mol_xtb_file, 'B')
    xtb_process.xtb_main(all_B_names, all_B_mols, dir_path=mol_xtb_file_, core=60, uhf=1)
    # xtb_process.shift_to_sugan(mol_xtb_file_, 1, 0, 1)
    mol_xtb_file_ = os.path.join(mol_xtb_file, 'Other')
    xtb_process.xtb_main(all_Other_names, all_Other_mols, dir_path=mol_xtb_file_, core=60, uhf=0)
    # xtb_process.shift_to_sugan(mol_xtb_file_, 1, 0, 0)

def smiles_DFT_calc(root_dir='first_xtb', 
                    mol_dir='mol', 
                    dft_dir='mol_dft', 
                    method="opt freq b3lyp/6-31g* em=gd3bj",
                    conf_limit=3,
                    rmsd_limit=1.5,
                    SpinMultiplicity = None
                    ):
    """通过Xtb结果，优化得到Gaussian优化输入文件

    Args:
        root_dir (str, optional): Xtb的根目录. Defaults to 'first_xtb'.
        mol_dir (str, optional): Mol分子的目录. Defaults to 'mol'.
        dft_dir (str, optional): 要存储Gaussian输入文件的目录. Defaults to 'mol_dft'.
        method (str, optional): Gaussian方法. Defaults to "opt freq b3lyp/6-31g* em=gd3bj".
        conf_limit (int, optional): Xtb读取结构的构象数量限制. Defaults to 3.
        rmsd_limit (float, optional): Xtb读取结构的RMSD限制. Defaults to 1.5.
        SpinMultiplicity (int, optional): 设定的自旋多重度. Defaults to None.
    """                    
    all_files = glob.glob(root_dir + "/*/*/*")
    for xtb_file in all_files:
        if ("crest.out" in xtb_file) or ("best" in xtb_file) or ("crest_conf" in xtb_file):
            pass
        else:
            if os.path.isdir(xtb_file):
                shutil.rmtree(xtb_file)
            else:
                os.remove((xtb_file))
    xtb_dirs = glob.glob(root_dir + "/*/*")
    for i, xtb_dir in enumerate(xtb_dirs):
        mol_name = os.path.split(xtb_dir)[-1][:-2]
        mol_file = mol_dir + f"/{mol_name}.mol" 
        mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
        title = "Singlemol"
        xtb_process.after_xtb(mol,xtb_dir=xtb_dir, save_dir=dft_dir, xtb_title=title, method=method, conf_limit=conf_limit, rmsd_limit=rmsd_limit, SpinMultiplicity=SpinMultiplicity)


def SPE_DFT_calc(target_dir, opt_name="Reactants", eng_name='Reactants_eng', mol_name="Mols", save_chk=None, method="b3lyp/6-311+g(d,p) em=gd3bj"):
    opt_file_dir = os.path.join(target_dir, opt_name)
    eng_dir = os.path.join(target_dir, eng_name)
    mol_files = glob.glob(os.path.join(target_dir, mol_name, "*.mol"))
    for mol_file in mol_files:
        log_files = glob.glob(opt_file_dir + "/" + os.path.split(mol_file)[-1].split(".")[0] + "*.log")
        if len(log_files) == 0:
            continue 
        for log_file in log_files:
            new_log_name = eng_dir + "/" + os.path.split(log_file)[-1].split('.')[0] + ".gjf" 
            opt_log = logfile_process.Logfile(log_file, mol_file_dir=mol_file)
            assert len(opt_log.running_positions) != 0
            title, charge, symbol_list, position,= opt_log.title, opt_log.charge, opt_log.symbol_list, opt_log.running_positions[-1]
            title = " ".join(str(each) for each in title)
            if save_chk:
                savechk = os.path.split(new_log_name.strip(".gjf"))[-1]
            else:
                savechk = None
            FormatConverter.block_to_gjf(symbol_list, position, new_log_name, charge, opt_log.multiplicity, title,
                        method=method, savechk=savechk)


def SPE_DFT_calc_wfn(target_dir, opt_name="Reactants", eng_name='Reactants_eng', mol_name="Mols", save_chk=None, method="b3lyp/6-311+g(d,p) em=gd3bj"):
    opt_file_dir = os.path.join(target_dir, opt_name)
    eng_dir = os.path.join(target_dir, eng_name)
    mol_files = glob.glob(os.path.join(target_dir, mol_name, "*.mol"))
    for mol_file in tqdm(mol_files):
        log_files = glob.glob(opt_file_dir + "/" + os.path.split(mol_file)[-1].split(".")[0] + "*.log")
        if len(log_files) == 0:
            continue 
        for log_file in log_files:
            new_log_name = eng_dir + "/" + os.path.split(log_file)[-1].split('.')[0] + ".gjf" 
            opt_log = logfile_process.Logfile(log_file, mol_file_dir=mol_file)
            assert len(opt_log.running_positions) != 0
            wfn_name = os.path.split(log_file)[-1].split('.')[0] + ".wfn"
            opt_log = logfile_process.Logfile(log_file, mol_file_dir=mol_file)
            assert len(opt_log.running_positions) != 0
            title, charge, symbol_list, position,= opt_log.title, opt_log.charge, opt_log.symbol_list, opt_log.running_positions[-1]
            title = " ".join(str(each) for each in title)
            FormatConverter.block_to_gjf(symbol_list, position, new_log_name, charge, opt_log.multiplicity, title,
                        method=method, final_line=wfn_name)


def collection_dft_single(result_path, mol_dir, dft_dir, spe_dir):
    # Bresult_path ='Data/All_Data/reactants_B.csv'
    # Nresult_path = 'Data/All_Data/reactants_N.csv'
    # Clresult_path = 'Data/All_Data/reactants_Cl.csv'
    error_reason, E_energy, G_energy, conf_idxs = [],[],[], []
    result_file = pd.read_csv(result_path)
    for line_id, line in tqdm(result_file.iterrows()):
        dft_dir_ = dft_dir
        spe_dir_ = spe_dir
        smiles = line['Smiles']
        Index = line['Index']
        atom_id = line['Atomid']
        mol_file = glob.glob(os.path.join(mol_dir, f'{Index:05}.mol'))
        if len(mol_file) == 0:
            print(smiles, Index, "Is Error")
            # return None
        mol_file = mol_file[0]
        opt_files = glob.glob(os.path.join(dft_dir_, f'{Index:05}*.log'))
        temp_idx, temp_E, temp_G = [], [], []
        for opt_file in opt_files:
            opt = mol_manipulation.logfile_process.Logfile(opt_file)
            spe_files = glob.glob(os.path.join(spe_dir_, os.path.split(opt_file)[-1]))
            if len(spe_files) == 0:
                print(smiles, Index, "OPT Error")
                continue
            conf_id = int(opt_file.split('.')[0].split("_")[-1])
            spe_file = spe_files[0]
            spe = mol_manipulation.logfile_process.Logfile(spe_file)
            electric_energy = spe.all_engs[0]
            G_cor = opt.all_engs[-1]
            temp_idx.append(conf_id)
            temp_E.append(electric_energy)
            temp_G.append(G_cor + electric_energy)
        if len(temp_G) == 0:
            error_reason.append("DFT Error")
            E_energy.append(np.nan)
            G_energy.append(np.nan)
            conf_idxs.append(np.nan)
        else:
            min_index = np.argmin(temp_G)
            E_energy.append(temp_E[min_index])
            G_energy.append(temp_G[min_index])
            conf_idxs.append(temp_idx[min_index])
            error_reason.append(np.nan)
    result_file[f"E_energy"] = E_energy
    result_file[f"G_energy"] = G_energy
    result_file[f"conf_idxs"] = conf_idxs
    result_file[f"error_reason"] = error_reason
    result_file.to_csv(result_path, index=False)

def descriptor_generator(csv_file, index_map, Cl_des_map):
    if type(csv_file) == str:
        result_csv = pd.read_excel(csv_file).dropna()
    else:
        result_csv = csv_file.dropna()
    all_X = []
    yields = []
    for row_id, row in result_csv.iterrows():
        desc = []
        B = int(row['B_Index'])
        B_eq = row['eqs']
        ini = int(row['ini_Index'])
        sol_1 = int(row['sol_Index'])
        S = int(row['S_Index'])
        if "Cl_Index" in result_csv.columns:
            Cl = int(row['Cl_Index'])
            Cl_atomid = int(row['Cl_Atomid'])
            if Cl in [5001, 5002, 5003, 5004]:
                desc += Cl_des_map[f"Cl_{Cl:05d}_Claid_{Cl_atomid:05d}"]
            else:
                desc += Cl_des_map[f"Cl_{Cl:05d}"]
        desc += index_map[B]
        # desc += index_map[B + 1000]
        desc += index_map[ini]
        desc += index_map[sol_1]
        desc += index_map[S]
        desc += [B_eq]
        # desc += [1 if each == B_eq else 0 for each in [1.5, 2, 2.5,3]]
        all_X.append(desc)
        yields.append(row['yield'])
    yields = np.array(yields)
    all_X_qm = np.array(all_X)
    return all_X_qm, yields

def draw_correlation_map(X, y=None, figure_size=(5, 5), colors='coolwarm', 
                         useSVG=False, save_name='test', annot=True, show_label=False):
    """
    绘制相关性热图
    X: 特征矩阵 (pandas DataFrame 或 numpy array)
    y: 目标变量 (可选, pandas Series 或 numpy array)。如果传入y，则对角线显示 |corr(X_i, y)|
    """
    df = pd.DataFrame(X)
    
    # 计算特征之间的绝对相关性矩阵
    correlation_matrix = np.abs(df.corr())
    
    # 如果传入了 y，则把对角线替换为每个特征与 y 的绝对相关性
    if y is not None:
        y = pd.Series(y).reset_index(drop=True)   # 确保对齐
        corr_with_y = np.abs(df.corrwith(y))
        # 把相关性矩阵的对角线替换为与 y 的相关性
        np.fill_diagonal(correlation_matrix.values, corr_with_y.values)
    
    # 生成上三角掩码（不显示对称部分）
    mask = np.tril(np.ones_like(correlation_matrix, dtype=bool))
    
    print("最大相关系数（不含对角线）:", 
          np.max(np.nan_to_num(correlation_matrix.to_numpy()[mask], 0)))
    
    # 绘图
    f, ax = plt.subplots(figsize=figure_size, dpi=300)
    annot_kws = {"fontsize": 10}
    
    ax = sns.heatmap(correlation_matrix, 
                     mask=~mask,
                     cmap=colors,          # 使用传入的 colors 参数
                     annot=annot,         
                     fmt='.2f',            # 建议改成 .2f，更清晰
                     center=0,           
                     cbar=True,
                     annot_kws=annot_kws,
                     linewidths=0.5,       # 增加轻微网格线，更美观
                     linecolor='white')
    
    # 颜色条字体大小
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10)
    
    if not show_label:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    else:
        # 可选：旋转标签避免重叠
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    
    plt.tight_layout()
    
    if useSVG:
        plt.savefig(f"{save_name}.svg", format="svg", bbox_inches='tight', dpi=300)
    else:
        plt.savefig(f"{save_name}.png", dpi=300, bbox_inches='tight')
    
    # plt.close()   # 防止在循环中显示过多图像

def calc_distribution2(y, eachsize=0.01, title=None, xlab=None, ylab="Count", y_max=None, y_min=None, figure_size = (4,3)):
    if y_max == None:    y_max = np.max(y)
    if y_min == None:    y_min = np.min(y)
    X = np.arange(y_min, y_max + eachsize, eachsize)
    des = [0 for each in X]
    z = (y - y_min)/eachsize
    for each in z:
        try:
            des[int(each)] += 1
        except:
            continue
    des = np.array(des)
    # des = des / len(y)
    
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111)
    ax.patch.set_alpha(0.0)
    plt.bar(X, des, width=eachsize/2, color="green")
    plt.xlim(y_min - eachsize, y_max + eachsize)
    plt.ylim(0, np.max(des) * 1.2)
    plt.xlabel(xlab, fontsize=30)
    plt.ylabel(ylab, fontsize=30)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if title != None:
        plt.title = title
    plt.tight_layout()
    plt.savefig('test.svg', format='svg')
    plt.show()
    return des 


# Model 
def normalize_axis(arr, axis=0, mean=[], std=[]):
    """
    对数组中的某一维进行标准化（z-score normalization）
    
    参数：
    arr: ndarray，输入的数组
    axis: int，标准化的维度
    
    返回值：
    normalized_arr: ndarray，标准化后的数组
    """
    if len(mean) == 0 or len(mean) == 0:
        mean = np.mean(arr, axis=axis, keepdims=True)  # 计算均值
        std = np.std(arr, axis=axis, keepdims=True)  # 计算标准差
    normalized_arr = (arr - mean) / std  # 标准化
    normalized_arr = np.nan_to_num(normalized_arr, 0)
    return normalized_arr, mean, std