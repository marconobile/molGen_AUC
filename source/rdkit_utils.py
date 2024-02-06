import os
from rdkit import Chem
from typing import List, Type, Callable, Union, Iterable
from source.io_utils import get_dir_filename_ext, delete_file
from pathlib import Path
from source.logger import Logger


def mols_from_file(pathfile: Union[str, Path], drop_none: bool = False) -> list[Chem.rdchem.Mol]:
    '''
    :param pathfile: path/to/file.ext as str or as pathlib.Path obj
    ext can be: [.sdf, .csv, .txt, .smiles]; if ext not in supported extension: raise TypeError
    :param drop_none: drops mols non valid for rdkit
    :return list of mols from file.ext
    '''
    if isinstance(pathfile, Path):
        pathfile = str(pathfile)
    filename_ext = os.path.splitext(pathfile)[-1].lower()
    if filename_ext in ['.sdf']:
        suppl = Chem.SDMolSupplier(pathfile)
    elif filename_ext in ['.csv', '.txt', '.smiles']:
        suppl = Chem.SmilesMolSupplier(pathfile, titleLine=False)
    else:
        raise TypeError(f"{filename_ext} not supported")
    if drop_none:
        return [x for x in suppl if x is not None]
    return [x for x in suppl]


def keep_valid_mols(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    return [m for m in mols if (m and validate_rdkit_mol(m))]


def validate_rdkit_mol(mol: Chem.Mol) -> bool:
    """
    Sanitizes an RDKit molecules and returns True if the molecule is chemically
    valid.
    :param mol: an RDKit molecule 
    :return: True if the molecule is chemically valid, False otherwise
    """
    if Chem is None:
        raise ImportError('`validate_rdkit_mol` requires RDkit.')
    if len(Chem.GetMolFrags(mol)) > 1:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except ValueError:
        return False


def read_smiles_from_file(pathfile: str):
    _, _, ext = get_dir_filename_ext(pathfile)
    if ext not in ["txt", "smiles"]:
        raise TypeError(f"extension {ext} not valid")
    with open(pathfile) as file:
        return [line.rstrip() for line in file]


def smi2mols(smiles: Iterable[str]) -> List[Chem.Mol]:
    return [Chem.MolFromSmiles(smi) for smi in smiles]


def mols2smi(mols: List[Chem.Mol]) -> List[str]:
    return [Chem.MolToSmiles(mol) for mol in mols]


def save_smiles(smiles: Iterable[str], path: str = ".", filename: str = "smiles", ext: str = '.smiles') -> str:
    '''
    :param smiles: iterable of smiles
    saves smiles in path/filename.ext
    extension can be provided in filename or as separate arg
    args:
        - smiles str iterable 
        - path directory where to save smiles 
        - filename name of the file, must not have extension
    '''
    path_to_file = os.path.join(path, filename)
    filename_ext = os.path.splitext(path_to_file)[-1].lower()
    if not filename_ext:
        if ext not in ['.txt', '.smiles']:
            raise TypeError(f"extension {ext} not valid")
        path_to_file += ext

    with open(path_to_file, "w+") as f:
        f.writelines("%s\n" % smi for smi in smiles)
    return path_to_file


def drop_duplicates_with_openbabel(in_file: str):
    '''
    Given an input .smiles file, generate an out_file.smiles with no duplicates or invalid mols
    '''
    dir, filename, ext = get_dir_filename_ext(in_file)
    out = Path(dir)/f"{filename}_no_duplicates.{ext}"
    cmd = f"obabel -ismiles {in_file} -osmiles -O {str(out)} --unique"
    os.system(cmd)  # synchronous call, the result is waited
    return out


def get_valid_non_duplicates(ref_path, return_dtype='smiles'):
    '''
    Given an input file of mols or smiles:
    1) filters out valid and sanitizable mols keep_valid_mols
    2) saves the smiles of valid mols to file via save_smiles
    3) uses obabel to drop duplicate mols
    4) deletes in/out temp files used by obabel
    5) returns list containing valid_non_duplicates smiles if smiles=True, else List[Chem.Mol]
    '''
    if not isinstance(ref_path, str) and not isinstance(ref_path, Path):
        raise TypeError(
            f"get_valid_non_duplicates() works only on files containing SMILES")
    if not isinstance(return_dtype, str) or return_dtype not in ["smiles", "mols"]:
        raise TypeError(
            f"return_dtype passed is not of type string. It must be either \"smiles\" or \"mols\"")
    dir, name, ext = get_dir_filename_ext(ref_path)
    ref_mols = mols_from_file(ref_path)
    log = Logger(dir, name=name+"_log.txt")
    log.write(f"Mols in {ref_path}: {len(ref_mols)}")
    valid_mols = keep_valid_mols(ref_mols)
    log.write(f"Valid mols in {ref_path}: {len(valid_mols)}")
    valid_smiles = mols2smi(valid_mols)
    path_to_saved_valid_smiles = save_smiles(
        valid_smiles, dir, f'{name}_valid')
    valid_non_duplicates_path = drop_duplicates_with_openbabel(
        path_to_saved_valid_smiles)
    valid_non_duplicates = read_smiles_from_file(valid_non_duplicates_path)
    log.write(
        f"Non duplicate mols in {valid_non_duplicates_path}: {len(valid_non_duplicates)}")
    delete_file(path_to_saved_valid_smiles)
    if return_dtype.lower()=='smiles': 
        out = read_smiles_from_file(valid_non_duplicates_path)
    else: 
        out = mols_from_file(valid_non_duplicates_path)
    delete_file(valid_non_duplicates_path)
    return out