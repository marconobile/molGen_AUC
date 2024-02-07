from pathlib import Path
from rdkit import DataStructs, Chem, RDLogger
import pandas as pd
from functools import partial
from source.argparser import argparser
from source.io_utils import get_dir_filename_ext
from source.rdkit_utils import get_valid_non_duplicates_mols, save_smiles
from source.multiprocessing_utils import apply_f_parallelized_batched
RDLogger.DisableLog('rdApp.*')

# Globals
nels_per_batch = 5000
workers = 4


def get_smi_natoms_fp(batch_mols, compute_fingerprint=Chem.RDKFingerprint): return [(Chem.MolToSmiles(m), m.GetNumAtoms(), compute_fingerprint(m)) for m in batch_mols]
def get_df(path, f=get_smi_natoms_fp):
    '''
    :param path: Path or str, filepath pointing at mols file
    :oaram f: function to be applied to a list of mols to build the out df here get_smi_natoms_fp given the requested cols of the df returned
    :return pandas df with cols "SMILES", "NATOMS", "FPs" for each smi in input file
    '''
    valid_non_duplicates_mols = get_valid_non_duplicates_mols(path)
    out = apply_f_parallelized_batched(
        f, valid_non_duplicates_mols, nels_per_batch=nels_per_batch, workers=workers)
    return pd.DataFrame(out, columns=["SMILES", "NATOMS", "FPs"])


def compute_exact_matches_batched(list_of_rows, ref_df=None): return [compute_exact_matches(row, ref_df) for row in list_of_rows]
def compute_exact_matches(other_row, ref_df):
    filter_ = ref_df[f"NATOMS"] == other_row["NATOMS"]
    possible_matches = ref_df[filter_]
    possible_matches["tanimoto"] = DataStructs.BulkTanimotoSimilarity(
        other_row["FPs"], possible_matches["FPs"].tolist())
    filter_ = possible_matches[f"tanimoto"] == 1
    exact_matches = possible_matches[filter_]
    return exact_matches



if __name__ == "__main__":

    # get paths from cmd line
    args = argparser("ref_path", "other_path")
    ref_path = Path(args.ref_path)
    other_path = Path(args.other_path)
    
    # preprocess ref data
    ref_df = get_df(ref_path)
    # step1: define function that is able to compute_exact_matches given a batch of other_df rows
    f = partial(compute_exact_matches_batched, ref_df=ref_df)
    
    other_df = get_df(other_path)

    # step2: preprocessing of other_df
    l = list(zip(*other_df.iterrows()))[1]

    # step3: apply compute_exact_matches in parallelized manner over a set of batches
    out = apply_f_parallelized_batched(f, l, nels_per_batch=nels_per_batch, workers=workers)

    # step4: postprocess of return to make it a list of smiles that are present on ref_df AND other_df
    exact_matches_smiles = [df["SMILES"].values[0] for df in out if not df.empty]

    _, filename_ref, _ = get_dir_filename_ext(ref_path)
    dir, filename_other, _ = get_dir_filename_ext(other_path)

    filename_out = f"{filename_ref}_intersection_{filename_other}"
    save_smiles(exact_matches_smiles, dir, filename_out)
