from pathlib import Path
from source.rdkit_utils import get_valid_non_duplicates, save_smiles
from source.io_utils import delete_file, get_dir_filename_ext, get_all_ext_files_in_path
from source.argparser import argparser


# get paths from cmd line
args = argparser("ref_path", "train_path", "gen_folder")
ref_path = Path(args.ref_path)
train_path = Path(args.train_path)
gen_folder = Path(args.gen_folder)

# get path of all files in gen folder
files = get_all_ext_files_in_path(Path(args.gen_folder), ".smiles")
if len(files) != 5:
    raise ValueError(f"Number of files in {args.gen_folder} must be 5 and they must all have .smiles extension, filenames must do not have blank spaces")

print(files)

filelist = {"ref": ref_path}
for i, file in enumerate(files):
    filelist[f"file_{i}"] = Path(file)

# step 1: get val-unique for gen data and ref
for k, v in filelist.items():
    filelist[k] = get_valid_non_duplicates(v)

# step 2: get intersection ref/train
    








