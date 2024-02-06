import os
from pathlib import Path
from typing import List, Type, Callable, Union


def get_dir_filename_ext(pathfile: Union[str, Path], ext: bool = False) -> (str, str, str):
    '''
    :param pathfile a path/to/file.ext
    :return (path/to/, file, ext)    
    '''
    if isinstance(pathfile, Path):
        pathfile = str(pathfile)
    dir, filename = os.path.split(pathfile)
    filename, ext = filename.split(".")
    return dir, filename, ext.lower()


def delete_file(path_to_file):
    if os.path.isfile(path_to_file):
        try:
            os.remove(path_to_file)
        except OSError:
            raise f"{path_to_file} already existing and could not be removed"

        
def get_all_ext_files_in_path(path, ext):
    if isinstance(path, Path): 
        path = str(path)
    file_list = []
    for root, _, files in os.walk(path):
        for file in files:
            if file.endswith(ext):
                file_list.append(os.path.join(root, file))
    return file_list