from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
from rdkit.Chem import DataStructs
import os
import sys
import time
import shutil
import multiprocessing as mp
from tqdm import tqdm
import gzip
from rdkit import RDLogger
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT, GetUSR
from urllib.request import urlopen
import zipfile

param_gen = AllChem.ETKDGv3()
param_gen.pruneRmsThresh=0.50
param_gen.numThreads=1
param_gen.maxIters=8192
param_gen.onlyHeavyAtomsForRMS=True

logfile = open("Log_file.log", 'at')
print("Start LOGGING:", file = logfile)
logfile.close()

chunk_size = 10*1024

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def gen_for_worker(file1, file2):
    n = 0
    with open(file1, 'rb', buffering=chunk_size) as f:
        while True:
            lines = f.readlines(chunk_size)
            if not lines:
                break
            n += 1
            yield lines, n, file2


def worker(line_n):
    lines, n, folder_out = line_n
    logging = list()
    failed = list()
    file_out = os.path.join(folder_out, "USRCAT_" + f'{n:08d}.bin.gz')
    print("FILE_OUT", file_out)
    size = 0
    with gzip.open(file_out, 'wb') as f:
        for line in lines:
            try:
                size += len(line)
                line = line.decode().rstrip('\n')
                smiles = line.split('\t', 1)[0]
                id = line.split('\t', 1)[1]
                smiles = sorted(smiles.split('.'), key=lambda x: len(x), reverse=True)[0]
                moll = Chem.MolFromSmiles(smiles)
                nr = int(AllChem.CalcNumRotatableBonds(moll))
                if nr <= 3:
                    nc = 50
                elif nr > 6:
                    nc = 300
                else:
                    nc = nr**3

                moll = Chem.AddHs(moll)
                cids_moll = AllChem.EmbedMultipleConfs(moll, numConfs = nc , params = param_gen)
                res_moll = AllChem.MMFFOptimizeMoleculeConfs(moll, maxIters = 8192, numThreads = 1, mmffVariant = 'MMFF94')
                logging.append('Mol ID: ' + str(id) + ' and SMILES: ' + str(smiles) + '\n')
                logging.append('Initial nc: ' + str(nc) + ' and ultimate pruned nc: ' + str(len(res_moll)) + '\n')
                moll = Chem.RemoveHs(moll)
                moll_confs = moll.GetConformers()
                pickle.dump((line, tuple(GetUSRCAT(moll, confId=i.GetId()) for i in moll_confs)), f)
            except:
                failed.append(line + '\n')
                print("Failed at: ", line)
    return size, logging, failed

def human_format(num: int) -> str:
    num = float(f'{num:.3g}')
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    num = f'{num:.1f}'.rstrip('0').rstrip('.')
    let = ['', 'K', 'M', 'B', 'T'][magnitude]
    return f'{num}{let}'

def add_count_to_file(file: str, count: int) -> str:
    name, ext = os.path.splitext(os.path.basename(file))
    if ext.lower() in ['.gz', '.bz2']:
        ext = '.' + name.split('.')[-1] + ext
        name = '.'.join(name.split('.')[:-1])
    folder = os.path.dirname(file)
    new_file = os.path.join(folder, f'{name}_{human_format(count)}{ext}')
    os.replace(file, new_file)
    return new_file

def download_sets(web_path, path_to_base):
    zipurl = ''
    zip_name = ''
    if web_path == 2:
        zipurl = ''
        zip_name = ''
    elif web_path == 3:
        zipurl = 'https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeAYk44LOtJJqQwsVPqKDx8BeaCQBizh155_tUfPj5l8aw?download=1'
    elif web_path == 4:
        zipurl = 'https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EZniFi17gP9OtUHokSINFA8BNgynKLXRWdocIt6eRK3tgA?download=1'
    elif web_path == 5:
        zipurl = 'https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EUyd0ux19QtMkzQta991MUABm7eHqkqREOc0ePFmaZZ9OQ?download=1'
    elif web_path == 6:
        zipurl = 'https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeSKrA8uYwZOgo35QNwLgQ8BCDlFHEY7f-TU7pb8OkTE5g?download=1'


    os.system('wget %s -O temp_usrcat.zip' %zipurl)
    with zipfile.ZipFile('temp_usrcat.zip', 'r') as zip_ref:
        zip_ref.extractall()
    os.remove('temp_usrcat.zip')

def do_base(files_in, base_folder, fail_file):
    global file_in
    global folder_out
    global file_failed

    file_in = files_in
    folder_out = base_folder
    file_failed = fail_file

    start = time.time()
    if not os.path.exists(folder_out):
        print(folder_out)
        os.makedirs(folder_out)
    err_counter = 0
    res_counter = 0
    with mp.Pool(processes=1) as pool, open(file_failed, 'at') as f, open("Log_file.log", 'at') as l, \
        tqdm(unit_scale=True, mininterval=1.0, total=os.path.getsize(file_in), unit='B', unit_divisor=1024) as t:
        for res in pool.imap(worker, gen_for_worker(file_in, folder_out)):
            t.update(res[0])
            l.writelines(res[1])
            res_counter += res[0]
            err_counter += len(res[2])
            f.writelines(res[2])
    add_count_to_file(file_failed, err_counter)
    end = time.time()
    total_time = (end - start) / 3600

    logfile = open("Log_file.log", 'at')
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours', file = logfile)
    logfile.close()
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours')

if __name__ == '__main__':

    start = time.time()
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    err_counter = 0
    res_counter = 0
    with mp.Pool(processes=58) as pool, open(file_failed, 'at') as f, open("Log_file.log", 'at') as l, \
        tqdm(unit_scale=True, mininterval=1.0, total=os.path.getsize(file_in), unit='B', unit_divisor=1024) as t:
        for res in pool.imap(worker, gen_for_worker(file_in)):
            t.update(res[0])
            l.writelines(res[1])
            res_counter += res[0]
            err_counter += len(res[2])
            f.writelines(res[2])
    add_count_to_file(file_failed, err_counter)
    end = time.time()
    total_time = (end - start) / 3600

    logfile = open("Log_file.log", 'at')
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours', file = logfile)
    logfile.close()
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours')
