import os
import sys
from glob import glob
from vaspio.incar import Incar
from vaspio.cluster_data import *


def submit_job(job_path, job_name, job_type='native'):
    # todo: add functions to change memory, wall-time limit, number of nodes ... etc
    init_dir = os.getcwd()
    os.chdir(job_path)
    if job_type == 'native':
        os.system(f"cp {path_submission_vasp_native}/{name_submission_vasp} {job_path}")
    elif job_type == 'vtst':
        os.system(f"cp {path_submission_vasp_vtst}/{name_submission_vasp} {job_path}")
    elif job_type == 'ase':
        os.system(f"cp {path_submission_vasp_ase}/{name_submission_vasp} {job_path}")
    else:
        sys.exit(f"Invalid job type: {job_type}. Choose 'native', 'vtst' or 'ase' ")
    if job_scheduler == 'sge':
        os.system(f'qsub -N {job_name} {name_submission_vasp}')
    elif job_scheduler == 'slurm':
        os.system(f'sbatch --job-name {job_name} {name_submission_vasp}')
    else:
        sys.exit('Invalid job_scheduler')
    os.chdir(init_dir)


def continue_job(job_path, job_name, dict_new_tags=None, job_type='native'):
    if os.path.getsize(f'{job_path}/CONTCAR') == 0:
        print(f'CHECK: Cannot continue {job_path}: empty CONTCAR')
    if os.path.isfile(f"{job_path}/CENTCAR"):
        if os.path.getsize(f'{job_path}/CENTCAR') == 0:
            print(f'CHECK: Cannot continue {job_path}: empty CENTCAR')
    else:
        if dict_new_tags is not None:
            incar = Incar.from_file(job_path)
            for tag in dict_new_tags:
                incar.update_tag(key=tag, value=dict_new_tags[tag])
            incar.write(job_path)
        num_previous_refines = str(len(glob(f'{job_path}/ref*/')))
        os.system(f"mkdir {job_path}/ref{num_previous_refines}")
        os.system(f"cp {job_path}/* {job_path}/ref{num_previous_refines}")
        os.system(f"cp {job_path}/CONTCAR {job_path}/POSCAR")
        if os.path.isfile(f"{job_path}/CENTCAR"):
            os.system(f"cp {job_path}/CENTCAR {job_path}/POSCAR")
            os.system(f"cp {job_path}/NEWMODECAR {job_path}/MODECAR")
        submit_job(job_path, job_name, job_type)


def check_queue(job_name):
    in_queue = False
    status = None
    with open(f'{path_qstat_list}/qstat_list.txt') as infile:
        qstat_list = infile.readlines()
    for line in qstat_list:
        if line.split()[2] == job_name:
            in_queue = True
            status = line.split()[4]
            break
    return in_queue, status


def rm_vasp_outputs(job_path):
    files_list = [f for f in os.listdir(job_path) if os.path.isfile(os.path.join(job_path, f))]
    save_list = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'MODECAR', name_submission_vasp]
    if os.path.isfile(f'{job_path}/WAVECAR'):
        if os.path.getsize(f'{job_path}/WAVECAR') != 0:
            save_list += ['WAVECAR']
    if os.path.isfile(f'{job_path}/CHGCAR'):
        if os.path.getsize(f'{job_path}/CHGCAR') != 0:
            save_list += ['CHGCAR']
    for file in files_list:
        if file not in save_list:
            os.remove(f"{job_path}/{file}")
