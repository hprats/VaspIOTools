import os
import sys
from glob import glob
from vaspio.incar import Incar


def submit_job(job_path, job_name, job_scheduler, path_submission_vasp_native, name_submission_vasp_native):
    # todo: add functions to change memory, wall-time limit, number of nodes ... etc
    init_dir = os.getcwd()
    os.chdir(job_path)
    os.system(f"cp {path_submission_vasp_native}/{name_submission_vasp_native} {job_path}")
    if job_scheduler == 'sge':
        os.system(f'qsub -N {job_name} {name_submission_vasp_native}')
    elif job_scheduler == 'slurm':
        os.system(f'sbatch --job-name {job_name} {name_submission_vasp_native}')
    else:
        sys.exit('Invalid job_scheduler')
    os.chdir(init_dir)


def continue_job(job_path, job_name, job_scheduler, path_submission_vasp_native, name_submission_vasp_native,
                 dict_new_tags=None):
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
        submit_job(job_path, job_name, job_scheduler, path_submission_vasp_native, name_submission_vasp_native)


def check_queue(job_name, path_qstat_list):
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


def rm_vasp_outputs(self, name_submission_vasp_ase, name_submission_vasp_native):
    files_list = [f for f in os.listdir(self.path) if os.path.isfile(os.path.join(self.path, f))]
    save_list = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'MODECAR', name_submission_vasp_ase, name_submission_vasp_native]
    if os.path.isfile(f'{self.path}/WAVECAR'):
        if os.path.getsize(f'{self.path}/WAVECAR') != 0:
            save_list += ['WAVECAR']
    if os.path.isfile(f'{self.path}/CHGCAR'):
        if os.path.getsize(f'{self.path}/CHGCAR') != 0:
            save_list += ['CHGCAR']
    for file in files_list:
        if file not in save_list:
            os.remove(f"{self.path}/{file}")
