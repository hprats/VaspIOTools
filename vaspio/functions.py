import os
import sys


# todo: add functions to change memory, walltime limit, number of nodes ... etc
def submit_job(job_path, job_name, job_scheduler, submission_script_path, submission_script_name):
    init_dir = os.getcwd()
    os.chdir(job_path)
    os.system(f"cp {submission_script_path}/{submission_script_name} {job_path}")
    if job_scheduler == 'sge':
        os.system(f'qsub -N {job_name} {submission_script_name}')
    elif job_scheduler == 'slurm':
        os.system(f'sbatch --job-name {job_name} {submission_script_name}')
    else:
        sys.exit('Invalid job_scheduler')
    os.chdir(init_dir)