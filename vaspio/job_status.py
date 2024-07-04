import os
import sys
import subprocess
from glob import glob

from vaspio.incar import Incar
from vaspio.cluster_functions import check_queue
from vaspio.cluster_data import *


def get_job_status(job_path, job_name):
    in_queue, status = check_queue(job_name)
    if in_queue:
        job_status = status
    elif os.path.isfile(f"{job_path}/README"):
        with open(f"{job_path}/README") as f:
            readme_info = f.readlines()[0].strip()
        job_status = f'README: {readme_info}'
    elif os.path.getsize(f'{job_path}/POSCAR') == 0:
        job_status = 'empty POSCAR'
    elif not os.path.isfile(f"{job_path}/{name_vasp_std_output}"):
        job_status = 'not submitted'
    elif converged(job_path):
        job_status = 'converged'
    elif vasp_bin_not_loaded(job_path):
        job_status = 'VASP binary file not loaded'
    elif len(glob(f'{job_path}/core.*')) > 0:
        job_status = 'core file generated'
    else:
        incar = Incar.from_file(job_path)
        if incar.tags['IBRION'] == '5':
            job_status = 'max wallclock or unclassified error'
        else:
            if bracketing_error(job_path):
                if get_dE_last_two_steps(job_path) <= 0.01:
                    job_status = 'converged'
                else:
                    job_status = 'bracketing'
            elif bad_termination(job_path):
                if get_dE_last_two_steps(job_path) <= 0.01:
                    job_status = 'converged'
                else:
                    job_status = 'bad termination'
            elif nsw_reached(job_path):
                job_status = 'NSW reached'
            elif nelm_reached(job_path):
                job_status = 'max wallclock, NELM reached'
            elif other_error(job_path):
                job_status = 'unclassified error'
            else:
                job_status = 'max wallclock'
    return job_status


def converged(job_path):
    incar = Incar.from_file(job_path)
    if 'NSW' not in incar.tags:  # is SPE
        return ' 1 F= ' in \
               str(subprocess.check_output(f"tail -n4 {job_path}/OSZICAR", shell=True))
    elif incar.tags['IBRION'] == '5':  # is vibrational analysis
        return 'Voluntary context switches' in \
               str(subprocess.check_output(f"tail -n4 {job_path}/OUTCAR", shell=True))
    else:
        return 'reached required accuracy' in \
           str(subprocess.check_output(f"tail -n4 {job_path}/{name_vasp_std_output}", shell=True))


def vasp_bin_not_loaded(job_path):
    try:
        if job_scheduler == 'sge':
            std_error_file = glob(f"{job_path}/*.e*")[0]
            return 'cannot be loaded' in str(subprocess.check_output(f"tail {std_error_file}", shell=True))
        elif job_scheduler == 'slurm':
            std_error_file = glob(f"{job_path}/slurm-*.out")[0]
            return 'mpirun: command not found' in str(subprocess.check_output(f"tail {std_error_file}", shell=True))
        else:
            sys.exit('Invalid job_scheduler')
    except subprocess.CalledProcessError:
        return False


def bracketing_error(job_path):
    return 'bracketing' in \
           str(subprocess.check_output(f"tail -n7 {job_path}/{name_vasp_std_output}", shell=True))


def get_dE_last_two_steps(job_path):
    output = str(subprocess.check_output(f"grep F {job_path}/OSZICAR | tail -n2", shell=True))
    last = float(output.split('\\n')[1].split('  d E')[0].split(' ')[-1])
    previous = float(output.split('\\n')[0].split('  d E')[0].split(' ')[-1])
    diff = abs(last - previous)
    return diff


def bad_termination(job_path):
    output = str(subprocess.check_output(f"tail -n10 {job_path}/{name_vasp_std_output}", shell=True))
    if 'BAD TERMINATION' in output:
        return True
    else:
        return False


def nsw_reached(job_path):
    incar = Incar.from_file(job_path)
    if 'NSW' not in incar.tags:
        return False
    else:
        if os.path.isfile(f"{job_path}/OSZICAR"):
            try:
                output = str(subprocess.check_output(f"grep F {job_path}/OSZICAR", shell=True))
                return f"{incar.tags['NSW']} F=" in output
            except subprocess.CalledProcessError:
                return False
        else:
            return False


def nelm_reached(job_path):
    if os.path.isfile(f"{job_path}/OSZICAR"):
        incar = Incar.from_file(job_path)
        nelm = incar.tags['NELM']
        try:
            output = str(subprocess.check_output(f"tail -n{int(nelm) + 15} {job_path}/OSZICAR", shell=True))
            return f"RMM: {nelm}" in output
        except subprocess.CalledProcessError:
            return False
    else:
        return False


def other_error(job_path):
    output = str(subprocess.check_output(f"tail -n10 {job_path}/{name_vasp_std_output}", shell=True))
    if 'error' in output and 'errors must be expected' not in output:
        return True
    else:
        return False
