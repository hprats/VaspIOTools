import os
import subprocess
from glob import glob
import numpy as np

from vaspio.cluster_functions import check_queue
from vaspio.read_functions import get_mlneb_image_energies
from vaspio.cluster_data import *


def get_mlneb_status(job_path, job_name):
    """Execute it on the cluster."""
    in_queue, status = check_queue(job_name)
    if in_queue:
        job_status = status
    elif os.path.isfile(f"{job_path}/README"):
        with open(f"{job_path}/README") as f:
            readme_info = f.readlines()[0].strip()
        job_status = f'README: {readme_info}'
    elif len(glob(f'{job_path}/core.*')) > 0:
        job_status = 'core file'  # see std error
    elif not os.path.isfile(f"{job_path}/{name_vasp_std_output}"):
        job_status = 'not submitted'
    elif mlneb_converged(job_path):
        if new_minimum(job_path):
            job_status = 'converged (new minimum)'
        else:
            job_status = 'converged'
    elif get_num_occurrences_timetxt(job_path=job_path, occurrence="final state") == 0:
        job_status = 'initial state'
    elif get_num_occurrences_timetxt(job_path=job_path, occurrence="ML-NEB") == 0:
        job_status = 'final state'
    elif os.path.getsize(f'{job_path}/initial.traj') == 0:
        job_status = 'empty initial.traj'
    elif os.path.getsize(f'{job_path}/final.traj') == 0:
        job_status = 'empty final.traj'
    elif not os.path.isfile(f"{job_path}/last_predicted_path.traj"):
        job_status = 'last_predicted_path.traj missing'
    else:
        job_status = 'max wallclock'
    return job_status


def mlneb_converged(job_path):
    return 'Energy barrier' in \
           str(subprocess.check_output(f"tail -n4 {job_path}/{name_vasp_std_output}", shell=True))


def new_minimum(job_path):
    image_energies = get_mlneb_image_energies(job_path)
    if np.argmax(image_energies) == 0 or np.argmax(image_energies) == len(image_energies):
        return True
    for i in range(1, np.argmax(image_energies)):
        if image_energies[i] + 0.20 < image_energies[0]:
            return True
    for i in range(np.argmax(image_energies) + 1, len(image_energies) - 1):
        if image_energies[i] + 0.20 < image_energies[-1]:
            return True
    return False


def get_num_occurrences_timetxt(job_path, occurrence):
    with open(f"{job_path}/time.txt", "r") as infile:
        lines = infile.readlines()
    num_occurrences = 0
    for line in lines:
        if occurrence in line:
            num_occurrences += 1
    return num_occurrences
