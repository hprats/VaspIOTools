import os
import subprocess
import numpy as np
from ase import io


def get_energy_oszicar(job_path):
    energy = None
    with open(f"{job_path}/OSZICAR") as infile:
        lines = infile.readlines()
    final = len(lines) - 1
    for i in range(final, 0, -1):
        if ' F= ' in lines[i]:
            energy = float(lines[i].split()[4])
            break
    return energy


def get_magmom(self):
    magmom = 0.0
    if 'ISPIN' in self.incar.tags:
        if int(self.incar.tags['ISPIN']) == 2:
            with open(f"{self.path}/OSZICAR", 'r') as infile:
                lines = infile.readlines()
            final = len(lines) - 1
            for i in range(final, 0, -1):
                if ' F= ' in lines[i]:
                    magmom = float(lines[i].split()[-1])
                    break
    return magmom


def get_num_imaginary_frequencies(self):
    if os.path.isfile(f"{self.path}/vibrations.txt"):
        output = str(subprocess.check_output(f"grep cm-1 {self.path}/vibrations.txt", shell=True))
    else:
        output = str(subprocess.check_output(f"grep cm-1 {self.path}/OUTCAR", shell=True))
    return output.count('f/i')


def get_displacement_vector_from_outcar(path):
    with open(f"{path}/OUTCAR", 'r') as infile:
        lines = infile.readlines()
    # Find start of imaginary vibration output
    line_start = 0
    i = 0
    for i in range(len(lines) - 1, 0, -1):
        if ' f/i= ' in lines[i]:
            line_start = i + 2
            break
    # Get length of imaginary vibration
    line = lines[i + 2]
    num_displacements = 0
    while 'Finite differences POTIM=' not in line:
        num_displacements += 1
        i += 1
        line = lines[i + 2]
    # ignore white line at the end of displacements
    num_displacements -= 1
    # Save vibration to numpy array
    displacement_vector = np.zeros((num_displacements, 3))
    for i in range(num_displacements):
        line = lines[line_start + i]
        displacement_vector[i][0] = line.split()[-3]
        displacement_vector[i][1] = line.split()[-2]
        displacement_vector[i][2] = line.split()[-1]
    return displacement_vector


def get_mlneb_image_energies(job_path):
    if os.path.isfile(f"{job_path}/ML-NEB.traj"):
        images = io.read(f"{job_path}/ML-NEB.traj", index=":")
        image_energies = np.zeros(len(images))
        for i in range(len(images)):
            image_energies[i] = images[i].get_potential_energy()
    elif os.path.isfile(f"{job_path}/last_predicted_path.traj"):
        images = io.read(f"{job_path}/last_predicted_path.traj", index=":")
        image_energies = np.zeros(len(images))
        for i in range(len(images)):
            image_energies[i] = images[i].get_potential_energy()
    else:
        image_energies = None
    return image_energies


def get_mlneb_energy_barrier(image_energies):
    return np.max(image_energies) - image_energies[0]


def get_mlneb_reaction_energy(image_energies):
    return image_energies[-1] - image_energies[0]
