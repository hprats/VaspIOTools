import os
import subprocess

import ase
from ase.io import read as ase_read

from vaspio.incar import Incar
from vaspio.kpoints import Kpoints


class NewJob:
    """A class that represents a new VASP job.

    Attributes:
        job_name (str): The name of the job. Will be used as the name of the folder.
        incar (vaspio.Incar): An object that represents an INCAR file.
        kpoints (vaspio.Kpoints): An object that represents a KPOINTS file
        atoms (ase.atoms.Atoms): An ASE Atoms object that contains the geometry for the system
        submission_script (list): A list of strings. Each element represents a line of the submission script.

    Examples:
        >>> tags = {'IBRION': 1, 'EDIFF': 1E-05, 'EDIFFG': -0.01}
        >>> submission_script = [
        >>>    '#!/bin/bash -l',
        >>>    '#$ -S /bin/bash',
        >>>    'module load vasp/5.4.1/intel-2015',
        >>>    'gerun vasp_std > vasp.out'
        >>> ]
        >>> my_job = NewJob(
        >>>    job_name='N2_gas',
        >>>    incar=Incar(tags),
        >>>    kpoints=Kpoints(numbers=[1, 1, 1]),
        >>>    atoms=ase.atoms.Atoms('N2', [(0, 0, 0), (0, 0, 1.12)],
        >>>    submission_script=[submission_script]
        >>> )
    """

    POTCARS_DIR = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potcars'
    NAME_SUBMISSION_SCRIPT = 'vasp_sub'

    def __init__(self, job_name, incar, kpoints, atoms, submission_script):
        if isinstance(job_name, str):
            self.job_name = job_name
            self.incar = incar
            self.kpoints = kpoints
            self.atoms = atoms
            self.submission_script = submission_script

    def write_submission_script(self):
        """Prints the submission script into a new file named vasp_sub"""
        f = open(f'{NewJob.NAME_SUBMISSION_SCRIPT}', 'w')
        for line in self.submission_script:
            f.write(line + '\n')
        f.close()

    def create_job_dir(self, path='.'):
        """Creates a new directory, with the name job_name, and writes there the VASP input files
        i.e. INCAR, POSCAR, KPOINTS and POTCAR, and the submission script"""
        initial_directory = os.getcwd()
        os.chdir(path)
        os.mkdir(self.job_name)
        os.chdir(self.job_name)
        self.incar.write()
        self.kpoints.write()
        self.atoms.ase_write(filename='POSCAR', vasp5=True)
        self.write_potcar(self.atoms)
        self.write_submission_script()
        os.chdir(initial_directory)

    @staticmethod
    def write_potcar(atoms):
        """Writes the POTCAR file to the current directory
        First, a list of elements is created (elements_for_potcar), which determines which
        atomic POTCAR files will be concatenated to make the full POTCAR"""
        elements_for_potcar = []
        current_element = atoms.get_chemical_symbols()[0]
        for element in atoms.get_chemical_symbols()[1:]:
            if element != current_element:
                elements_for_potcar.append(element)
                current_element = element
        with open('POTCAR', 'w') as outfile:
            for element in elements_for_potcar:
                with open(f"{NewJob.POTCARS_DIR}/POTCAR-{element}") as infile:
                    for line in infile:
                        outfile.write(line)


class Job:
    """A class that represents a VASP job, typically a finished one.

    Attributes:
        job_name (str): The name of the folder containing the job files.
        path (str): The full path to the job folder.
        incar (vaspio.Incar): An object that represents an INCAR file.
        kpoints (vaspio.Kpoints): An object that represents a KPOINTS file.
        atoms_poscar (ase.atoms.Atoms): An ASE Atoms object that contains the geometry of the POSCAR file.
        atoms_contcar (ase.atoms.Atoms): An ASE Atoms object that contains the geometry of the CONTCAR file.
        energy (float): The total energy of the system.
        num_restarts (int): The number of times that the job has been restarted (due to CPU time limit).
        status (str): A string describing the status of the job (e.g. qw, r, fine, ...). See get_job_status().
        name_std_output (str): The name of the VASP standard output file.

    Examples:
        >>> job = Job(job_name='Ni_111#CO-top', path=f'path-to-folder/VC_001##CO2I_0-tM')
        >>> job.read()
        >>> print(job.energy)
    """

    POTCARS_DIR = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potcars'

    def __init__(self, path,
                 incar=None, kpoints=None, atoms_poscar=None, atoms_contcar=None,
                 energy=None, num_restarts=None, status=None, name_std_output='vasp.out'):

        self.job_name = path.split('/')[-1]
        self.path = path
        self.incar = incar
        self.kpoints = kpoints
        self.atoms_poscar = atoms_poscar
        self.atoms_contcar = atoms_contcar
        self.energy = energy
        self.num_restarts = num_restarts
        self.status = status
        self.name_std_output = name_std_output

    def read(self):
        """Read INCAR, KPOINTS, POSCAR and CONTCAR files and store their information"""
        if os.path.isfile(f'{self.path}/INCAR'):
            self.incar = self.read_incar()
        if os.path.isfile(f'{self.path}/KPOINTS'):
            self.kpoints = self.read_kpoints()
        if os.path.isfile(f'{self.path}/POSCAR'):
            self.atoms_poscar = ase_read('POSCAR')
        if os.path.isfile(f'{self.path}/CONTCAR'):
            self.atoms_contcar = ase_read('CONTCAR')
        self.status = self.get_job_status()
        if self.status == 'fine':
            self.energy = self.get_energy_oszicar()

    def read_incar(self):
        """Read INCAR file and store its information as a vaspio.Incar object"""
        with open(f'{self.path}/INCAR') as infile:
            lines = infile.readlines()
        incar_tags = {}
        for line in lines:
            tag = line.strip().split(' = ')[0]
            value = line.strip().split(' = ')[1]
            incar_tags[tag] = value
        return Incar(incar_tags)

    def read_kpoints(self):
        """Read KPOINTS file and store its information as a vaspio.Kpoints object"""
        with open(f'{self.path}/KPOINTS') as infile:
            lines = infile.readlines()
        num_x = lines[3].split(' ')[0]
        num_y = lines[3].split(' ')[1]
        num_z = lines[3].split(' ')[2]
        if num_x != 1 and num_y != 1 and num_z == 1:
            is_slab = True
        else:
            is_slab = False
        return Kpoints(is_slab=is_slab, numbers=[num_x, num_y, num_z])

    def converged(self):
        return 'reached required accuracy' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/{self.name_std_output}", shell=True))

    def bracketing_error(self):
        return 'fatal error in bracketing' in \
               str(subprocess.check_output(f"tail -n7 {self.path}/{self.name_std_output}", shell=True))

    def get_bracketing_diff(self):
        output = str(subprocess.check_output(f"grep F {self.path}/OSZICAR | tail -n2", shell=True))
        last = float(output.split('\\n')[1].split('  d E')[0].split(' ')[-1])
        previous = float(output.split('\\n')[0].split('  d E')[0].split(' ')[-1])
        diff = abs(last - previous)
        return diff

    def nsw_reached(self):
        try:
            output = str(subprocess.check_output(f"grep F {self.path}/OSZICAR", shell=True))
            return f"{self.incar.tags['NSW']} F=" in output
        except:
            return False

    def nelm_reached(self):
        nelm = self.incar.tags['NELM']
        return f"RMM: {nelm}" in \
               str(subprocess.check_output(f"tail -n{int(nelm) + 15} {self.path}/OSZICAR", shell=True))

    def other_error(self):
        return 'error' in \
               str(subprocess.check_output(f"tail -n10 {self.path}/{self.name_std_output}", shell=True))

    def get_job_status(self):
        """Check the status of the job."""
        if not os.path.isfile(f"{self.path}/{self.name_std_output}"):
            job_status = 'qw'
        elif os.stat(f"{self.path}/{self.name_std_output}").st_size == 0:  # standard output is empty
            job_status = 'VASP bin not loaded'
        elif self.converged():
            job_status = 'fine'
        else:  # not converged or maybe still running, check
            if self.bracketing_error():
                if self.get_bracketing_diff() <= 0.01:
                    job_status = 'fine'
                else:
                    job_status = 'bracketing'
            elif self.nsw_reached():
                job_status = 'NSW reached'
            elif self.nelm_reached():
                job_status = 'NELM reached (CPU time)'
            elif self.other_error():
                job_status = 'error'
            else:
                job_status = 'not converged (CPU time) or still running'
        return job_status

    def get_energy_oszicar(self):
        energy = None
        with open(f"{self.path}/OSZICAR") as infile:
            lines = infile.readlines()
        final = len(lines) - 1
        for i in range(final, 0, -1):
            if ' F= ' in lines[i]:
                energy = float(lines[i].split()[4])
                break
        return energy
