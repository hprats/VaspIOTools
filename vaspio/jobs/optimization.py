import os
import sys
import subprocess
from copy import deepcopy
from glob import glob
import json

import ase.neb
from ase.io import read

from vaspio.input_files.incar import Incar
from vaspio.input_files.kpoints import Kpoints
from vaspio.input_files.potcar import write_potcar

from vaspio.cluster_data import *


class NewJobNative:
    """A class that represents a new VASP job.

    Attributes:
        path (str): The path of the job including the job name. Will be used as the name of the folder.
        incar (vaspio.Incar): An object that represents an INCAR file.
        kpoints (vaspio.Kpoints): An object that represents a KPOINTS file
        atoms (ase.atoms.Atoms): An ASE Atoms object that contains the geometry for the system

    Examples:
        >>> tags = {'IBRION': 1, 'EDIFF': 1E-05, 'EDIFFG': -0.01}
        >>> submission_script = [
        >>>    '#!/bin/bash -l',
        >>>    '#$ -S /bin/bash',
        >>>    'module load vasp/5.4.1/intel-2015',
        >>>    'gerun vasp_std > vasp.out'
        >>> ]
        >>> my_job = NewJobNative(
        >>>    path='/home/test/N2_gas',
        >>>    incar=Incar(tags),
        >>>    kpoints=Kpoints(numbers=[1, 1, 1]),
        >>>    atoms=ase.atoms.Atoms('N2', [(0, 0, 0), (0, 0, 1.12)]
        >>> )
    """

    def __init__(self, path, name, incar, kpoints, atoms, pp_dict, pp_path):
        if isinstance(path, str):
            self.path = path
            self.name = name
            self.incar = incar
            self.kpoints = kpoints
            self.atoms = atoms
            self.pp_dict = pp_dict
            self.pp_path = pp_path

    def create_job_dir(self):
        """Creates a new director and writes there the VASP input files i.e. INCAR, POSCAR, KPOINTS and POTCAR"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.atoms.write(filename=f'{self.path}/POSCAR', vasp5=True)
            write_potcar(job_path=self.path,
                         poscar_elements=self.atoms.get_chemical_symbols(),
                         pp_dict=self.pp_dict,
                         pp_path=self.pp_path)
        else:
            print(f'{self.path} already exists (nothing done)')


class JobNative:
    """A class that represents a VASP job, typically a finished one.

    Attributes:
        path (str): The full path to the job folder.
        energy (float): The total energy of the system.
        status (str): A string describing the status of the job (e.g. qw, r, fine, ...). See get_job_status().

    Examples:
        >>> job = JobNative(path=f'path-to-folder/VC_001##CO2I_0-tM')
        >>> job.read_from_cluster()
        >>> print(job.energy)
    """

    def __init__(self, path, energy=None, status=None, incar=None, magmom=None):

        self.path = path
        self.name = path.split('/')[-1]
        self.energy = energy
        self.status = status
        self.incar = incar
        self.magmom = magmom

    def write_json(self):
        """Write JSON file."""
        new_job = deepcopy(self)
        json_data = json.dumps(new_job.__dict__, default=lambda o: o.__dict__, indent=4, sort_keys=True)  # serialize
        with open(f"{self.path}/{self.name}.json", 'w') as outfile:
            print(json_data, file=outfile)

    @classmethod
    def read_from_json(cls, path):
        """Don't include attributes that will be set by __init__, e.g. adsorbate_config."""
        if len(glob(f"{path}/*.json")) == 0:
            print('JSON file not found')
            return None
        elif len(glob(f"{path}/*.json")) > 1:
            print('More than one JSON file found')
            return None
        else:
            with open(glob(f"{path}/*.json")[0], 'r') as json_file:
                dct = json.loads(json_file.read())
            incar = Incar.from_dict(dct['incar'])
            kpoints = Kpoints.from_dict(dct['kpoints'])
            energy = dct['energy']
            status = dct['status']
            magmom = dct['magmom']

            job = cls(path=path, energy=energy, status=status, incar=incar, magmom=magmom)
            return job

    @classmethod
    def read_from_cluster(cls, path):
        job = cls(path=path)
        job.incar = Incar.from_file(path=path)
        job.kpoints = Kpoints.from_file(path=path)
        job.status = job.get_job_status()
        if job.status == 'converged':
            job.energy = job.get_energy_oszicar()
            job.magmom = job.get_magmom()
        return job

    @classmethod
    def read_converged(cls, path):
        """Read a converged calculation."""
        job = cls(path=path)
        job.incar = Incar.from_file(path=path)
        job.kpoints = Kpoints.from_file(path=path)
        job.atoms_poscar = read(f'{path}/POSCAR')
        if os.path.isfile(f"{path}/CONTCAR"):
            job.atoms_contcar = read(f'{path}/CONTCAR')
        else:
            print(f"{job.name}: no CONTCAR file")
        if os.path.isfile(f"{path}/OSZICAR"):
            job.energy = job.get_energy_oszicar()
            job.magmom = job.get_magmom()
        else:
            print(f"{job.name}: no OSZICAR file")
        return job

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

    def converged(self):
        if 'NSW' not in self.incar.tags:  # is SPE
            return ' 1 F= ' in \
                   str(subprocess.check_output(f"tail -n4 {self.path}/OSZICAR", shell=True))
        else:
            return 'reached required accuracy' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/{name_vasp_std_output}", shell=True))

    def bracketing_error(self):
        return 'bracketing' in \
               str(subprocess.check_output(f"tail -n7 {self.path}/{name_vasp_std_output}", shell=True))

    def get_dE_last_two_steps(self):
        output = str(subprocess.check_output(f"grep F {self.path}/OSZICAR | tail -n2", shell=True))
        last = float(output.split('\\n')[1].split('  d E')[0].split(' ')[-1])
        previous = float(output.split('\\n')[0].split('  d E')[0].split(' ')[-1])
        diff = abs(last - previous)
        return diff

    def nsw_reached(self):
        if 'NSW' not in self.incar.tags:
            return False
        else:
            if os.path.isfile(f"{self.path}/OSZICAR"):
                try:
                    output = str(subprocess.check_output(f"grep F {self.path}/OSZICAR", shell=True))
                    return f"{self.incar.tags['NSW']} F=" in output
                except subprocess.CalledProcessError:
                    return False
            else:
                return False

    def nelm_reached(self):
        if os.path.isfile(f"{self.path}/OSZICAR"):
            nelm = self.incar.tags['NELM']
            try:
                output = str(subprocess.check_output(f"tail -n{int(nelm) + 15} {self.path}/OSZICAR", shell=True))
                return f"RMM: {nelm}" in output
            except subprocess.CalledProcessError:
                return False
        else:
            return False

    def bad_termination(self):
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_vasp_std_output}", shell=True))
        if 'BAD TERMINATION' in output:
            return True
        else:
            return False

    def other_error(self):
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_vasp_std_output}", shell=True))
        if 'error' in output and 'errors must be expected' not in output:
            return True
        else:
            return False

    def vasp_bin_not_loaded(self):
        try:
            if job_scheduler == 'sge':
                std_error_file = glob(f"{self.path}/*.e*")[0]
                return 'cannot be loaded' in str(subprocess.check_output(f"tail {std_error_file}", shell=True))
            elif job_scheduler == 'slurm':
                std_error_file = glob(f"{self.path}/slurm-*.out")[0]
                return 'mpirun: command not found' in str(subprocess.check_output(f"tail {std_error_file}", shell=True))
            else:
                sys.exit('Invalid job_scheduler')
        except subprocess.CalledProcessError:
            return False

    def get_job_status(self):
        """Execute it on the cluster."""
        """If modified here, modify also in vasp_tools.py."""
        in_queue, status = check_queue(job_name=self.name)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif os.path.getsize(f'{self.path}/POSCAR') == 0:
            job_status = 'empty poscar'
        elif not os.path.isfile(f"{self.path}/{name_vasp_std_output}"):
            job_status = 'not submitted'
        elif self.converged():
            job_status = 'converged'
        elif self.vasp_bin_not_loaded():
            job_status = 'VASP bin not loaded'
        elif self.bracketing_error():
            if self.get_dE_last_two_steps() <= 0.01:
                job_status = 'converged'
            else:
                job_status = 'bracketing'
        elif len(glob(f'{self.path}/core.*')) > 0:
            job_status = 'core file'
        elif self.bad_termination():
            if self.get_dE_last_two_steps() <= 0.01:
                job_status = 'converged'
            else:
                job_status = 'bad termination'
        elif self.nsw_reached():
            job_status = 'NSW reached'
        elif self.nelm_reached():
            job_status = 'max wallclock, NELM reached'
        elif not os.path.isfile(f"{self.path}/OSZICAR"):
            job_status = 'error'
        elif self.other_error():
            job_status = 'error'
        else:
            job_status = 'max wallclock'
        return job_status

    def rm_vasp_outputs(self):
        files_list = [f for f in os.listdir(self.path) if os.path.isfile(os.path.join(self.path, f))]
        save_list = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', name_submission_vasp_ase, name_submission_vasp_native]
        if os.path.isfile(f'{self.path}/WAVECAR'):
            if os.path.getsize(f'{self.path}/WAVECAR') != 0:
                save_list += ['WAVECAR']
        if os.path.isfile(f'{self.path}/CHGCAR'):
            if os.path.getsize(f'{self.path}/CHGCAR') != 0:
                save_list += ['CHGCAR']
        for file in files_list:
            if file not in save_list:
                os.remove(f"{self.path}/{file}")

    # todo: use submit_job in functions.py and add new tags ouside this function
    def submit(self, dict_new_tags=None):
        self.rm_vasp_outputs()
        if dict_new_tags is not None:  # e.g. re-submit with ALGO = Normal when NELM is reached
            for tag in dict_new_tags:
                self.incar.update_tag(key=tag, value=dict_new_tags[tag])
            self.incar.write(self.path)
        init_dir = os.getcwd()
        os.chdir(self.path)
        if job_scheduler == 'sge':
            os.system(f'qsub -N {self.name} {name_submission_vasp_native}')
        elif job_scheduler == 'slurm':
            os.system(f'sbatch --job-name {self.name} {name_submission_vasp_native}')
        else:
            sys.exit('Invalid job_scheduler')
        os.chdir(init_dir)

    def restart(self, dict_new_tags=None):
        num_previous_refines = str(len(glob(f'{self.path}/ref*/')))
        if os.path.getsize(f'{self.path}/CONTCAR') == 0:
            print(f'CHECK: Cannot refine {self.name}: empty CONTCAR')
        else:
            os.system(f"mkdir {self.path}/ref{num_previous_refines}")
            os.system(f"cp {self.path}/* {self.path}/ref{num_previous_refines}")
            os.system(f"cp {self.path}/CONTCAR {self.path}/POSCAR")
            # Add new tags if needed
            if dict_new_tags is not None:
                for tag in dict_new_tags:
                    self.incar.update_tag(key=tag, value=dict_new_tags[tag])
            # Start from WAVECAR and CHGCAR if needed
            if os.path.isfile(f'{self.path}/WAVECAR'):
                if os.path.getsize(f'{self.path}/WAVECAR') != 0:
                    self.incar.update_tag(key='ISTART', value=1)
            if os.path.isfile(f'{self.path}/CHGCAR'):
                if os.path.getsize(f'{self.path}/CHGCAR') != 0:
                    self.incar.update_tag(key='ICHARG', value=1)
            self.incar.write(self.path)
            self.submit()
