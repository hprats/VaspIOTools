import os
import sys
import subprocess

from ase.constraints import FixAtoms

from vaspio.input_files.incar import Incar
from vaspio.input_files.kpoints import Kpoints
from vaspio.variables import *

from cluster_data import *


class NewVibrationsNative:

    def __init__(self, path, name, incar, kpoints, atoms, num_free_atoms, potim=0.03, potcars_dir=potcars_dir_local):
        self.path = path
        self.name = name
        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.num_free_atoms = num_free_atoms
        self.potcars_dir = potcars_dir

        self.incar.update_tag(key='IBRION', value='5')
        self.incar.update_tag(key='POTIM', value=potim)
        self.incar.remove_tag(key='NPAR')

        c = FixAtoms(indices=list(range(len(self.atoms) - num_free_atoms)))
        self.atoms.set_constraint(c)

    def create_job_dir(self):
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.atoms.write(filename=f'{self.path}/POSCAR', vasp5=True)
            self.write_potcar()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_potcar(self):
        current_element = self.atoms.get_chemical_symbols()[0]
        elements_for_potcar = [current_element]
        for element in self.atoms.get_chemical_symbols()[1:]:
            if element != current_element:
                elements_for_potcar.append(element)
                current_element = element
        cmd = 'cat'
        for element in elements_for_potcar:
            cmd += f' {self.potcars_dir}/{project_PP_dict[element]}/POTCAR'
        os.system(f'{cmd} > {self.path}/POTCAR')

    def submit(self):
        init_dir = os.getcwd()
        os.chdir(self.path)
        if job_scheduler == 'sge':
            os.system(f'qsub -N {self.name} {name_submission_script}')
        elif job_scheduler == 'slurm':
            os.system(f'sbatch --job-name {self.name} {name_submission_script}')
        else:
            sys.exit('Invalid job_scheduler')
        os.chdir(init_dir)


class VibrationsNative:

    def __init__(self, path, name):

        self.path = path
        self.name = name

    @classmethod
    def read_from_cluster(cls, path, name):
        job = cls(path=path, name=name)
        job.incar = Incar.from_file(path=path)
        job.kpoints = Kpoints.from_file(path=path)
        job.status = job.get_job_status()
        if job.status == 'converged':
            if not os.path.isfile(f"{path}/vibrations.txt"):
                os.system(f"grep cm-1 {path}/OUTCAR > {path}/vibrations.txt")
            job.num_imaginary_frequencies = job.get_num_imaginary_frequencies()
        return job

    def converged(self):
        return 'Voluntary context switches' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/OUTCAR", shell=True))

    def get_num_imaginary_frequencies(self):
        output = str(subprocess.check_output(f"grep cm-1 {self.path}/OUTCAR", shell=True))
        return output.count('f/i')

    def get_job_status(self):
        in_queue, status = check_queue(job_name=self.name)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif self.converged():
            job_status = 'converged'
        else:
            job_status = 'max wallclock or error'
        return job_status

    def submit(self):
        init_dir = os.getcwd()
        os.chdir(self.path)
        if job_scheduler == 'sge':
            os.system(f'qsub -N {self.name} {name_submission_script}')
        elif job_scheduler == 'slurm':
            os.system(f'sbatch --job-name {self.name} {name_submission_script}')
        else:
            sys.exit('Invalid job_scheduler')
        os.chdir(init_dir)
