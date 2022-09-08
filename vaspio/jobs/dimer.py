import os
import sys
import numpy as np
import subprocess
from glob import glob

from vaspio.input_files.incar import Incar
from vaspio.input_files.kpoints import Kpoints
from vaspio.variables import *
from vaspio.cluster_data import *


class NewDimerVTST:

    def __init__(self, path, name, incar, kpoints, atoms, outcar_path=None, displacement_vector=None, nsw=500, ediffg=-0.01,
                 potcars_dir=potcars_dir_local, PP_dict=project_PP_dict):
        self.path = path
        self.name = name
        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.potcars_dir = potcars_dir
        self.PP_dict = PP_dict

        incar.add_tag(key='ICHAIN', value=2)
        incar.update_tag(key='IBRION', value=3)
        incar.update_tag(key='POTIM', value=0)
        incar.add_tag(key='IOPT', value=2)

        incar.add_tag(key='NSW', value=nsw)
        incar.add_tag(key='EDIFFG', value=ediffg)

        if displacement_vector is None:
            if outcar_path is None:
                print(f" {self.name}: cannot run Dimer, no displacement vector provided")
            else:
                self.displacement_vector = self.get_displacement_vector_from_outcar(path=outcar_path)
        else:
            self.displacement_vector = displacement_vector

    def create_job_dir(self):
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.atoms.write(filename=f'{self.path}/POSCAR', vasp5=True)
            self.write_potcar()
            self.write_modecar()
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
            cmd += f' {self.potcars_dir}/{self.PP_dict[element]}/POTCAR'
        os.system(f'{cmd} > {self.path}/POTCAR')

    def write_modecar(self):
        with open(f"{self.path}/MODECAR", 'w') as outfile:
            for i in range(len(self.displacement_vector)):
                outfile.write(f"{self.displacement_vector[i][0]} {self.displacement_vector[i][1]} {self.displacement_vector[i][2]}\n")

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

    @staticmethod
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


class DimerVTST:

    def __init__(self, path, name, incar=None, energy=None, status=None):

        self.path = path
        self.name = name
        self.energy = energy
        self.status = status
        self.incar = incar

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

    @classmethod
    def read_from_cluster(cls, path, name):
        """Don't need to read POSCAR and CONTCAR, since they won't be written to JSON"""
        job = cls(path=path, name=name)
        job.incar = Incar.from_file(path=path)
        job.kpoints = Kpoints.from_file(path=path)
        job.status = job.get_job_status()
        if job.status == 'converged':
            job.energy = job.get_energy_oszicar()
        return job

    def converged(self):
        return 'reached required accuracy' in \
           str(subprocess.check_output(f"tail -n4 {self.path}/{name_std_output}", shell=True))

    def nsw_reached(self):
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
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_std_output}", shell=True))
        if 'BAD TERMINATION' in output:
            return True
        else:
            return False

    def other_error(self):
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_std_output}", shell=True))
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
        in_queue, status = check_queue(job_name=self.name)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif os.path.getsize(f'{self.path}/POSCAR') == 0:
            job_status = 'empty poscar'
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif self.converged():
            job_status = 'converged'
        elif self.vasp_bin_not_loaded():
            job_status = 'VASP bin not loaded'
        elif len(glob(f'{self.path}/core.*')) > 0:
            job_status = 'core file'  # see std error
        elif self.bad_termination():
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
        for file in files_list:
            if file not in ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'MODECAR', name_ase_script, name_submission_script]:
                os.remove(f"{self.path}/{file}")

    def submit(self):
        init_dir = os.getcwd()
        os.chdir(self.path)
        self.rm_vasp_outputs()
        if job_scheduler == 'sge':
            os.system(f'qsub -N {self.name} {name_submission_script}')
        elif job_scheduler == 'slurm':
            os.system(f'sbatch --job-name {self.name} {name_submission_script}')
        else:
            sys.exit('Invalid job_scheduler')
        os.chdir(init_dir)

    def restart(self, dict_new_tags=None):
        num_previous_refines = str(len(glob(f'{self.path}/ref*/')))
        if os.path.getsize(f'{self.path}/CENTCAR') == 0:
            print(f'CHECK: Cannot restart {self.name}: empty CENTCAR')
        else:
            os.system(f"mkdir {self.path}/ref{num_previous_refines}")
            os.system(f"cp {self.path}/* {self.path}/ref{num_previous_refines}")
            os.system(f"cp {self.path}/CENTCAR {self.path}/POSCAR")
            os.system(f"cp {self.path}/NEWMODECAR {self.path}/MODECAR")
            if dict_new_tags is not None:
                for tag in dict_new_tags:
                    self.incar.update_tag(key=tag, value=dict_new_tags[tag])
                self.incar.write(self.path)
            self.submit()
