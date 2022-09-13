import os
import sys
import numpy as np
import subprocess
from glob import glob

import ase.neb
from ase import io

from vaspio.variables import *
from vaspio.input_files.incar import Incar
from vaspio.input_files.kpoints import Kpoints
from vaspio.cluster_data import *  # If executed locally, this could be an empty file


class NewNebNative:

    def __init__(self, path, images, incar, kpoints, atoms_initial, atoms_final,
                 energy_initial, energy_final, potcars_dir=potcars_dir_local, PP_dict=project_PP_dict):
        if isinstance(path, str):
            self.path = path
            self.images = images
            self.name = path.split('/')[-1]
            self.incar = incar
            self.kpoints = kpoints
            self.atoms_initial = atoms_initial
            self.atoms_final = atoms_final
            self.energy_initial = energy_initial
            self.energy_final = energy_final
            self.potcars_dir = potcars_dir
            self.PP_dict = PP_dict

    def create_job_dir(self):
        """Creates a new directory, with the name job_name, and writes there the VASP input files
        i.e. INCAR, POSCAR, KPOINTS and POTCAR, and the submission script"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            # Write input files
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.write_potcar()
            self.write_energies()
            # Create all folders
            for i in range(self.images + 2):
                if i <= 9:
                    os.mkdir(f'{self.path}/0{i}')
                else:
                    os.mkdir(f'{self.path}/{i}')
            # Write POSCAR files for reactants and products
            self.atoms_initial.write(filename=f'{self.path}/00/POSCAR', vasp5=True)
            if self.images <= 9:
                self.atoms_final.write(filename=f'{self.path}/0{self.images + 1}/POSCAR', vasp5=True)
            else:
                self.atoms_final.write(filename=f'{self.path}/{self.images + 1}/POSCAR', vasp5=True)
            # Write POSCAR files for images
            constraints = self.atoms_initial.constraints
            list_images = [self.atoms_initial]
            for i in range(self.images):
                image = self.atoms_initial.copy()
                image.set_constraint(constraints)
                list_images.append(image)
            list_images.append(self.atoms_final)
            neb = ase.neb.NEB(list_images, climb=False, k=0.5)
            neb.interpolate('idpp')
            for i in range(1, self.images + 1):
                if i <= 9:
                    list_images[i].write(filename=f'{self.path}/0{i}/POSCAR', vasp5=True)
                else:
                    list_images[i].write(filename=f'{self.path}/{i}/POSCAR', vasp5=True)
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_energies(self):
        with open(f"{self.path}/energies.txt", 'w') as outfile:
            outfile.write(f'{self.energy_initial}\n')
            outfile.write(f'{self.energy_final}\n')

    def write_potcar(self):
        """Writes the POTCAR file to the current directory
        First, a list of elements is created (elements_for_potcar), which determines which
        atomic POTCAR files will be concatenated to make the final POTCAR"""
        current_element = self.atoms_initial.get_chemical_symbols()[0]
        elements_for_potcar = [current_element]
        for element in self.atoms_initial.get_chemical_symbols()[1:]:
            if element != current_element:
                elements_for_potcar.append(element)
                current_element = element
        cmd = 'cat'
        for element in elements_for_potcar:
            cmd += f' {self.potcars_dir}/{self.PP_dict[element]}/POTCAR'
        os.system(f'{cmd} > {self.path}/POTCAR')


class NewNebML:

    def __init__(self, path, n_images, fmax, incar, kpoints, atoms_initial, atoms_final):
        self.path = path
        self.n_images = n_images
        self.fmax = fmax
        self.name = path.split('/')[-1]
        incar.add_tag(key='n_images', value=n_images)
        incar.add_tag(key='fmax', value=fmax)
        self.incar = incar
        self.kpoints = kpoints
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final

    def create_job_dir(self):
        """Creates a new directory, with the name job_name, and writes there the VASP input files
        i.e. INCAR, POSCAR, KPOINTS and POTCAR, and the submission script"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.write_trajectory_files()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_trajectory_files(self):
        os.mkdir(f"{self.path}/optimized_structures")
        io.write(f"{self.path}/optimized_structures/initial.traj", self.atoms_initial)
        io.write(f"{self.path}/optimized_structures/final.traj", self.atoms_final)
        list_images = [self.atoms_initial]
        constraints = self.atoms_initial.constraints
        for i in range(self.n_images - 2):
            image = self.atoms_initial.copy()
            image.set_constraint(constraints)
            list_images.append(image)
        list_images.append(self.atoms_final)
        neb = ase.neb.NEB(list_images, climb=False, k=0.5)
        neb.interpolate('idpp')
        io.write(f"{self.path}/optimized_structures/list_images.traj", list_images)


class NebML:

    def __init__(self, path):
        self.path = path
        self.name = path.split('/')[-1]
        self.incar = None
        self.kpoints = None
        self.energies = None
        self.status = None

    def get_energies(self):
        if os.path.isfile(f"{self.path}/ML-NEB.traj"):
            images = io.read(f"{self.path}/ML-NEB.traj", index=":")
            energies = np.zeros(len(images))
            for i in range(len(images)):
                energies[i] = images[i].get_potential_energy()
        elif os.path.isfile(f"{self.path}/last_predicted_path.traj"):
            images = io.read(f"{self.path}/last_predicted_path.traj", index=":")
            energies = np.zeros(len(images))
            for i in range(len(images)):
                energies[i] = images[i].get_potential_energy()
        else:
            energies = None
        return energies

    def get_energy_barrier(self):
        if self.energies is None:
            return None
        else:
            return np.max(self.energies) - self.energies[0]

    def get_reaction_energy(self):
        return self.energies[-1] - self.energies[0]

    @classmethod
    def read_from_cluster(cls, path):
        job = cls(path=path)
        if os.path.isfile(f"{job.path}/INCAR"):
            job.incar = Incar.from_file(path=path)
        if os.path.isfile(f"{job.path}/KPOINTS"):
            job.kpoints = Kpoints.from_file(path=path)
        job.energies = job.get_energies()
        job.status = job.get_job_status()
        job.energy_barrier = job.get_energy_barrier()
        return job

    def converged(self):
        return 'Energy barrier' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/{name_std_output}", shell=True))

    def new_minimum(self):
        if np.argmax(self.energies) == 0 or np.argmax(self.energies) == len(self.energies):
            return True
        for i in range(1, np.argmax(self.energies)):
            if self.energies[i] + 0.20 < self.energies[0]:
                return True
        for i in range(np.argmax(self.energies) + 1, len(self.energies) - 1):
            if self.energies[i] + 0.20 < self.energies[-1]:
                return True
        return False

    def get_num_occurrences_timetxt(self, occurrence="ML-NEB"):
        with open(f"{self.path}/time.txt", "r") as infile:
            lines = infile.readlines()
        num_occurrences = 0
        for line in lines:
            if occurrence in line:
                num_occurrences += 1
        return num_occurrences

    def get_job_status(self):
        """Execute it on the cluster."""
        in_queue, status = check_queue(job_name=self.name)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif len(glob(f'{self.path}/core.*')) > 0:
            job_status = 'core file'  # see std error
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif self.converged():
            if self.new_minimum():
                job_status = 'converged (new minimum)'
            else:
                job_status = 'converged'
        elif self.get_num_occurrences_timetxt(occurrence="final state") == 0:
            job_status = 'initial state'
        elif self.get_num_occurrences_timetxt(occurrence="ML-NEB") == 0:
            job_status = 'final state'
        elif os.path.getsize(f'{self.path}/initial.traj') == 0:
            job_status = 'empty initial.traj'
        elif os.path.getsize(f'{self.path}/final.traj') == 0:
            job_status = 'empty final.traj'
        else:
            job_status = 'max wallclock'
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

    @staticmethod
    def get_import_lines():
        lines = [
            "from ase.io import read\n",
            "from ase.optimize import BFGS\n",
            "from ase.calculators.vasp import Vasp\n",
            "import shutil\n",
            "import copy\n",
            "from catlearn.optimize.mlneb import MLNEB\n",
            "from datetime import datetime\n",
            " \n"
        ]
        return lines

    @staticmethod
    def get_calculator_lines(lines):
        index_start = 0
        index_end = 0
        for i, line in enumerate(lines):
            if "ase_calculator = " in line:
                index_start = i
            if "# Kpoints" in line:
                index_end = i + 5
        return lines[index_start:index_end]

    @staticmethod
    def get_initial_lines(lines):
        index_fmax = lines.index("shutil.copy('./initial.traj', './optimized_structures/initial.traj')\n") - 1
        fmax = lines[index_fmax].split("=")[1].split(")")[0]
        lines = [
            "# Optimize initial state:\n",
            "now = datetime.now()\n",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
            "with open('time.txt', 'a') as outfile:\n",
            "\toutfile.write(f'{dt_string}: starting optimization of initial state ...\\n')\n",
            "slab = read('./initial.traj')\n",
            "slab.set_calculator(copy.deepcopy(ase_calculator))\n",
            "qn = BFGS(slab, trajectory='initial.traj')\n",
            f"qn.run(fmax={fmax})\n",
            "shutil.copy('./initial.traj', './optimized_structures/initial.traj')\n",
            " \n"
        ]
        return lines

    @staticmethod
    def get_final_lines(lines, restart=False):
        index_fmax = lines.index("shutil.copy('./final.traj', './optimized_structures/final.traj')\n") - 1
        fmax = lines[index_fmax].split("=")[1].split(")")[0]
        lines = [
            "# Optimize final state:\n",
            "now = datetime.now()\n",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
            "with open('time.txt', 'a') as outfile:\n",
            "\toutfile.write(f'{dt_string}: starting optimization of final state ...\\n')\n"
        ]
        if restart:
            lines += ["slab = read('./final.traj')\n"]
        else:
            lines += ["slab = read('./optimized_structures/final.traj')\n"]
        lines += [
            "slab.set_calculator(copy.deepcopy(ase_calculator))\n",
            "qn = BFGS(slab, trajectory='final.traj')\n",
            f"qn.run(fmax={fmax})\n",
            "shutil.copy('./final.traj', './optimized_structures/final.traj')\n",
            " \n"
        ]
        return lines

    @staticmethod
    def get_MLneb_lines(lines, restart=False):
        index_n_images = lines.index("neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n") + 2
        n_images = lines[index_n_images].split("=")[1].split(",")[0]
        index_fmax = lines.index("neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n") + 6
        fmax = lines[index_fmax].split("=")[1].split(",")[0]
        lines = [
            "# Run ML-NEB:\n",
            "now = datetime.now()\n",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
            "with open('time.txt', 'a') as outfile:\n",
            "\toutfile.write(f'{dt_string}: starting ML-NEB ...\\n')\n",
            "neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n",
            "\tase_calc=copy.deepcopy(ase_calculator),\n",
            f"\tn_images={n_images},\n",
            "\tinterpolation='idpp',\n"
        ]
        if restart:
            lines += ["\trestart=True,\n"]
        else:
            lines += ["\t#restart=True,\n"]
        lines += [
            "\t)\n",
            f"neb_catlearn.run(fmax={fmax}, trajectory='ML-NEB.traj')\n",
            " \n"
        ]
        return lines

    @staticmethod
    def get_print_lines():
        lines = [
            "# Print results:\n",
            "now = datetime.now()\n",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
            "with open('time.txt', 'a') as outfile:\n",
            "\toutfile.write(f'{dt_string}: ML-NEB finished\\n')\n",
            " \n",
            "from catlearn.optimize.tools import plotneb\n",
            "plotneb(trajectory='ML-NEB.traj', view_path=True)\n"
        ]
        return lines

    def restart(self):
        if self.status == "max wallclock" and self.get_num_occurrences_timetxt(occurrence="ML-NEB") > 1:
            pass
        else:
            with open(f"{self.path}/run.py", "r") as infile:
                lines = infile.readlines()
            ase_script = self.get_import_lines()
            ase_script += self.get_calculator_lines(lines=lines)
            if self.status == "max wallclock":
                ase_script += self.get_MLneb_lines(lines=lines, restart=True)
                ase_script += self.get_print_lines()
            elif self.status == "final state":
                ase_script += self.get_final_lines(lines=lines, restart=True)
                ase_script += self.get_MLneb_lines(lines=lines, restart=False)
                ase_script += self.get_print_lines()
            elif self.status == "initial state":
                ase_script += self.get_initial_lines(lines=lines)
                ase_script += self.get_final_lines(lines=lines, restart=False)
                ase_script += self.get_MLneb_lines(lines=lines, restart=False)
                ase_script += self.get_print_lines()
            else:
                sys.exit(f"{self.name}: ERROR, cannot determine where to restart the job, check manually")
            with open(f"{self.path}/{name_ase_script}", 'w') as outfile:
                for line in ase_script:
                    outfile.write(line)
        self.submit()

    def write_ase_script(self):
        ase_script = [
            # Imports
            "from ase.io import read",
            "from ase.optimize import BFGS",
            "from ase.calculators.vasp import Vasp",
            "import shutil",
            "import copy",
            "from catlearn.optimize.mlneb import MLNEB",
            "from datetime import datetime",
            " ",
            # Initial time
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'w') as outfile:",
            "\toutfile.write(f'{dt_string}: starting job ...\\n')",
            " ",
            # Ase calculator
            "ase_calculator = Vasp(setups={'base': 'recommended', 'W': '_pv'},"  # todo: udpate so that W_pw is not used by default
        ]
        for tag in [tag for tag in self.incar.tags if tag not in ['ediffg', 'n_images', 'fmax']]:
            ase_script.append(f"\t{tag}={self.incar.tags[tag]},")
        ase_script += ["\t# Kpoints",
                       f"\tkpts=({self.kpoints.num_x}, {self.kpoints.num_y}, {self.kpoints.num_z}),",
                       f"\tgamma=True,",
                       "\t)",
                       " ",
                       ]
        ase_script += [
            # Initial state optimization
            "# Optimize initial state:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting optimization of initial state ...\\n')",
            "slab = read('./optimized_structures/initial.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='initial.traj')",
            f"qn.run(fmax={abs(float(self.incar.tags['ediffg']))})",
            "shutil.copy('./initial.traj', './optimized_structures/initial.traj')",
            " ",
            # Final state optimization
            "# Optimize final state:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting optimization of final state ...\\n')",
            "slab = read('./optimized_structures/final.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='final.traj')",
            f"qn.run(fmax={abs(float(self.incar.tags['ediffg']))})",
            "shutil.copy('./final.traj', './optimized_structures/final.traj')",
            " ",
            # NEB
            "# Run ML-NEB:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting ML-NEB ...\\n')",
            "neb_catlearn = MLNEB(start='initial.traj', end='final.traj',",
            "\tase_calc=copy.deepcopy(ase_calculator),",
            f"\tn_images={self.incar.tags['n_images']},",
            "\tinterpolation='idpp',",
            "\t#restart=True,",
            "\t)",
            f"neb_catlearn.run(fmax={self.incar.tags['fmax']}, trajectory='ML-NEB.traj')",
            " ",
            # Print results
            "# Print results:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: ML-NEB finished\\n')",
            " ",
            "from catlearn.optimize.tools import plotneb",
            "plotneb(trajectory='ML-NEB.traj', view_path=True)"
        ]
        with open(f"{self.path}/{name_ase_script}", 'w') as outfile:
            for line in ase_script:
                outfile.write(line + '\n')
