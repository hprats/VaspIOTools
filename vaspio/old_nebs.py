import os
import sys
import numpy as np
import subprocess
from glob import glob
import json

import ase.neb
from ase import io
from ase.constraints import FixAtoms
from pymatgen.io.ase import AseAtomsAdaptor

from vaspio.variables import *
from vaspio.incar import Incar
from vaspio.kpoints import Kpoints
from vaspio.jobs import NewJobNative

from cluster_data import *  # If executed locally, this could be an empty file


class NewNebNative:

    def __init__(self, path, images, incar, kpoints, atoms_initial, atoms_final,
                 energy_initial, energy_final):
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
            cmd += f' {potcars_dir_local}/{project_PP_dict[element]}/POTCAR'
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

    def __init__(self, path, incar=None, kpoints=None, status=None, energies=None, vibrations=None,
                 num_atoms_adsorbate=None):
        self.path = path
        self.name = path.split('/')[-1]
        self.incar = incar
        self.kpoints = kpoints
        self.status = status
        self.energies = energies
        self.vibrations = vibrations
        self.num_atoms_adsorbate = num_atoms_adsorbate

    def __len__(self):
        return len(self.energies)

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
            status = dct['status']
            energies = dct['energies']
            vibrations = dct['vibrations']
            num_atoms_adsorbate = dct['num_atoms_adsorbate']

            job = cls(path=path, status=status, energies=energies, vibrations=vibrations,
                      num_atoms_adsorbate=num_atoms_adsorbate)
            return job

    def get_energies(self):
        if os.path.isfile(f"{self.path}/ML-NEB.traj"):
            images = io.read(f"{self.path}/ML-NEB.traj", index=":")
            energies = np.zeros(len(images))
            for i in range(len(images)):
                energies[i] = images[i].get_potential_energy()
        else:
            energies = None
        return energies

    def get_vibrations(self):
        pass

    def get_energy_barrier(self):
        return np.max(self.energies) - self.energies[0]

    def get_reaction_energy(self):
        return self.energies[-1] - self.energies[0]

    @classmethod
    def read_from_cluster(cls, path, num_atoms_adsorbate):
        job = cls(path=path, num_atoms_adsorbate=num_atoms_adsorbate)
        if os.path.isfile(f"{job.path}/INCAR"):
            job.incar = Incar.from_file(path=path)
        if os.path.isfile(f"{job.path}/KPOINTS"):
            job.kpoints = Kpoints.from_file(path=path)
        job.energies = job.get_energies()
        job.status = job.get_job_status()
        return job

    def NEB_converged(self):
        return 'Energy barrier' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/{name_std_output}", shell=True))

    def freq_done(self):
        return 'Voluntary context switches' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/vibrations/OUTCAR", shell=True))

    def is_ts(self):
        output = str(subprocess.check_output(f"grep cm-1 {self.path}/vibrations/OUTCAR", shell=True))
        return output.count('f/i') == 1

    def empty_file(self):
        if os.path.isfile(f"{self.path}/initial.traj"):
            try:
                with open(f"{self.path}/initial.traj", "r") as infile:
                    lines = infile.readlines()
                return True
            except UnicodeDecodeError:
                return False
        if os.path.isfile(f"{self.path}/final.traj"):
            try:
                with open(f"{self.path}/final.traj", "r") as infile:
                    lines = infile.readlines()
                return True
            except UnicodeDecodeError:
                return False
        return False

    def new_minimum(self):
        if np.argmax(self.energies) == 0 or np.argmax(self.energies) == len(self):
            return True
        for i in range(1, np.argmax(self.energies)):
            if self.energies[i] + 0.20 < self.energies[0]:
                return True
        for i in range(np.argmax(self.energies) + 1, -1):
            if self.energies[i] + 0.20 < self.energies[-1]:
                return True
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
        elif len(glob(f'{self.path}/core.*')) > 0:
            job_status = 'core (see std error)'
        elif self.empty_file():
            job_status = 'empty file'
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif os.path.isdir(f"{self.path}/vibrations"):
            if self.freq_done():
                if self.is_ts():
                    if self.new_minimum():
                        job_status = 'done (new minimum)'
                    else:
                        job_status = 'done'
                else:
                    job_status = 'not TS'
            else:
                job_status = 'vibrations not converged'
        elif self.NEB_converged():
            job_status = 'NEB converged'
        elif not os.path.isfile(f"{self.path}/results_neb.csv"):
            job_status = 'results_neb.csv not found; maybe initial.traj not optimized?'
        else:
            job_status = 'max wallclock'
        return job_status

    @staticmethod
    def submit(path, name):
        init_dir = os.getcwd()
        os.chdir(path)
        if job_scheduler == 'sge':
            os.system(f'qsub -N {name} {name_submission_script}')
        elif job_scheduler == 'slurm':
            os.system(f'sbatch --job-name {name} {name_submission_script}')
        else:
            sys.exit('Invalid job_scheduler')
        os.chdir(init_dir)

    def restart(self):
        """ Continuation job"""
        # Check if it was already a continuation job
        continuation_job = True
        with open(f"{self.path}/{name_ase_script}", "r") as infile:
            lines = infile.readlines()
        for line in lines:
            if "Optimize initial state" in line:
                continuation_job = False
        if not continuation_job:
            # Modify ase file
            ase_script_lines = []
            for i in range(len(lines)):
                if 'now = datetime' in lines[i]:
                    break
                else:
                    ase_script_lines.append(lines[i])
            ase_script_lines += ["now = datetime.now()\n",
                                 "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
                                 "with open('time.txt', 'a') as outfile:\n",
                                 "\toutfile.write(f'{dt_string}: restarting job ...\\n')\n",
                                 " \n"]
            for i in range(len(lines)):
                if 'ase_calculator' in lines[i]:
                    while 'Optimize initial state' not in lines[i]:
                        ase_script_lines.append(lines[i])
                        i += 1
                    break
            for i in range(len(lines)):
                if 'neb_catlearn' in lines[i]:
                    for line in lines[i:]:
                        ase_script_lines.append(line)
                    break
            for i in range(len(ase_script_lines)):
                if "\t#restart=True," in ase_script_lines[i]:
                    ase_script_lines[i] = "\trestart=True" + "\n"
                    break
            with open(f"{self.path}/{name_ase_script}", 'w') as outfile:
                for line in ase_script_lines:
                    outfile.write(line)
        self.submit(path=self.path, name=self.name)

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
            "ase_calculator = Vasp(setups={'base': 'recommended', 'W': '_pv'},"
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

    def run_vibrations_ase(self):
        # todo: not tested
        ts = io.read(f"{self.path}/ML-NEB.traj", index=np.argmax(self.energies))
        io.write(f"{self.path}/ts.traj", ts)  # optional
        with open(f"{self.path}/{name_ase_script}", "r") as infile:
            lines = infile.readlines()
        ase_script_lines = ["from ase.io import read\n",
                            "from ase.calculators.vasp import Vasp\n",
                            "from ase.vibrations import Vibrations\n",
                            "from datetime import datetime\n",
                            " \n",
                            "now = datetime.now()\n",
                            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
                            "with open('time.txt', 'w') as outfile:\n",
                            "\toutfile.write(f'{dt_string}: starting job ...\\n')\n",
                            " \n"]
        for i in range(len(lines)):
            if 'ase_calculator' in lines[i]:
                while lines[i] != ' \n':
                    ase_script_lines.append(lines[i])
                    i += 1
                break
        ase_script_lines += [" " + "\n",
                             f"ts = read('../ML-NEB.traj', index={np.argmax(self.energies)})\n",
                             "ts.set_calculator(ase_calculator)\n",
                             f"vib = Vibrations(ts, indices='-{self.num_atoms_adsorbate}:', name='vib')\n",
                             "vib.run()\n",
                             " \n",
                             "vib.summary(log='vib.txt')\n",
                             f"for mode in range({self.num_atoms_adsorbate}*3):\n",
                             "\tvib.write_mode(mode)\n",
                             " \n",
                             "now = datetime.now()\n",
                             "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
                             "with open('time.txt', 'a') as outfile:\n",
                             "\toutfile.write(f'{dt_string}: finished \\n')\n"]
        # Run vibrations
        init_dir = os.getcwd()
        os.mkdir(f"{self.path}/vibrations")
        with open(f"{self.path}/vibrations/{name_ase_script}", 'w') as outfile:
            for line in ase_script_lines:
                outfile.write(line)
        os.chdir(f"{self.path}/vibrations")
        os.system(f'cp ../{name_submission_script} .')
        self.submit(path=f"{self.path}/vibrations", name=f"vib_{self.name}")
        os.chdir(init_dir)

    def run_vibrations_native(self):
        """ Much faster than using ase.vibrations """
        #todo: add if for neb or dimer
        if len(self.energies) > 0:
            ts = io.read(f"{self.path}/ML-NEB.traj", index=np.argmax(self.energies))
            c = FixAtoms(indices=list(range(len(ts) - self.num_atoms_adsorbate)))
            ts.set_constraint(c)
            incar = Incar.from_file(path=self.path)
            incar.update_tag(key='IBRION', value='5')
            incar.update_tag(key='POTIM', value='0.03')
            incar.remove_tag(key='NPAR')
            kpoints = Kpoints.from_file(path=self.path)
            job = NewJobNative(
                path=f"{self.path}/vibrations",
                incar=incar,
                kpoints=kpoints,
                atoms=ts,
                potcars_dir=potcars_dir_cluster
            )
            job.create_job_dir()
            os.system(f"cp {path_submission_script_freq_native}/{name_submission_script} {self.path}/vibrations")
            self.submit(path=f"{self.path}/vibrations", name=f"vib_{self.name}")
        else:
            print(f"{self.name}: can't run vibrations, energy array is empty")

    def get_displacement_vector(self):
        line_start = 0
        i = 0
        with open(f"{self.path}/vibrations/OUTCAR", 'r') as infile:
            lines = infile.readlines()
        # Find start of imaginary vibration output
        for i in range(len(lines) - 1, 0, -1):
            if ' f/i= ' in lines[i]:
                line_start = i + 2
                break
        # Get length of imaginary vibration
        line = lines[i + 2]
        num_displacements = 0
        while line != ' \n':
            num_displacements += 1
            i += 1
            line = lines[i + 2]
        # Save vibration to numpy array
        displacement_vector = np.zeros((num_displacements, 3))
        for i in range(num_displacements):
            line = lines[line_start + i]
            displacement_vector[i][0] = line.split()[-3]
            displacement_vector[i][1] = line.split()[-2]
            displacement_vector[i][2] = line.split()[-1]
        return displacement_vector

    def run_dimer_vtst(self, ediffg=-0.02):
        ts = io.read(f"{self.path}/ML-NEB.traj", index=np.argmax(self.energies))
        incar = Incar.from_file(path=self.path)
        incar.update_tag(key='IBRION', value=3)
        incar.update_tag(key='POTIM', value=0)
        incar.add_tag(key='NSW', value=600)
        incar.add_tag(key='IOPT', value=2)
        incar.add_tag(key='ICHAIN', value=2)
        incar.add_tag(key='EDIFFG', value=ediffg)
        slab_pmg = AseAtomsAdaptor().get_structure(ts)
        job = NewJobNative(
            path=f"{self.path}/vibrations/dimer",
            incar=incar,
            kpoints=Kpoints(density=60, lat_x=slab_pmg.lattice.a, lat_y=slab_pmg.lattice.b),
            atoms=ts,
            potcars_dir=potcars_dir_cluster
        )
        job.create_job_dir()
        self.create_modecar()
        os.system(f"cp {path_submission_script_vtst_native}/{name_submission_script} {self.path}/vibrations/dimer")
        self.submit(path=f"{self.path}/vibrations/dimer", name=f"{self.name}")

    def create_modecar(self):
        displacement_vector = self.get_displacement_vector()
        with open(f"{self.path}/vibrations/dimer/MODECAR", 'w') as outfile:
            for i in range(len(displacement_vector)):
                outfile.write(f"{displacement_vector[i][0]} {displacement_vector[i][1]} {displacement_vector[i][2]}\n")



