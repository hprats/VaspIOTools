import os
import numpy as np
import subprocess
from glob import glob
import json

import ase.neb
from ase import io
from ase.constraints import FixAtoms

from vaspio.variables import *
from vaspio.incar import Incar
from vaspio.kpoints import Kpoints
from vaspio.jobs import NewJobNative


class NewNebNative:

    def __init__(self, path, images, incar, kpoints, atoms_initial, atoms_final, submission_script,
                 energy_initial, energy_final):
        if isinstance(path, str):
            self.path = path
            self.images = images
            self.name = path.split('/')[-1]
            self.incar = incar
            self.kpoints = kpoints
            self.atoms_initial = atoms_initial
            self.atoms_final = atoms_final
            self.submission_script = submission_script
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
            self.write_submission_script()
            self.write_energies()
            # Create all folders
            for i in range(self.images+2):
                if i <= 9:
                    os.mkdir(f'{self.path}/0{i}')
                else:
                    os.mkdir(f'{self.path}/{i}')
            # Write POSCAR files for reactants and products
            self.atoms_initial.write(filename=f'{self.path}/00/POSCAR', vasp5=True)
            if self.images <= 9:
                self.atoms_final.write(filename=f'{self.path}/0{self.images+1}/POSCAR', vasp5=True)
            else:
                self.atoms_final.write(filename=f'{self.path}/{self.images+1}/POSCAR', vasp5=True)
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
            for i in range(1, self.images+1):
                if i <= 9:
                    list_images[i].write(filename=f'{self.path}/0{i}/POSCAR', vasp5=True)
                else:
                    list_images[i].write(filename=f'{self.path}/{i}/POSCAR', vasp5=True)
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_submission_script(self):
        """Prints the submission script into a new file named vasp_sub"""
        with open(f"{self.path}/{name_submission_script}", 'w') as infile:
            for line in self.submission_script:
                infile.write(line + '\n')

    def write_energies(self):
        """Prints the submission script into a new file named vasp_sub"""
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
            cmd += f' {potcars_dir}/{project_PP_dict[element]}/POTCAR'
        os.system(f'{cmd} > {self.path}/POTCAR')


class NewNebML:

    def __init__(self, path, n_images, incar_tags, n_kpoints, atoms_initial, atoms_final, submission_script):
        if isinstance(path, str):
            self.path = path
            self.n_images = n_images
            self.name = path.split('/')[-1]
            self.incar_tags = incar_tags
            self.n_kpoints = n_kpoints
            self.atoms_initial = atoms_initial
            self.atoms_final = atoms_final
            self.submission_script = submission_script

    def create_job_dir(self):
        """Creates a new directory, with the name job_name, and writes there the VASP input files
        i.e. INCAR, POSCAR, KPOINTS and POTCAR, and the submission script"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.write_submission_script()
            self.write_ase_script()
            self.write_trajectory_files()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_submission_script(self):
        """Prints the submission script into a new file named vasp_sub"""
        with open(f"{self.path}/{name_submission_script}", 'w') as outfile:
            for line in self.submission_script[:-1]:
                outfile.write(line + '\n')
            outfile.write('module load python3/3.9\n')
            outfile.write(f'python {name_ase_script} > vasp.out\n')

    def write_ase_script(self):
        ase_script = [
            "from ase.io import read",
            "from ase.optimize import BFGS",
            "from ase.calculators.vasp import Vasp",
            "import shutil",
            "import copy",
            "from catlearn.optimize.mlneb import MLNEB",
            "from datetime import datetime",
            " ",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'w') as outfile:",
            "\toutfile.write(f'{dt_string}: starting job ...\\n')",
            " ",
            "ase_calculator = Vasp(setups={'base': 'recommended', 'W': '_pv'},",
            f"\tediff={self.incar_tags['EDIFF']},",
            f"\tnelm={self.incar_tags['NELM']},",
            f"\tismear=1,",
            f"\tsigma=0.2,",
            "\tlwave=False,",
            "\tlcharg=False,",
            f"\tencut=415,",
            f"\talgo='Fast',",
            f"\tlreal='Auto',",
            f"\txc='pbe',",
            f"\tgga='PE',",
            f"\tivdw=11,",
            f"\tldipol=True,",
            f"\tidipol=3,",
            f"\tdipol=[0.5, 0.5, 0.5],",
            f"\tlasph=True,",
            f"\tispin={self.incar_tags['ISPIN']},",
            f"\tnpar={self.incar_tags['NPAR']},",
            "\t# Kpoints",
            f"\tkpts=({self.n_kpoints[0]}, {self.n_kpoints[1]}, {self.n_kpoints[2]}),",
            f"\tgamma=True,",
            "\t)",
            " ",
            "# Optimize initial state:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting optimization of initial state ...\\n')",
            "slab = read('./optimized_structures/initial.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='initial.traj')",
            f"qn.run(fmax={abs(self.incar_tags['EDIFFG'])})",
            "shutil.copy('./initial.traj', './optimized_structures/initial.traj')",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: initial state optimized; starting optimization of final state ...\\n')",
            " ",
            "# Optimize final state:",
            "slab = read('./optimized_structures/final.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='final.traj')",
            f"qn.run(fmax={abs(self.incar_tags['EDIFFG'])})",
            "shutil.copy('./final.traj', './optimized_structures/final.traj')",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: final state optimized, starting ML-NEB ...\\n')",
            " ",
            "neb_catlearn = MLNEB(start='initial.traj', end='final.traj',",
            "\tase_calc=copy.deepcopy(ase_calculator),",
            f"\tn_images={self.n_images},",
            "\tinterpolation='idpp',",
            "\t#restart=True,",
            "\t)",
            " ",
            f"neb_catlearn.run(fmax=0.05, trajectory='ML-NEB.traj')",
            " ",
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

    def write_trajectory_files(self):
        os.mkdir(f"{self.path}/optimized_structures")
        io.write(f"{self.path}/optimized_structures/initial.traj", self.atoms_initial)
        io.write(f"{self.path}/optimized_structures/final.traj", self.atoms_final)
        list_images = [self.atoms_initial]
        constraints = self.atoms_initial.constraints
        for i in range(self.n_images-2):
            image = self.atoms_initial.copy()
            image.set_constraint(constraints)
            list_images.append(image)
        list_images.append(self.atoms_final)
        neb = ase.neb.NEB(list_images, climb=False, k=0.5)
        neb.interpolate('idpp')
        io.write(f"{self.path}/optimized_structures/list_images.traj", list_images)


class NebML:

    def __init__(self, path, status=None, energies=None, vibrations=None, num_atoms_adsorbate=None):
        self.path = path
        self.name = path.split('/')[-1]
        self.status = status
        self.energies = energies
        self.vibrations = vibrations
        self.num_atoms_adsorbate = num_atoms_adsorbate

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
        images = io.read(f"{self.path}/ML-NEB.traj", index=":")
        energies = np.zeros(len(images))
        for i in range(len(images)):
            energies[i] = images[i].get_potential_energy()
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
        job.status = job.get_job_status()
        if job.status == 'done':
            job.energies = job.get_energies()
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

    def new_minimum(self):
        new_minimum = False
        for i in range(1, np.argmax(self.energies)):
            if self.energies[i] < self.energies[0]:
                new_minimum = True
        for i in range(np.argmax(self.energies)+1, -1):
            if self.energies[i] < self.energies[-1]:
                new_minimum = True
        return new_minimum

    def check_queue(self, qstat_list):
        in_queue = False
        status = None
        for line in qstat_list:
            job_name = line.split()[2]
            if job_name == self.name:
                in_queue = True
                status = line.split()[4]
                break
        return in_queue, status

    def get_job_status(self):
        """Execute it on the cluster."""
        with open(f'{path_qstat_list}/qstat_list.txt') as infile:
            qstat_list = infile.readlines()
        in_queue, status = self.check_queue(qstat_list)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif os.path.isdir(f"{self.path}/vibrations"):
            if self.freq_done():
                if self.is_ts():
                    job_status = 'done'
                else:
                    job_status = 'not TS'
            else:
                job_status = 'vibrations not converged'
        elif self.NEB_converged():
            if self.new_minimum():
                job_status = 'new minimum'
            else:
                job_status = 'NEB converged'
        else:
            job_status = 'max wallclock'
        return job_status

    def submit(self):
        init_dir = os.getcwd()
        os.chdir(self.path)
        os.system(f'qsub -N {self.name} {name_submission_script}')
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
        init_dir = os.getcwd()
        os.chdir(self.path)
        os.system(f'qsub -N {self.name} {name_submission_script}')
        os.chdir(init_dir)

    def run_vibrations_ase(self):
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
        os.system(f'qsub -N vib_{self.name} {name_submission_script}')
        os.chdir(init_dir)

    def run_vibrations_native(self, submission_script):
        """ Much faster than using ase.vibrations """
        ts = io.read(f"{self.path}/ML-NEB.traj", index=np.argmax(self.energies))
        c = FixAtoms(indices=list(range(len(ts)-self.num_atoms_adsorbate)))
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
            submission_script=submission_script
        )
        job.create_job_dir()

