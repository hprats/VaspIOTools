import os
import subprocess
from glob import glob
import json

from copy import deepcopy
import ase.neb
from ase import io

from vaspio.incar import Incar
from vaspio.kpoints import Kpoints

potcars_dir = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potpaw_PBE'
name_mlneb_script = 'ase_neb.py'
name_submission_script = 'vasp_native'
name_submission_script_ase = 'vasp_ase'
name_std_output = 'vasp.out'
path_qstat_list = '/home/ucechpr'

VASP_recommended_PP = {
    # Groups 1-2 & 13-18
    'H': 'H', 'He': 'He', 'Li': 'Li_sv', 'Be': 'Be', 'B': 'B', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F', 'Ne': 'Ne',
    'Na': 'Na_pv', 'Mg': 'Mg', 'Al': 'Al', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar',
    'K': 'K_sv', 'Ca': 'Ca_sv', 'Ga': 'Ga_d', 'Ge': 'Ge_d', 'As': 'As', 'Se': 'Se', 'Br': 'Br', 'Kr': 'Kr',
    'Rb': 'Rb_sv', 'Sr': 'Sr_sv', 'In': 'In_d', 'Sn': 'Sn_d', 'Sb': 'Sb', 'Te': 'Te', 'I': 'I', 'Xe': 'Xe',
    'Cs': 'Cs_sv', 'Ba': 'Ba_sv', 'Tl': 'Tl_d', 'Pb': 'Pb_d', 'Bi': 'Bi_d', 'Po': 'Po_d', 'At': 'At', 'Rn': 'Rn',
    # d block
    'Sc': 'Sc_sv', 'Ti': 'Ti_sv', 'V': 'V_sv', 'Cr': 'Cr_pv', 'Mn': 'Mn_pv', 'Fe': 'Fe', 'Co': 'Co', 'Ni': 'Ni', 'Cu': 'Cu', 'Zn': 'Zn',
    'Y': 'Y_sv', 'Zr': 'Zr_sv', 'Nb': 'Nb_sv', 'Mo': 'Mo_sv', 'Tc': 'Tc_pv', 'Ru': 'Ru_pv', 'Rh': 'Rh_pv', 'Pd': 'Pd', 'Ag': 'Ag', 'Cd': 'Cd',
    'La': 'La', 'Hf': 'Hf_pv', 'Ta': 'Ta_pv', 'W': 'W_sv', 'Re': 'Re', 'Os': 'Os', 'Ir': 'Ir', 'Pt': 'Pt', 'Au': 'Au', 'Hg': 'Hg',
    # f block
    'Ce': 'Ce'
}

MP_recommended_PP = deepcopy(VASP_recommended_PP)
MP_recommended_PP.update({'Be': 'Be_sv', 'Mg': 'Mg_pv', 'Ti': 'Ti_pv', 'Fe': 'Fe_pv', 'Ni': 'Ni_pv', 'Cu': 'Cu_pv',
                          'Nb': 'Nb_pv', 'Mo': 'Mo_pv', 'W': 'W_pv', 'Re': 'Re_pv'})

project_PP_dict = deepcopy(VASP_recommended_PP)  # for TMC4MPO project
project_PP_dict.update({'W': 'W_pv'})


class NewJobNative:
    """A class that represents a new VASP job.

    Attributes:
        path (str): The path of the job including the job name. Will be used as the name of the folder.
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
        >>> my_job = NewJobNative(
        >>>    path='/home/test/N2_gas',
        >>>    incar=Incar(tags),
        >>>    kpoints=Kpoints(numbers=[1, 1, 1]),
        >>>    atoms=ase.atoms.Atoms('N2', [(0, 0, 0), (0, 0, 1.12)],
        >>>    submission_script=[submission_script]
        >>> )
    """

    def __init__(self, path, incar, kpoints, atoms, submission_script):
        if isinstance(path, str):
            self.path = path
            self.name = path.split('/')[-1]
            self.incar = incar
            self.kpoints = kpoints
            self.atoms = atoms
            self.submission_script = submission_script

    def create_job_dir(self):
        """Creates a new directory, with the name job_name, and writes there the VASP input files
        i.e. INCAR, POSCAR, KPOINTS and POTCAR, and the submission script"""
        if not os.path.exists(self.path):
            os.mkdir(self.path)
            self.incar.write(self.path)
            self.kpoints.write(self.path)
            self.atoms.write(filename=f'{self.path}/POSCAR', vasp5=True)
            self.write_potcar()
            self.write_submission_script()
        else:
            print(f'{self.path} already exists (nothing done)')

    def write_submission_script(self):
        """Prints the submission script into a new file named vasp_sub"""
        with open(f"{self.path}/{name_submission_script}", 'w') as infile:
            for line in self.submission_script:
                infile.write(line + '\n')

    def write_potcar(self):
        """Writes the POTCAR file to the current directory
        First, a list of elements is created (elements_for_potcar), which determines which
        atomic POTCAR files will be concatenated to make the final POTCAR"""
        current_element = self.atoms.get_chemical_symbols()[0]
        elements_for_potcar = [current_element]
        for element in self.atoms.get_chemical_symbols()[1:]:
            if element != current_element:
                elements_for_potcar.append(element)
                current_element = element
        cmd = 'cat'
        for element in elements_for_potcar:
            cmd += f' {potcars_dir}/{project_PP_dict[element]}/POTCAR'
        os.system(f'{cmd} > {self.path}/POTCAR')


class JobNative:
    """A class that represents a VASP job, typically a finished one.

    Attributes:
        name (str): The name of the folder containing the job files.
        path (str): The full path to the job folder.
        energy (float): The total energy of the system.
        status (str): A string describing the status of the job (e.g. qw, r, fine, ...). See get_job_status().

    Examples:
        >>> job = JobNative(job_name='Ni_111#CO-top', path=f'path-to-folder/VC_001##CO2I_0-tM')
        >>> job.read_from_cluster()
        >>> print(job.energy)
    """

    def __init__(self, path, incar=None, energy=None, status=None, magmom=None):

        self.name = path.split('/')[-1]
        self.path = path
        self.energy = energy
        self.status = status
        self.magmom = magmom
        self.incar = incar

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
            energy = dct['energy']
            status = dct['status']
            magmom = dct['magmom']

            job = cls(path=path, incar=incar, energy=energy, status=status, magmom=magmom)
            return job

    def read_incar(self):
        """Read INCAR file and store its information as a Incar object"""
        with open(f'{self.path}/INCAR') as infile:
            lines = infile.readlines()
        incar_tags = {}
        for line in lines:
            tag = line.strip().split(' = ')[0]
            value = line.strip().split(' = ')[1]
            incar_tags[tag] = value
        return Incar(incar_tags)

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

    @classmethod
    def read_from_cluster(cls, path):
        """Don't need to read POSCAR and CONTCAR, since they won't be written to JSON"""
        job = cls(path=path)
        job.incar = job.read_incar()
        job.status = job.get_job_status()
        if job.status == 'fine':
            job.energy = job.get_energy_oszicar()
            job.magmom = job.get_magmom()
        return job

    def converged(self):
        if 'NSW' not in self.incar.tags:  # is SPE
            return ' 1 F= ' in \
                   str(subprocess.check_output(f"tail -n4 {self.path}/{name_std_output}", shell=True))
        else:
            return 'reached required accuracy' in \
               str(subprocess.check_output(f"tail -n4 {self.path}/{name_std_output}", shell=True))

    def bracketing_error(self):
        return 'bracketing' in \
               str(subprocess.check_output(f"tail -n7 {self.path}/{name_std_output}", shell=True))

    def get_dE_last_two_steps(self):
        output = str(subprocess.check_output(f"grep F {self.path}/OSZICAR | tail -n2", shell=True))
        last = float(output.split('\\n')[1].split('  d E')[0].split(' ')[-1])
        previous = float(output.split('\\n')[0].split('  d E')[0].split(' ')[-1])
        diff = abs(last - previous)
        return diff

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

    def core(self):
        if len(glob(f'{self.path}/core.*')) > 0:
            return True
        else:
            return False

    def bad_termination(self):
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_std_output}", shell=True))
        if 'BAD TERMINATION' in output:
            return True
        else:
            return False

    def empty_poscar(self):
        file = open(f"{self.path}/POSCAR")
        lines = [line for line in file]
        if len(lines) == 0:
            return True
        else:
            return False

    def empty_contcar(self):
        file = open(f"{self.path}/CONTCAR")
        lines = [line for line in file]
        if len(lines) == 0:
            return True
        else:
            return False

    def other_error(self):
        output = str(subprocess.check_output(f"tail -n10 {self.path}/{name_std_output}", shell=True))
        if 'error' in output and 'errors must be expected' not in output:
            return True
        else:
            return False

    def vasp_bin_not_loaded(self):  # if this doesnt work, check other version in jobwp2.py
        try:
            std_error_file = glob(f"{self.path}/*.e*")[0]
            return 'cannot be loaded' in str(subprocess.check_output(f"tail {std_error_file}", shell=True))
        except subprocess.CalledProcessError:
            return False

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
        """If modified here, modify also in vasp_tools.py."""
        with open(f'{path_qstat_list}/qstat_list.txt') as infile:
            qstat_list = infile.readlines()
        in_queue, status = self.check_queue(qstat_list)
        if in_queue:
            job_status = status
        elif os.path.isfile(f"{self.path}/README"):
            with open(f"{self.path}/README") as f:
                readme_info = f.readlines()[0].strip()
            job_status = f'README: {readme_info}'
        elif self.empty_poscar():
            job_status = 'empty poscar'
        elif not os.path.isfile(f"{self.path}/{name_std_output}"):
            job_status = 'not submitted'
        elif self.converged():
            job_status = 'fine'
        elif not os.path.isfile(f"{self.path}/OSZICAR"):
            job_status = 'error'
        else:
            if self.vasp_bin_not_loaded():
                job_status = 'VASP bin not loaded'
            elif self.bracketing_error():
                if self.get_dE_last_two_steps() <= 0.01:
                    job_status = 'fine'
                else:
                    job_status = 'bracketing'
            elif self.core():
                job_status = 'core files generated'
            elif self.bad_termination():
                if self.get_dE_last_two_steps() <= 0.01:
                    job_status = 'fine'
                else:
                    job_status = 'bad termination'
            elif self.nsw_reached():
                job_status = 'NSW reached'
            elif self.nelm_reached():
                job_status = 'max wallclock, NELM reached'
            elif self.other_error():
                job_status = 'error'
            else:
                job_status = 'max wallclock'
        return job_status

    def rm_vasp_outputs(self):
        files_list = [f for f in os.listdir(self.path) if os.path.isfile(os.path.join(self.path, f))]
        for file in files_list:
            if file not in ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'vasp_native', 'submitted']:
                os.remove(f"{self.path}/{file}")

    def restart(self):
        init_dir = os.getcwd()  # does not conflict with external variables with same name
        self.rm_vasp_outputs()
        os.chdir(self.path)
        os.system(f'qsub -N {self.name} vasp_native')
        os.chdir(init_dir)

    def refine(self, dict_new_tags):
        init_dir = os.getcwd()
        os.chdir(self.path)
        num_previous_refines = str(len(glob('ref*/')))
        if self.empty_contcar():
            print(f'CHECK: Cannot refine {self.name}: empty CONTCAR')
        else:
            os.system(f'mkdir ref{num_previous_refines}; cp * ref{num_previous_refines}; cp CONTCAR POSCAR')
            self.rm_vasp_outputs()
            for tag in dict_new_tags:
                self.incar.update_tag(key=tag, value=dict_new_tags[tag])
                self.incar.write(self.path)
            os.system(f'qsub -N {self.name} vasp_native')
        os.chdir(init_dir)