import os
import subprocess
from glob import glob

from ase.io import read as ase_read

from vaspio.input_files.incar import Incar
from vaspio.input_files.poscar import Poscar
from vaspio.input_files.kpoints import Kpoints
from vaspio.input_files.potcar import Potcar


class NewJob:

    def __init__(self, job_name=None, incar=None, poscar=None, kpoints=None, potcar=None, submission_script=None):
        if isinstance(job_name, str):
            self.job_name = job_name
        else:
            print("job_name must be a string")
        if isinstance(incar, Incar):
            self._incar = incar
        else:
            print("incar must be an Incar object")
        if isinstance(poscar, Poscar):
            self._poscar = poscar
        else:
            print("poscar must be a Poscar object")
        if isinstance(kpoints, Kpoints):
            self._kpoints = kpoints
        else:
            print("kpoints must be a Kpoints object")
        if isinstance(potcar, Potcar):
            self._potcar = potcar
        else:
            print("potcar must be a Potcar object")
        if isinstance(submission_script, list):
            self._submission_script = submission_script
        else:
            print("submission_script must be a list of strings")

    @property
    def incar(self):
        return self._incar

    @incar.setter
    def incar(self, new_incar):
        if isinstance(new_incar, Incar):
            self._incar = new_incar
        else:
            print("new_incar for NewJob is not an Incar object")

    @property
    def poscar(self):
        return self._poscar

    @poscar.setter
    def poscar(self, new_poscar):
        if isinstance(new_poscar, Poscar):
            self._poscar = new_poscar
        else:
            print("new_poscar for NewJob is not an Poscar object")

    @property
    def kpoints(self):
        return self._kpoints

    @kpoints.setter
    def kpoints(self, new_kpoints):
        if isinstance(new_kpoints, Kpoints):
            self._kpoints = new_kpoints
        else:
            print("new_kpoints for NewJob is not an Kpoints object")

    @property
    def potcar(self):
        return self._potcar

    @potcar.setter
    def potcar(self, new_potcar):
        if isinstance(new_potcar, Potcar):
            self._potcar = new_potcar
        else:
            print("new_potcar for NewJob is not an Potcar object")

    def write_submission_script(self):
        f = open('vasp_sub', 'w')
        for line in self._submission_script:
            f.write(line + '\n')
        f.close()

    def create_job_dir(self, path='.'):
        if self._incar is None:
            print("incar is missing")
        elif self._poscar is None:
            print("poscar is missing")
        elif self._kpoints is None:
            print("kpoints is missing")
        elif self._potcar is None:
            print("potcar is missing")
        else:
            initial_directory = os.getcwd()
            os.chdir(path)
            os.mkdir(self.job_name)
            os.chdir(self.job_name)
            self._incar.write()
            self._poscar.write()
            self._kpoints.write()
            self._potcar.write()
            self.write_submission_script()
            os.chdir(initial_directory)


class Job:

    def __init__(self, job_name, path=None,
                 incar=None, poscar=None, kpoints=None, potcar=None, contcar=None,
                 energy=None, num_restarts=None, status=None, name_std_output='vasp.out'):
        if not isinstance(job_name, str):
            print("job_name must be a string")
        if not isinstance(path, str):
            print("path must be a string")

        self.job_name = job_name
        self.path = path
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints
        self.potcar = potcar
        self.contcar = contcar
        self.energy = energy
        self.num_restarts = num_restarts
        self.status = status
        self.name_std_output = name_std_output

    def read(self, deep_read=False):
        os.chdir(self.path)
        if os.path.isfile('INCAR'):
            self.incar = self.read_incar()
        if os.path.isfile('POSCAR') and deep_read:
            self.poscar = self.read_poscar()
        if os.path.isfile('KPOINTS') and deep_read:
            self.kpoints = self.read_kpoints()
        if os.path.isfile('POTCAR') and deep_read:
            self.potcar = self.read_potcar()
        if os.path.isfile('CONTCAR') and deep_read:
            self.contcar = self.read_contcar()
        self.status = self.get_job_status()
        if self.status == 'fine':
            self.energy = self.get_energy_oszicar()

    def read_incar(self):
        with open(f'{self.path}/INCAR') as infile:
            lines = infile.readlines()
        incar_tags = {}
        for line in lines:
            tag = line.strip().split(' = ')[0]
            value = line.strip().split(' = ')[1]
            incar_tags[tag] = value
        return Incar(incar_tags)

    def read_poscar(self):
        atoms = ase_read(f'{self.path}/POSCAR')
        return Poscar(atoms)

    def read_contcar(self):  # TODO: think if I want Job.contcar to be a Poscar object or not
        atoms = ase_read(f'{self.path}/CONTCAR')
        return Poscar(atoms)

    def read_kpoints(self):
        pass

    def read_potcar(self):
        pass

    def converged(self):
        return 'reached required accuracy' in str(subprocess.check_output(f"tail -n4 {self.name_std_output}", shell=True))

    def bracketing_error(self):
        return 'fatal error in bracketing' in str(subprocess.check_output(f"tail -n7 {self.name_std_output}", shell=True))

    def get_bracketing_diff(self):
        output = str(subprocess.check_output('grep F OSZICAR | tail -n2', shell=True))
        last = float(output.split('\\n')[1].split('  d E')[0].split(' ')[-1])
        previous = float(output.split('\\n')[0].split('  d E')[0].split(' ')[-1])
        diff = abs(last - previous)
        return diff

    def nsw_reached(self):
        try:
            output = str(subprocess.check_output("grep F OSZICAR", shell=True))
            return f"{self.incar.tags['NSW']} F=" in output
        except:
            return False

    def nelm_reached(self):
        nelm = self.incar.tags['NELM']
        return f"RMM: {nelm}" in str(subprocess.check_output(f"tail -n{int(nelm) + 15} OSZICAR", shell=True))

    def other_error(self):
        return 'error' in str(subprocess.check_output(f"tail -n10 {self.name_std_output}", shell=True))

    def get_job_status(self):
        if not os.path.isfile(self.name_std_output):
            job_status = 'qw'
        elif os.stat(self.name_std_output).st_size == 0:  # standard output is empty
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
        with open('OSZICAR') as infile:
            lines = infile.readlines()
        final = len(lines) - 1
        for i in range(final, 0, -1):
            if ' F= ' in lines[i]:
                energy = float(lines[i].split()[4])
                break
        return energy
