import os

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
