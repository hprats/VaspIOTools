import os
import subprocess

from ase.io import read as ase_read

from vaspio.incar import Incar


class NewJob:

    POTCARS_DIR = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potcars'

    def __init__(self, job_name, incar, kpoints, atoms, submission_script):
        if isinstance(job_name, str):
            self.job_name = job_name
            self.incar = incar
            self.kpoints = kpoints
            self.atoms = atoms
            self.submission_script = submission_script

    def write_submission_script(self):
        f = open('vasp_sub', 'w')
        for line in self.submission_script:
            f.write(line + '\n')
        f.close()

    def create_job_dir(self, path='.'):
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

    POTCARS_DIR = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potcars'

    def __init__(self, job_name, path=None,
                 incar=None, kpoints=None, atoms_poscar=None, atoms_contcar=None,
                 energy=None, num_restarts=None, status=None, name_std_output='vasp.out'):
        if not isinstance(job_name, str):
            print("job_name must be a string")
        if not isinstance(path, str):
            print("path must be a string")

        self.job_name = job_name
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
        os.chdir(self.path)
        if os.path.isfile('INCAR'):
            self.incar = self.read_incar()
        if os.path.isfile('KPOINTS'):
            self.kpoints = self.read_kpoints()
        if os.path.isfile('POSCAR'):
            self.atoms_poscar = ase_read('POSCAR')
        if os.path.isfile('CONTCAR'):
            self.atoms_contcar = ase_read('CONTCAR')
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

    def read_kpoints(self):
        pass

    def converged(self):
        return 'reached required accuracy' in str(subprocess.check_output(f"tail -n4 {self.name_std_output}", shell=True))

    def bracketing_error(self):
        return 'fatal error in bracketing' in str(subprocess.check_output(f"tail -n7 {self.name_std_output}", shell=True))

    @staticmethod
    def get_bracketing_diff():
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

    @staticmethod
    def get_energy_oszicar():
        energy = None
        with open('OSZICAR') as infile:
            lines = infile.readlines()
        final = len(lines) - 1
        for i in range(final, 0, -1):
            if ' F= ' in lines[i]:
                energy = float(lines[i].split()[4])
                break
        return energy
