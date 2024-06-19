import os

from vaspio.input_files.potcar import write_potcar


class NewJobNative:
    """A class that represents a new VASP job."""
    def __init__(self, incar, kpoints, atoms, pp_dict, pp_path):
        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.pp_dict = pp_dict
        self.pp_path = pp_path

    def create_job_dir(self, job_path):
        """Creates a new director and writes there the VASP input files i.e. INCAR, POSCAR, KPOINTS and POTCAR"""
        if not os.path.exists(job_path):
            os.mkdir(job_path)
            self.incar.write(job_path)
            self.kpoints.write(job_path)
            self.atoms.write(filename=f'{job_path}/POSCAR', vasp5=True)
            write_potcar(job_path=job_path,
                         poscar_elements=self.atoms.get_chemical_symbols(),
                         pp_dict=self.pp_dict,
                         pp_path=self.pp_path)
        else:
            print(f'{job_path} already exists (nothing done)')
