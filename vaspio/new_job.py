import os

from ase.constraints import FixAtoms
from ase.neb import NEB

from vaspio.write_functions import write_potcar, write_modecar


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


class NewVibrationAnalysisNative:
    """A class that represents a new VASP job."""
    def __init__(self, incar, kpoints, atoms, pp_dict, pp_path, num_free_atoms, potim=0.03):

        incar.update_tag(key='IBRION', value='5')
        incar.update_tag(key='POTIM', value=potim)
        incar.remove_tag(key='NPAR')

        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.pp_dict = pp_dict
        self.pp_path = pp_path
        self.num_free_atoms = num_free_atoms

        c = FixAtoms(indices=list(range(len(self.atoms) - num_free_atoms)))
        self.atoms.set_constraint(c)

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


class NewDimerVTST:
    """A class that represents a new VASP job."""
    def __init__(self, incar, kpoints, atoms, pp_dict, pp_path, displacement_vector):

        incar.add_tag(key='ICHAIN', value=2)
        incar.update_tag(key='IBRION', value=3)
        incar.update_tag(key='POTIM', value=0)
        incar.add_tag(key='IOPT', value=2)

        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.pp_dict = pp_dict
        self.pp_path = pp_path
        self.displacement_vector = displacement_vector

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
            write_modecar(job_path, self.displacement_vector)
        else:
            print(f'{job_path} already exists (nothing done)')


class NewNebNative:
    """A class that represents a new VASP job."""
    def __init__(self, incar, kpoints, atoms, pp_dict, pp_path, images, atoms_initial, atoms_final, energy_initial,
                 energy_final):
        self.incar = incar
        self.kpoints = kpoints
        self.atoms = atoms
        self.pp_dict = pp_dict
        self.pp_path = pp_path
        self.images = images
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final
        self.energy_initial = energy_initial
        self.energy_final = energy_final

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
            # Write initial and final energies
            with open(f"{job_path}/energies.txt", 'w') as outfile:
                outfile.write(f'{self.energy_initial}\n')
                outfile.write(f'{self.energy_final}\n')
            # Create all folders
            for i in range(self.images + 2):
                if i <= 9:
                    os.mkdir(f'{job_path}/0{i}')
                else:
                    os.mkdir(f'{job_path}/{i}')
            # Write POSCAR files for reactants and products
            self.atoms_initial.write(filename=f'{job_path}/00/POSCAR', vasp5=True)
            if self.images <= 9:
                self.atoms_final.write(filename=f'{job_path}/0{self.images + 1}/POSCAR', vasp5=True)
            else:
                self.atoms_final.write(filename=f'{job_path}/{self.images + 1}/POSCAR', vasp5=True)
            # Write POSCAR files for images
            constraints = self.atoms_initial.constraints
            list_images = [self.atoms_initial]
            for i in range(self.images):
                image = self.atoms_initial.copy()
                image.set_constraint(constraints)
                list_images.append(image)
            list_images.append(self.atoms_final)
            neb = NEB(list_images, climb=False, k=0.5)
            neb.interpolate('idpp')
            for i in range(1, self.images + 1):
                if i <= 9:
                    list_images[i].write(filename=f'{job_path}/0{i}/POSCAR', vasp5=True)
                else:
                    list_images[i].write(filename=f'{job_path}/{i}/POSCAR', vasp5=True)
        else:
            print(f'{job_path} already exists (nothing done)')
