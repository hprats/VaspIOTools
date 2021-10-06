import os
import ase


def get_elements(list_elements_repeated):
    elements = []
    current_element = list_elements_repeated[0]
    elements.append(current_element)
    for element in list_elements_repeated[1:]:
        if element != current_element:
            elements.append(element)
            current_element = element
    return elements


class Poscar:
    """A class that represents a POSCAR file.

    Attributes:
        atoms (ase.atoms.Atoms object): an Atoms object representing the system.

    Examples:
    >>> atoms = ase.atoms.Atoms('N2', [(0, 0, 0), (0, 0, 1.1)])
    >>> my_poscar = Poscar(atoms)
    """

    def __init__(self, atoms=None):
        if isinstance(atoms, ase.atoms.Atoms):
            self.atoms = atoms
        else:
            print("atoms passed is not a valid ase.atoms.Atoms object")

        self.elements = get_elements(atoms.get_chemical_symbols())

    @classmethod
    def from_dict(cls, dct):
        atoms = dct['atoms']

        incar = cls(atoms=atoms)
        return incar


    def write(self, path='.'):
        """Write the POSCAR file."""
        initial_directory = os.getcwd()
        os.chdir(path)
        self.atoms.write(filename='POSCAR', vasp5=True)
        os.chdir(initial_directory)
