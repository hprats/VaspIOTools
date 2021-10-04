import os


class Potcar:
    """A class that represents a POTCAR file.

    Attributes:
        elements (list): List of elements needed to create the POTCAR.
            It must be consistent with the POSCAR file, and can be created as follows:
                elements = Poscar(atoms).elements

    Examples:
    >>> my_potcar = Potcar(elements)
    """

    POTCARS_DIR = '/Users/hectorpratsgarcia/PycharmProjects/tmc4mpo/potcars'

    def __init__(self, elements=None):
        if elements is None:
            self._elements = []
        else:
            self._elements = elements

    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, new_elements):
        if isinstance(new_elements, list):
            self._elements = new_elements
        else:
            print("New elements variable for POTCAR file is not a list")

    def write(self, path='.'):
        """Write the POTCAR file."""
        initial_directory = os.getcwd()
        os.chdir(path)
        with open('POTCAR', 'w') as outfile:
            for element in self._elements:
                with open(f"{Potcar.POTCARS_DIR}/POTCAR-{element}") as infile:
                    for line in infile:
                        outfile.write(line)
        os.chdir(initial_directory)
