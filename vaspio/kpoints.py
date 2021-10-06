import os
import numpy as np


class Kpoints:
    """A class that represents a KPOINTS file.

    Attributes:
        is_slab (bool): True for slab systems, False for bulk or gas.
            If is_slab = True, the number of kpoints in the z direction is one.
        density (int): The density of kpoints along the x and y direction (is_slab = True)
            or x, y and z (is_slab = False), in kpoints / Å.
        numbers (int): A list of integers representing the number of kpoints along the
            x, y and z directions.
        lat_x (float): length of the supercell in the x direction, in Å.
        lat_y (float): length of the supercell in the y direction, in Å.
        lat_z (float): length of the supercell in the z direction, in Å.

    Examples:

        Slab calculation with density:
        >>> my_kpoints = Kpoints(density=60, lat_x=6.5, lat_y=5.7)

        Slab calculation with number of kpoints:
        >>> my_kpoints = Kpoints(numbers=[6, 6, 1])

        Bulk calculation with density:
        >>> my_kpoints = Kpoints(density=60, lat_x=3.5, lat_y=3.5, lat_z=2.9)

        Isolated molecule calculation:
        >>> my_kpoints = Kpoints(numbers=[1, 1, 1])
    """

    def __init__(self, is_slab=True, density=None, numbers=None,
                 lat_x=None, lat_y=None, lat_z=None):
        self._is_slab = is_slab
        self._density = density
        if numbers is None:   # Density specified
            if isinstance(density, int):
                self.num_x = int(np.ceil(density / lat_x))
                self.num_y = int(np.ceil(density / lat_y))
                if is_slab:
                    self.num_z = 1
                else:
                    self.num_z = int(np.ceil(density / lat_z))
            else:
                print("K-points density must be int")
        else:  # Numbers specified
            if isinstance(numbers, list):
                self.num_x = numbers[0]
                self.num_y = numbers[1]
                self.num_z = numbers[2]
            else:
                print("K-point numbers must be a list")

    @classmethod
    def from_dict(cls, dct):
        is_slab = dct['_is_slab']
        density = dct['_density']
        numbers = dct['numbers']
        lat_x = dct['lat_x']
        lat_y = dct['lat_y']
        lat_z = dct['lat_z']

        kpoints = cls(is_slab=is_slab, density=density, numbers=numbers,
                      lat_x=lat_x, lat_y=lat_y, lat_z=lat_z)
        return kpoints

    @property
    def is_bulk(self):
        return self._is_slab

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, new_density):
        if isinstance(new_density, int):
            self._density = new_density
        else:
            print("New density for KPOINTS file is not an integer")

    def write(self, path='.'):
        """Write the KPOINTS file."""
        initial_directory = os.getcwd()
        os.chdir(path)
        f = open('KPOINTS', 'w')
        f.write('Automatic mesh\n')
        f.write('0\n')
        f.write('Gamma\n')
        f.write(f"{self.num_x} {self.num_y} {self.num_z}\n")
        f.write('0 0 0\n')
        f.close()
        os.chdir(initial_directory)

    @property
    def is_slab(self):
        return self._is_slab
