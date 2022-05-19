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
        self.is_slab = is_slab
        self.density = density
        if numbers is None:   # Density specified
            if isinstance(density, int):
                self.num_x = int(np.ceil(density / lat_x))
                self.num_y = int(np.ceil(density / lat_y))
                if is_slab:
                    self.num_z = 1
                else:
                    self.num_z = int(np.ceil(density / lat_z))
            else:
                print(f"K-points density must be int, is {type(self.density)}")
        else:  # Numbers specified
            if isinstance(numbers, list):
                self.num_x = numbers[0]
                self.num_y = numbers[1]
                self.num_z = numbers[2]
            else:
                print("K-point numbers must be a list")

    def write(self, path):
        """Write the KPOINTS file."""
        with open(f"{path}/KPOINTS", 'w') as infile:
            infile.write('Automatic mesh\n')
            infile.write('0\n')
            infile.write('Gamma\n')
            infile.write(f"{self.num_x} {self.num_y} {self.num_z}\n")
            infile.write('0 0 0\n')

    @classmethod
    def from_dict(cls, dct):
        is_slab = dct['is_slab']
        if 'density' in dct:
            density = dct['density']
        else:
            density = None
        if 'num_x' in dct and 'num_y' in dct and 'num_z' in dct:
            numbers = [int(dct['num_x']), int(dct['num_y']), int(dct['num_z'])]
        else:
            numbers = None
        kpoints = cls(is_slab=is_slab, density=density, numbers=numbers)
        return kpoints

    @classmethod
    def from_file(cls, path):
        """Read an existing KPOINTS file and store it as a KPOINTS object"""
        with open(f"{path}/KPOINTS", 'r') as infile:
            lines = infile.readlines()
        num_x = lines[-2].split(' ')[0].strip()
        num_y = lines[-2].split(' ')[1].strip()
        num_z = lines[-2].split(' ')[2].strip()
        kpoints = cls(numbers=[num_x, num_y, num_z])
        return kpoints
