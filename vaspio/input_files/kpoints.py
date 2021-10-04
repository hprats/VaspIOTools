import os
import numpy as np


class Kpoints:

    def __init__(self, is_bulk=False, density=None, numbers=None,
                 lat_x=None, lat_y=None, lat_z=None):
        self._is_bulk = is_bulk
        self._density = density
        if numbers is None:   # Density specified
            if isinstance(density, int):
                self.num_x = int(np.ceil(density / lat_x))
                self.num_y = int(np.ceil(density / lat_y))
                if is_bulk:
                    self.num_z = int(np.ceil(density / lat_z))
                else:
                    self.num_z = 1
            else:
                print("K-points density must be int")
        else:  # Numbers specified
            if isinstance(numbers, list):
                self.num_x = numbers[0]
                self.num_y = numbers[1]
                self.num_z = numbers[2]
            else:
                print("K-point numbers must be a list")

    @property
    def is_bulk(self):
        return self._is_bulk

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
        initial_directory = os.getcwd()
        os.chdir(path)
        f = open('KPOINTS', 'w')
        f.write('Automatic mesh\n')
        f.write('0\n')
        f.write('Gamma\n')
        f.write(f"{self.num_x} {self.num_y} {self.num_z}")
        f.write('0 0 0\n')
        f.close()
        os.chdir(initial_directory)
