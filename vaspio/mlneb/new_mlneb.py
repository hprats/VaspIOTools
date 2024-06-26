import os
import ase.mep
from ase import io


class NewNebML:   # todo: update to select specific PP_dict in run.py

    def __init__(self, incar, kpoints, n_images, fmax, atoms_initial, atoms_final, pp_dict, pp_path):
        incar.add_tag(key='n_images', value=n_images)
        incar.add_tag(key='fmax', value=fmax)
        self.incar = incar   # todo: replace by incar_tags?
        self.kpoints = kpoints
        self.n_images = n_images
        self.fmax = fmax
        self.atoms_initial = atoms_initial
        self.atoms_final = atoms_final
        self.pp_dict = pp_dict   # todo: currently this is ignored
        self.pp_path = pp_path

    def create_job_dir(self, job_path):
        if not os.path.exists(job_path):
            os.mkdir(job_path)
            self.incar.write(job_path)
            self.kpoints.write(job_path)
            self.write_trajectory_files(job_path)
            self.write_ase_script(job_path)
        else:
            print(f'{job_path} already exists (nothing done)')

    def write_trajectory_files(self, job_path):
        os.mkdir(f"{job_path}/optimized_structures")
        io.write(f"{job_path}/optimized_structures/initial.traj", self.atoms_initial)
        io.write(f"{job_path}/optimized_structures/final.traj", self.atoms_final)
        list_images = [self.atoms_initial]
        constraints = self.atoms_initial.constraints
        for i in range(self.n_images - 2):
            image = self.atoms_initial.copy()
            image.set_constraint(constraints)
            list_images.append(image)
        list_images.append(self.atoms_final)
        neb = ase.mep.NEB(list_images, climb=False, k=0.5)
        neb.interpolate('idpp')
        io.write(f"{job_path}/optimized_structures/list_images.traj", list_images)

    def write_ase_script(self, job_path):
        print("Warning, currently the pp_dict and pp_path flags are ignored and setups={'base': 'recommended'} is used")
        ase_script = [
            # Imports
            "from ase.io import read",
            "from ase.optimize import BFGS",
            "from ase.calculators.vasp import Vasp",
            "import shutil",
            "import copy",
            "from catlearn.optimize.mlneb import MLNEB",
            "from datetime import datetime",
            " ",
            # Initial time
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'w') as outfile:",
            "\toutfile.write(f'{dt_string}: starting job ...\\n')",
            " ",
            # Ase calculator
            "ase_calculator = Vasp(setups={'base': 'recommended'},"  # todo: update to include custom PP_dict
        ]
        for tag in [tag for tag in self.incar.tags if tag not in ['ediffg', 'n_images', 'fmax']]:
            ase_script.append(f"\t{tag}={self.incar.tags[tag]},")
        ase_script += ["\t# Kpoints",
                       f"\tkpts=({self.kpoints.num_x}, {self.kpoints.num_y}, {self.kpoints.num_z}),",
                       f"\tgamma=True,",
                       "\t)",
                       " ",
                       ]
        ase_script += [
            # Initial state optimization
            "# Optimize initial state:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting optimization of initial state ...\\n')",
            "slab = read('./optimized_structures/initial.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='initial.traj')",
            f"qn.run(fmax={abs(float(self.incar.tags['ediffg']))})",
            "shutil.copy('./initial.traj', './optimized_structures/initial.traj')",
            " ",
            # Final state optimization
            "# Optimize final state:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting optimization of final state ...\\n')",
            "slab = read('./optimized_structures/final.traj')",
            "slab.set_calculator(copy.deepcopy(ase_calculator))",
            "qn = BFGS(slab, trajectory='final.traj')",
            f"qn.run(fmax={abs(float(self.incar.tags['ediffg']))})",
            "shutil.copy('./final.traj', './optimized_structures/final.traj')",
            " ",
            # NEB
            "# Run ML-NEB:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: starting ML-NEB ...\\n')",
            "neb_catlearn = MLNEB(start='initial.traj', end='final.traj',",
            "\tase_calc=copy.deepcopy(ase_calculator),",
            f"\tn_images={self.incar.tags['n_images']},",
            "\tinterpolation='idpp',",
            "\t#restart=True,",
            "\t)",
            f"neb_catlearn.run(fmax={self.incar.tags['fmax']}, trajectory='ML-NEB.traj')",
            " ",
            # Print results
            "# Print results:",
            "now = datetime.now()",
            "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')",
            "with open('time.txt', 'a') as outfile:",
            "\toutfile.write(f'{dt_string}: ML-NEB finished\\n')",
            " ",
            "from catlearn.optimize.tools import plotneb",
            "plotneb(trajectory='ML-NEB.traj', view_path=True)"
        ]
        with open(f"{job_path}/run.py", 'w') as outfile:
            for line in ase_script:
                outfile.write(line + '\n')
