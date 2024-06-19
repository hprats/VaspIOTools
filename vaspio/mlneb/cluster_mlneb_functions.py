import sys

from vaspio.cluster_functions import submit_job
from vaspio.mlneb.mlneb_status import get_mlneb_status

def get_import_lines():
    lines = [
        "from ase.io import read\n",
        "from ase.optimize import BFGS\n",
        "from ase.calculators.vasp import Vasp\n",
        "import shutil\n",
        "import copy\n",
        "from catlearn.optimize.mlneb import MLNEB\n",
        "from datetime import datetime\n",
        " \n"
    ]
    return lines


def get_calculator_lines(lines):
    index_start = 0
    index_end = 0
    for i, line in enumerate(lines):
        if "ase_calculator = " in line:
            index_start = i
        if "# Kpoints" in line:
            index_end = i + 5
    return lines[index_start:index_end]


def get_initial_lines(lines):
    index_fmax = lines.index("shutil.copy('./initial.traj', './optimized_structures/initial.traj')\n") - 1
    fmax = lines[index_fmax].split("=")[1].split(")")[0]
    lines = [
        "# Optimize initial state:\n",
        "now = datetime.now()\n",
        "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
        "with open('time.txt', 'a') as outfile:\n",
        "\toutfile.write(f'{dt_string}: starting optimization of initial state ...\\n')\n",
        "slab = read('./initial.traj')\n",
        "slab.set_calculator(copy.deepcopy(ase_calculator))\n",
        "qn = BFGS(slab, trajectory='initial.traj')\n",
        f"qn.run(fmax={fmax})\n",
        "shutil.copy('./initial.traj', './optimized_structures/initial.traj')\n",
        " \n"
    ]
    return lines


def get_final_lines(lines, restart=False):
    index_fmax = lines.index("shutil.copy('./final.traj', './optimized_structures/final.traj')\n") - 1
    fmax = lines[index_fmax].split("=")[1].split(")")[0]
    lines = [
        "# Optimize final state:\n",
        "now = datetime.now()\n",
        "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
        "with open('time.txt', 'a') as outfile:\n",
        "\toutfile.write(f'{dt_string}: starting optimization of final state ...\\n')\n"
    ]
    if restart:
        lines += ["slab = read('./final.traj')\n"]
    else:
        lines += ["slab = read('./optimized_structures/final.traj')\n"]
    lines += [
        "slab.set_calculator(copy.deepcopy(ase_calculator))\n",
        "qn = BFGS(slab, trajectory='final.traj')\n",
        f"qn.run(fmax={fmax})\n",
        "shutil.copy('./final.traj', './optimized_structures/final.traj')\n",
        " \n"
    ]
    return lines


def get_mlneb_lines(lines, restart=False):
    index_n_images = lines.index("neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n") + 2
    n_images = lines[index_n_images].split("=")[1].split(",")[0]
    index_fmax = lines.index("neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n") + 6
    fmax = lines[index_fmax].split("=")[1].split(",")[0]
    lines = [
        "# Run ML-NEB:\n",
        "now = datetime.now()\n",
        "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
        "with open('time.txt', 'a') as outfile:\n",
        "\toutfile.write(f'{dt_string}: starting ML-NEB ...\\n')\n",
        "neb_catlearn = MLNEB(start='initial.traj', end='final.traj',\n",
        "\tase_calc=copy.deepcopy(ase_calculator),\n",
        f"\tn_images={n_images},\n",
        "\tinterpolation='idpp',\n"
    ]
    if restart:
        lines += ["\trestart=True,\n"]
    else:
        lines += ["\t#restart=True,\n"]
    lines += [
        "\t)\n",
        f"neb_catlearn.run(fmax={fmax}, trajectory='ML-NEB.traj')\n",
        " \n"
    ]
    return lines


def get_print_lines():
    lines = [
        "# Print results:\n",
        "now = datetime.now()\n",
        "dt_string = now.strftime('%d/%m/%Y %H:%M:%S')\n",
        "with open('time.txt', 'a') as outfile:\n",
        "\toutfile.write(f'{dt_string}: ML-NEB finished\\n')\n",
        " \n",
        "from catlearn.optimize.tools import plotneb\n",
        "plotneb(trajectory='ML-NEB.traj', view_path=True)\n"
    ]
    return lines


def continue_mlneb(job_path, job_name, job_status, job_scheduler, path_submission_script, name_submission_script):
    with open(f"{job_path}/run.py", "r") as infile:
        lines = infile.readlines()
    ase_script = get_import_lines()
    ase_script += get_calculator_lines(lines=lines)
    if job_status == "max wallclock":
        ase_script += get_mlneb_lines(lines=lines, restart=True)
        ase_script += get_print_lines()
    elif job_status == "final state":
        ase_script += get_final_lines(lines=lines, restart=True)
        ase_script += get_mlneb_lines(lines=lines, restart=False)
        ase_script += get_print_lines()
    elif job_status == "initial state":
        ase_script += get_initial_lines(lines=lines)
        ase_script += get_final_lines(lines=lines, restart=False)
        ase_script += get_mlneb_lines(lines=lines, restart=False)
        ase_script += get_print_lines()
    else:
        sys.exit(f"{job_path}: ERROR, cannot determine where to restart the job, check manually")
    with open(f"{job_path}/run.py", 'w') as outfile:
        for line in ase_script:
            outfile.write(line)
    submit_job(job_path, job_name, job_scheduler, path_submission_script, name_submission_script)
