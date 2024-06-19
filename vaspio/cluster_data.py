from copy import deepcopy

# Cluster data
cluster_name = 'Young'
job_scheduler = 'sge'

# Names
ase_vasp_script_name = 'run.py'
native_script_name = 'vasp_sub'
std_output_name = 'vasp.out'

# Paths
path_qstat_list = '/home/ucechpr'
path_submission_vasp_native = '/home/ucechpr/scripts/submission/vasp_native'
path_submission_vasp_vtst = '/home/ucechpr/scripts/submission/vasp_vtst'
path_submission_vasp_ase = '/home/ucechpr/scripts/submission/vasp_ase'
pp_path = '/home/ucechpr/apps/vasp/vasp_PP_LIBRARY/potpaw_PBE'


# Potcar dict
VASP_recommended_PP = {
    # Groups 1-2 & 13-18
    'H': 'H', 'He': 'He', 'Li': 'Li_sv', 'Be': 'Be', 'B': 'B', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F', 'Ne': 'Ne',
    'Na': 'Na_pv', 'Mg': 'Mg', 'Al': 'Al', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar',
    'K': 'K_sv', 'Ca': 'Ca_sv', 'Ga': 'Ga_d', 'Ge': 'Ge_d', 'As': 'As', 'Se': 'Se', 'Br': 'Br', 'Kr': 'Kr',
    'Rb': 'Rb_sv', 'Sr': 'Sr_sv', 'In': 'In_d', 'Sn': 'Sn_d', 'Sb': 'Sb', 'Te': 'Te', 'I': 'I', 'Xe': 'Xe',
    'Cs': 'Cs_sv', 'Ba': 'Ba_sv', 'Tl': 'Tl_d', 'Pb': 'Pb_d', 'Bi': 'Bi_d', 'Po': 'Po_d', 'At': 'At', 'Rn': 'Rn',
    # d block
    'Sc': 'Sc_sv', 'Ti': 'Ti_sv', 'V': 'V_sv', 'Cr': 'Cr_pv', 'Mn': 'Mn_pv', 'Fe': 'Fe', 'Co': 'Co', 'Ni': 'Ni',
    'Cu': 'Cu', 'Zn': 'Zn',
    'Y': 'Y_sv', 'Zr': 'Zr_sv', 'Nb': 'Nb_sv', 'Mo': 'Mo_sv', 'Tc': 'Tc_pv', 'Ru': 'Ru_pv', 'Rh': 'Rh_pv', 'Pd': 'Pd',
    'Ag': 'Ag', 'Cd': 'Cd',
    'La': 'La', 'Hf': 'Hf_pv', 'Ta': 'Ta_pv', 'W': 'W_sv', 'Re': 'Re', 'Os': 'Os', 'Ir': 'Ir', 'Pt': 'Pt', 'Au': 'Au',
    'Hg': 'Hg',
    # f block
    'Ce': 'Ce'
}

MP_recommended_PP = deepcopy(VASP_recommended_PP)
MP_recommended_PP.update({'Be': 'Be_sv', 'Mg': 'Mg_pv', 'Ti': 'Ti_pv', 'Fe': 'Fe_pv', 'Ni': 'Ni_pv', 'Cu': 'Cu_pv',
                          'Nb': 'Nb_pv', 'Mo': 'Mo_pv', 'W': 'W_pv', 'Re': 'Re_pv'})


def check_queue(job_name):
    in_queue = False
    status = None
    with open(f'{path_qstat_list}/qstat_list.txt') as infile:
        qstat_list = infile.readlines()
    for line in qstat_list:
        if line.split()[2] == job_name:
            in_queue = True
            status = line.split()[4]
            break
    return in_queue, status
