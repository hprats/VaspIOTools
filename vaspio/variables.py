from copy import deepcopy

potcars_dir_local = '/Users/install/PycharmProjects/tmc4mpo/potpaw_PBE'
name_ase_script = 'run.py'
name_submission_script = 'vasp_sub'
name_std_output = 'vasp.out'

VASP_recommended_PP = {
    # Groups 1-2 & 13-18
    'H': 'H', 'He': 'He', 'Li': 'Li_sv', 'Be': 'Be', 'B': 'B', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F', 'Ne': 'Ne',
    'Na': 'Na_pv', 'Mg': 'Mg', 'Al': 'Al', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar',
    'K': 'K_sv', 'Ca': 'Ca_sv', 'Ga': 'Ga_d', 'Ge': 'Ge_d', 'As': 'As', 'Se': 'Se', 'Br': 'Br', 'Kr': 'Kr',
    'Rb': 'Rb_sv', 'Sr': 'Sr_sv', 'In': 'In_d', 'Sn': 'Sn_d', 'Sb': 'Sb', 'Te': 'Te', 'I': 'I', 'Xe': 'Xe',
    'Cs': 'Cs_sv', 'Ba': 'Ba_sv', 'Tl': 'Tl_d', 'Pb': 'Pb_d', 'Bi': 'Bi_d', 'Po': 'Po_d', 'At': 'At', 'Rn': 'Rn',
    # d block
    'Sc': 'Sc_sv', 'Ti': 'Ti_sv', 'V': 'V_sv', 'Cr': 'Cr_pv', 'Mn': 'Mn_pv', 'Fe': 'Fe', 'Co': 'Co', 'Ni': 'Ni', 'Cu': 'Cu', 'Zn': 'Zn',
    'Y': 'Y_sv', 'Zr': 'Zr_sv', 'Nb': 'Nb_sv', 'Mo': 'Mo_sv', 'Tc': 'Tc_pv', 'Ru': 'Ru_pv', 'Rh': 'Rh_pv', 'Pd': 'Pd', 'Ag': 'Ag', 'Cd': 'Cd',
    'La': 'La', 'Hf': 'Hf_pv', 'Ta': 'Ta_pv', 'W': 'W_sv', 'Re': 'Re', 'Os': 'Os', 'Ir': 'Ir', 'Pt': 'Pt', 'Au': 'Au', 'Hg': 'Hg',
    # f block
    'Ce': 'Ce'
}

MP_recommended_PP = deepcopy(VASP_recommended_PP)
MP_recommended_PP.update({'Be': 'Be_sv', 'Mg': 'Mg_pv', 'Ti': 'Ti_pv', 'Fe': 'Fe_pv', 'Ni': 'Ni_pv', 'Cu': 'Cu_pv',
                          'Nb': 'Nb_pv', 'Mo': 'Mo_pv', 'W': 'W_pv', 'Re': 'Re_pv'})

project_PP_dict = deepcopy(VASP_recommended_PP)  # for TMC4MPO project
project_PP_dict.update({'W': 'W_pv'})