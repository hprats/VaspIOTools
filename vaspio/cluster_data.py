# Cluster data test
cluster_name = 'Young'
job_scheduler = 'sge'

# Names
name_submission_vasp_ase = 'run.py'
name_submission_vasp_native = 'vasp_sub'
name_vasp_std_output = 'vasp.out'

# Paths
path_qstat_list = '/home/ucechpr'
path_submission_vasp_native = '/home/ucechpr/scripts/submission/vasp_native'
path_submission_vasp_vtst = '/home/ucechpr/scripts/submission/vasp_vtst'
path_submission_vasp_ase = '/home/ucechpr/scripts/submission/vasp_ase'
pp_path = '/home/ucechpr/apps/vasp/vasp_PP_LIBRARY/potpaw_PBE'


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
