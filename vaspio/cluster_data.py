cluster_name = 'Young'
job_scheduler = 'sge'
path_qstat_list = '/home/ucechpr'
path_submission_script_native = '/home/ucechpr/scripts/submission/native'
path_submission_script_vtst_native = '/home/ucechpr/scripts/submission/vtst_native'
path_submission_script_freq_native = '/home/ucechpr/scripts/submission/freq_native'
path_submission_script_ase = '/home/ucechpr/scripts/submission/ase'
potcars_dir_cluster = '/home/ucechpr/apps/vasp/vasp_PP_LIBRARY/potpaw_PBE'

npar = 10


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


