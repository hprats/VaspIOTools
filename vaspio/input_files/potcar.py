import os


def write_potcar(job_path, poscar_elements, pp_dict, pp_path):
    """Writes the POTCAR file to the current directory
    First, a list of elements is created (elements_for_potcar), which determines which
    atomic POTCAR files will be concatenated to make the final POTCAR"""
    current_element = poscar_elements[0]
    elements_for_potcar = [current_element]
    for element in poscar_elements[1:]:
        if element != current_element:
            elements_for_potcar.append(element)
            current_element = element
    cmd = 'cat'
    for element in elements_for_potcar:
        cmd += f' {pp_path}/{pp_dict[element]}/POTCAR'
    os.system(f'{cmd} > {job_path}/POTCAR')
