import numpy as np
def delete_vacancies( atoms ):
    """
    Delete all the vacancies in the atoms object
    """

    # Collect the index of all elements having symbol X
    del_list = [atom.index for atom in atoms if atom.symbol == "X"]
    del_list.sort()
    del_list.reverse()
    for index in del_list:
        del atoms[index]
    return atoms
