"""
In this file it is defined the function to perform random pairing given a selection (like the one given by roulette wheel).
    
Functions: 
    pairing: Given a selection (list with IDs), it performs a pairing of individuals.
"""
import random


def pairing(list_selected_ind):
    """ This function performs random pairing. NOTE that this function accepts repetitions and elements may be paired with themselves. However, this rarely happens and it can be used as a elite process.

    :param list_selected_ind: List with the IDs of the selected individuals. It is the output of the function roulette_selection.
    :return:
        * :paired_ind: List of tuples with the IDs of the paired individuals.
    """
    list_sel = list_selected_ind.copy()
    if len(list_sel) % 2 == 1:
        list_sel.pop()  # If odd length, pop last individual of the list
    random.shuffle(list_sel)
    paired_ind = []
    while len(list_sel) > 0:
        paired_ind.append((list_sel.pop(),
                           list_sel.pop()))  # There can be repetition of individuals (an individual may be paired with itself)
    return paired_ind

