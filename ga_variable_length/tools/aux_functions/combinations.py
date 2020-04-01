"""
In this file it is defined the auxiliary function to get all the combinations of some length of the elements of a list. It is exactly like the tool itertools.combinations, but, as said below, itertools goes prrr (just wanted to define this function in a recurrent and dynamic programming manner).

Function:
    combinations: function that returns the combinations of some length of the elements in a list.
"""


def combinations(list_get_comb, length_combination):
    """ Generator to get all the combinations of some length of the elements of a list.  ---> itertools.combinations goes prrrr.

    :param list_get_comb: (list) List from which it is wanted to get the combination of its elements.
    :param length_combination: (int) Length of the combinations of the elements of list_get_comb.
    :return:
        * :generator: Generator with the combinations of this list.
    """
    if type(list_get_comb) != list:
        raise TypeError(
            "The parameter 'list_get_comb' must be a list.")
    if type(length_combination) != int:
        raise TypeError(
            "The parameter 'length_combination' must be a positive integer smaller than the length of the given list.")
    if length_combination <= 0 or length_combination > len(list_get_comb):
        raise ValueError(
            "The parameter 'length_combination' must be a positive integer smaller than the length of the given list.")

    # Generator to get the combinations of the indices of the list
    def get_indices_combinations(sub_list_indices, max_index):
        """ Generator that returns the combinations of the indices

        :param sub_list_indices: (list) Sub-list from which to generate ALL the possible combinations.
        :param max_index: (int) Maximum index.
        :return:
        """
        if len(sub_list_indices) == 1:  # Last index of the list of indices
            for index in range(sub_list_indices[0], max_index + 1):
                yield [index]
        elif all([sub_list_indices[-i - 1] == max_index - i for i in
                  range(len(sub_list_indices))]):  # The current sublist has reached the end
            yield sub_list_indices
        else:
            for comb in get_indices_combinations(sub_list_indices[1:],
                                                 max_index):  # Get all the possible combinations of the sublist sub_list_indices[1:]
                yield [sub_list_indices[0]] + comb
            # Advance one position and check all possible combinations
            new_sub_list = []
            new_sub_list.extend([sub_list_indices[0] + i + 1 for i in range(len(sub_list_indices))])
            for new_comb in get_indices_combinations(new_sub_list, max_index):
                yield new_comb  # Return all the possible combinations of the new list

    # Start the algorithm:
    ini_list_indices = list(range(length_combination))
    for list_indices in get_indices_combinations(ini_list_indices, len(list_get_comb) - 1):
        yield [list_get_comb[i] for i in list_indices]
