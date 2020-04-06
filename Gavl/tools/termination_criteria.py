"""
In this file it is defined the functions to check the termination criteria.
    
Functions:
    max_num_generation_reached: Function that checks if the max number of generations is reached.
    goal_fitness_reached: Function that checks if the goal fitness is reached.
    check_termination_criteria: This function selects the termination criteria
"""


def max_num_generation_reached(generation_count, max_generations):
    """ This function returns True if the maximum number of generations is reached.

    :param generation_count: (int) Number (count) of the current generation.
    :param max_generations: (int) Maximum number of generations.
    :return:
        (bool) True if the maximum number of generations is reached.
    """
    return generation_count >= max_generations


def goal_fitness_reached(generation_best_fitness, goal_fitness, minimize):
    """ This function returns True if the goal fitness is reached.

    :param generation_best_fitness: (int) Current generation best fitness.
    :param goal_fitness: (int) Goal fitness.
    :param minimize: (bool) If 1, the goal is minimize the fitness function and viceversa.
    :return:
        (bool) True if the goal fitness is reached.
    """
    if minimize:
        return generation_best_fitness <= goal_fitness
    else:
        return generation_best_fitness >= goal_fitness


def check_termination_criteria(termination_criteria_args):
    """ This function checks the termination criteria.

    :param termination_criteria_args: (dictionary) This is a dictionary with the needed arguments to check the termination criteria. It can take the values:
        * {'termination_criteria': 'goal_fitness_reached', 'goal_fitness': _ , 'generation_fitness': _ , 'minimize': _ }
        * {'termination_criteria': 'max_num_generation_reached', 'generation_goal': _ , 'generation_count': _ }
    :return:
        (bool) True if the termination criteria is met. False otherwise.
    """
    termination_criteria = termination_criteria_args['termination_criteria']
    if termination_criteria == 'max_num_generation_reached':
        generation_count = termination_criteria_args['generation_count']
        max_generations = termination_criteria_args['generation_goal']
        return max_num_generation_reached(generation_count, max_generations)
    elif termination_criteria == 'goal_fitness_reached':
        generation_best_fitness = termination_criteria_args['generation_fitness']
        goal_fitness = termination_criteria_args['goal_fitness']
        minimize = termination_criteria_args['minimize']
        return goal_fitness_reached(generation_best_fitness, goal_fitness, minimize)
    else:
        raise ValueError("The termination criteria must be either 'max_num_generation_reached' or 'goal_fitness_reached'.")
