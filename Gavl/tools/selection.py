"""
In this file it is defined the function to perform a roulette wheel selection.
    
Functions: 
    roulette_selection: Given a list with the tuples (id_individual, normalized_fitness), this function calculates the a roulette wheel selection based in the normalized fitness.
"""
import random


def roulette_selection(population, minimize, num_selected_ind):
    """ This function returns a list with the ids of the individuals selected by the roulette wheel selection. Note that the normalized fitness of the population must be calculated before calling this function (call the method Gavl._Population__calculate_normalized_fitness).

    :param population: (list of Individuals) This is a list of individuals (see class Individual).
    :param minimize: (int) Int that represents if the goal is minimizing the fitness (minimize = 1) or maximizing it (minimize = 0).
    :param num_selected_ind: (int) number of individuals to be selected.
    :return:
        * :list_selected_individuals: (list of str) List with the ids of the selected individuals. Note that there can be repeated individuals.
    """
    if minimize:
        list_ids_normalizedfitness = [(ind._id, ind.inverse_normalized_fitness_value) for ind in population]
    else:
        list_ids_normalizedfitness = [(ind._id, ind.normalized_fitness_value) for ind in population]
    list_selected_individuals = []  # List that will contain the ids of the selected individuals.
    sum_fit = sum(map(lambda x: x[1], list_ids_normalizedfitness))  # Sum of the INVERSE normalized fitness
    for _ in range(num_selected_ind):
        population_cumulative_fitness = 0  # Population cumulative fitness
        selected_cumulative_fitness = random.random() * sum_fit  # Cumulative fitness of the selected individual
        for ind in list_ids_normalizedfitness:
            population_cumulative_fitness += ind[1]
            if selected_cumulative_fitness <= population_cumulative_fitness:
                list_selected_individuals.append(ind[0])  # Add the id of the individual to the population
                break
    return list_selected_individuals
