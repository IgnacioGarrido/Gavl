"""
In this file it is shown an example of the knapsack problem (the one explained in the README).
"""
import Gavl

ga = Gavl.Gavl()  # Initialize

# Compulsory hyperparameters

ga.set_hyperparameter('size_population', 50)  # Set the size of the population to 50 individuals
ga.set_hyperparameter('min_length_chromosome', 1)  # Set the minimum length of the individual to 1 genes
ga.set_hyperparameter('max_length_chromosome',
                      9)  # Set the maximum length of the individual to 9 genes (there are 9 items)

# Definition of the fitness:
# The fitness function receives the chromosome, which is a list with the genes. It must return the fitness value (number).

MAX_WEIGTH = 15  # Maximum allowed weight
PENALIZATION = 10  # Penalization for each kg over MAX_WEIGTH

# Dictionary ---> {'item': (price, weight)}
prices_weights = {
    'pen': (5, 3),
    'pencil': (4, 2),
    'food': (7, 6),
    'rubber': (3, 1),
    'book': (10, 9),
    'scissors': (6, 3),
    'glasses': (7, 5),
    'case': (7, 7),
    'sharpener': (2, 1)
}


def fun_fitness(chromosome):
    """ Definition of the fitness function.

    :param chromosome: (list of genes) Chromosome of the individual that is being analyzed.
    """
    fitness = 0  # Cumulative fitness
    sum_weights = 0  # Cumulative weights
    for item in chromosome:
        fitness += prices_weights[item][0]
        sum_weights += prices_weights[item][1]
    fitness -= PENALIZATION * (max(sum_weights, MAX_WEIGTH) - MAX_WEIGTH)  # Add penalization
    return fitness


ga.set_hyperparameter('fitness', fun_fitness)  # Set the fitness function

# Definition of the values that the genes can take

possible_genes = ['pen', 'pencil', 'food', 'rubber', 'book', 'scissors', 'glasses', 'case', 'sharpener']
ga.set_hyperparameter('possible_genes',
                      possible_genes)  # The possible values that the genes can take are the items' names

# # In case it is wanted to NOT allow chromosomes that overpass the maximum allowed weight it can be done the following (notice that it has preferred to define a penalization term instead):
#
# def check_valid_individual(chromosome):
#     """ Definition of a valid individual.
#
#     :param chromosome: (list of genes) Chromosome of the individual that is being analyzed.
#     """
#     sum_weights = 0  # Cumulative weights
#     for item in chromosome:
#         sum_weights += prices_weights[item][1]
#     if sum_weights > MAX_WEIGTH:
#         return False
#     else:
#         return True
#
# ga.set_hyperparameter('check_valid_individual', check_valid_individual)

# Set other parameters

ga.set_hyperparameter('termination_criteria', {
    'max_num_generation_reached': 30})  # Set the termination criteria to maximum number of generations (in this case 30 generations).
ga.set_hyperparameter('repeated_genes_allowed',
                      0)  # In this case the repeated genes are NOT allowed ---> Each item can be taken only once
ga.set_hyperparameter('minimize', 0)  # In this case the goal is to maximize the fitness
ga.set_hyperparameter('elitism_rate', 0.1)  # Set an elitism rate of 0.1
ga.set_hyperparameter('mutation_rate', 0.2)  # Set a mutation rate of 0.2
ga.set_hyperparameter('mutation_type', 'both')  # Set the mutation type
ga.set_hyperparameter('keep_diversity', 5)  # Set a mechanism to keep the diversity every 5 generations
ga.set_hyperparameter('show_progress',
                      1)  # Show progress of the genetic algorithm (this will print the generation number)

best_individual = ga.optimize()  # Optimize

best_individual, population, historic_fitness = ga.get_results()  # Get other results
