"""
In this file it is defined the functions to perform crossover.

Functions: 
    cross_individuals: this function performs the crossover of two given individuals.
    mating: this function calculates the crossover of ALL the paired individuals with the help of the function cross_ind.
"""
import random
from .aux_functions.combinations import combinations


def cross_individuals(chromosome_a, chromosome_b, min_length_chromosome, max_length_chromosome, repeated_genes_allowed,
                      check_valid_individual):
    """  This function calculates the crossover between two different individuals. It selects a random number of genes of the individuals a and b that are going to be interchanged (note that the number of genes from individual a to individual b may be different than the number of genes changed from individual b to individual a), and calculates if there is any possible crossover of that size (with the function check_valid_individual). If so, one of the possible crossovers (selected randomly) is performed and the resulting new chromosomes are returned. If there is no possible crossover of that size, it is tested another different combination of sizes. Notice that the function check_valid_individual is used to test the created individuals, and if it is done 2000 unsuccessful crossovers, it is taken as an impossible to couple pair of individuals and their original chromosomes are returned.

    :param chromosome_a: (Individual) Individual A's chromosome.
    :param chromosome_b: (Individual) Individual B's chromosome.
    :param min_length_chromosome: (int) Minimum number of genes of the chromosome.
    :param max_length_chromosome: (int) Maximum number of genes of the chromosome.
    :param repeated_genes_allowed: (int) Boolean that indicates if an individual can have repeated genes. It can take the values 1 (repeated genes allowed) or 0 (repeated genes not allowed).
    :param check_valid_individual: (function) Function that receives a chromosome and returns a boolean that indicates whether the chromosome makes a valid individual (True) or not (False) according to some criteria.
    :return:
        * :crossed_a: Crossed individual a.
        * :crossed_b: Crossed individual b.
    """
    if repeated_genes_allowed:
        genes_a = chromosome_a  # List of genes of chromosome a that can be crossed
        genes_b = chromosome_b  # List of genes of chromosome a that can be crossed
    else:
        genes_a = [gen for gen in chromosome_a if
                   gen not in chromosome_b]  # List of genes of chromosome_a that are not in chromosome_b.
        genes_b = [gen for gen in chromosome_b if
                   gen not in chromosome_a]  # List of genes of chromosome_b that are not in chromosome_a.
    random.shuffle(genes_a)  # randomize the order in which the genes of a are taken
    random.shuffle(genes_b)  # randomize the order in which the genes of a are taken
    # Start the iterations
    number_of_genes_to_change_a = list(range(1, len(genes_a) + 1))  # List with the number of possible genes to change of chromosome a
    random.shuffle(number_of_genes_to_change_a)  # randomize the number of genes that will be taken from a
    count_crossover_tried = 0
    for num_a in number_of_genes_to_change_a:  # Take a random number of genes of individual a to interchange
        number_of_genes_to_change_b = list(range(1, len(genes_b) + 1))  # List with the number of possible genes to change of chromosome a
        random.shuffle(number_of_genes_to_change_b)  # randomize the number of genes that will be taken from b
        for num_b in number_of_genes_to_change_b:  # Take a random number of genes of individual b to interchange
            if min_length_chromosome <= len(
                    chromosome_a) - num_a + num_b <= max_length_chromosome and min_length_chromosome <= len(
                    chromosome_b) - num_b + num_a <= max_length_chromosome:  # The resulting individuals are within the limits min_length_chromosome and max_length_chromosome.
                genes_a_combinations = combinations(genes_a, num_a)  # Genes of a to interchange
                for genes_change_a in genes_a_combinations:
                    genes_b_combinations = combinations(genes_b, num_b)  # Genes of a to interchange
                    for genes_change_b in genes_b_combinations:
                        count_crossover_tried += 1
                        crossed_a = chromosome_a.copy()
                        crossed_b = chromosome_b.copy()
                        for gen in genes_change_a:
                            crossed_a.remove(gen)
                            crossed_b.append(gen)
                        for gen in genes_change_b:
                            crossed_b.remove(gen)
                            crossed_a.append(gen)
                        if not check_valid_individual(crossed_b):
                            break  # No combination of genes of b will result in a valid chromosome for a ---> Get other combination of genes of a
                        if check_valid_individual(crossed_a):  # else try with other combination of genes
                            return crossed_a, crossed_b
                        if count_crossover_tried >= 2000:
                            return chromosome_a, chromosome_b  # If it has been tested 2000 unsuccessful crossovers, then break
            else:
                break  # The resulting individuals would NOT be within the limits min_length_chromosome and max_length_chromosome. ---> SOLUTION: simply try other combination.
    return chromosome_a, chromosome_b  # If there are not possible crossovers, just return the two individuals


def mating(list_of_paired_ind, min_length_chromosome, max_length_chromosome, repeated_genes_allowed,
           check_valid_individual):
    """ This function returns the mated couples of individuals when it is possible to make this mating. If it is not possible to make the mating because all the combinations give invalid individuals (if repeated_genes_allowed = 0 and the two individuals are exactly the same genes), then the two intended to be paired individuals are returned.

    :param list_of_paired_ind: (list of tuples of Individuals) List of tuples where each tuple represents two paired individuals. It has the form [(Individual_a, Individual_b), (Individual_c, Individual_d), ...]; having paired in this example the individual a with the individual b, and the individual c with the individual d. Note that the tuples contain objects of the class Individual.
    :param min_length_chromosome: (int) Minimum number of genes of the chromosome.
    :param max_length_chromosome: (int) Maximum number of genes of the chromosome.
    :param repeated_genes_allowed: (int) Boolean that indicates if an individual can have repeated genes. It can take the values 1 (repeated genes allowed) or 0 (repeated genes not allowed).
    :param check_valid_individual: (function) Function that receives a chromosome and returns a boolean that indicates whether the chromosome makes a valid individual (True) or not (False) according to some criteria.
    :return:
        * :crossed_individuals: (list of chromosomes) List with the crossed individuals' chromosomes.
    """
    crossed_individuals = []  # Output ---> List with the crossed individuals.
    for (chromosome_a, chromosome_b) in list_of_paired_ind:
        crossed_a, crossed_b = cross_individuals(chromosome_a=chromosome_a, chromosome_b=chromosome_b,
                                                 min_length_chromosome=min_length_chromosome,
                                                 max_length_chromosome=max_length_chromosome,
                                                 repeated_genes_allowed=repeated_genes_allowed,
                                                 check_valid_individual=check_valid_individual)
        crossed_individuals.append(crossed_a)
        crossed_individuals.append(crossed_b)
    return crossed_individuals
