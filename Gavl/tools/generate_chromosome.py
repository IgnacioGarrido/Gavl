"""
In this file it is defined the functions to create a new individual (chromosome).

Functions:
    generate_chromosome: Main function.
"""
import random


def generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes, repeated_genes_allowed):
    """ Function called to create a new individual (its chromosome). It randomly chooses its length (between min_length_chromosome and min_length_chromosome), and it randomly chooses genes among the list of possible_genes.

    :param min_length_chromosome: (int) Minimum allowed length of the chromosome.
    :param max_length_chromosome: (int) Maximum allowed length of the chromosome.
    :param possible_genes: (list of ...) List with the all the possible values that the genes can take.
    :param repeated_genes_allowed: (bool) It is a boolean that indicates whether the genes can be repeated in the chromosome (repeated_genes_allowed = 1) or they cannot be repeated (repeated_genes_allowed = 0).
    :return:
        * (list of genes) List that represents the chromosome.
    """
    # Choose a random number of genes
    number_of_genes = random.randrange(min_length_chromosome, max_length_chromosome + 1)
    # Create new chromosome:
    if repeated_genes_allowed:
        chromosome = random.choices(possible_genes, weights=None, k=number_of_genes)
        return chromosome
    else:
        possible_genes_aux = possible_genes.copy()
        random.shuffle(possible_genes_aux)
        return possible_genes_aux[:number_of_genes]
