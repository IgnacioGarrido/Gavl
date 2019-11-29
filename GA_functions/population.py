#!/usr/bin/env python3
"""
In this file it is defined the functions to create the population.

Functions:
    create_individual: Creation of a single individual.
    population: Creation of the whole population.
"""
from .aux_functions import check_valid_chromosome
from .Individual import Individual
from random import randrange
  
#CREATE_INDIVIDUAL
# Description: This function creates individuals with a random and variable 
#   length of the genes between the specified limits. Each individual is a set 
#   of ints.
#
#   @Inputs:
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes
#       mast_np: numpy array with the weights of each item.
#   @Outputs:
#       Individual(chrom, 0): individual of the class Individual.
def create_individual(min_number_of_genes, max_number_of_genes, mast_np):
    num_elem_to_choose = mast_np.shape[1] #Number of elements from which to choose
    pool = list(range(num_elem_to_choose))
    num_genes = randrange(min_number_of_genes,max_number_of_genes + 1,1)
    chrom=[pool.pop(randrange(len(pool))) for x in range(num_genes)]
    if check_valid_chromosome(chrom, mast_np): #Check if it is a valid individual
        return Individual(chrom, 0, 0) #Temporal value of the fitness = 0
    else:
        return create_individual(min_number_of_genes, max_number_of_genes, mast_np)

#POPULATION
# Description: This function creates a population of valid individuals.
#
#   @Inputs:
#       num_individuals: Number of individuals.
#       mast_np: numpy array with the weights of each item.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes
#   @Outputs:
#       list_individuals: List with the individuals.
def population(num_individuals, mast_np, min_number_of_genes, max_number_of_genes):
    return [create_individual(min_number_of_genes, max_number_of_genes, mast_np) for i in range(num_individuals)]