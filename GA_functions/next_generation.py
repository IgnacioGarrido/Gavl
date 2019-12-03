#!/usr/bin/env python3
"""
In this file it is defined the function to calculate the next generation.

Functions: 
    next_generation: given the population it calculates the next generation.

"""
from .selection import roulette_selection
from .pairing import pairing
from .crossover import mating
from .mutation import mutation
from random import shuffle

#NEXT_GENERATION
# Description: This function returns the chromosomes of the next generation, 
#   specified the elitism, and mutation rate.
#    
#   @Inputs:
#       pop: output of the function calculate_fitness_and_order.
#       mast_np: numpy array with the weights of each item.
#       min_number_of_genes: Minimum number of genes of the chromosome.
#       max_number_of_genes: Maximum number of genes of the chromosome.
#       max_num_gen_changed_crossover: Maximum number of genes that can be 
#           taken of an individual in order to make the crossover.
#       elitism_rate: percentage of elitism. 
#       mutation_rate: Mutation rate per individual.
#       mutation_type: It can be either:
#           'mut_gene': One of the genes of an individual is randomly changed.
#           'addsub_gene': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed_mutation: It is the number of genes that is changed, 
#           ie, the number of genes that are added, substracted or changed by
#           new ones in MUTATION.
#   @Outputs:
#       new_gen: List with the chromosomes of the new generation.
def next_generation(pop, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover = 2, elitism_rate = 0.05, mutation_rate = 0.4, mutation_type = 'both', num_gen_changed_mutation = 1):
    new_gen = [] #List with the chromosomes of the next generation
    num_elitism = int(len(pop)*elitism_rate) #Number of individuals chosen as elite
    num_cross = len(pop) - num_elitism #Number of individuals chosen for crossover
    if (num_cross % 2) == 1: #Odd number of Individuals for crossover
        num_elitism += 1 #Add one individual to elitism
        num_cross -= 1 #Remove one from crossover
    #Elite:
    new_gen.extend(x.chromosome for x in pop[:num_elitism]) #Add elite
    #Crosover:
    selected_ind = roulette_selection(pop, num_cross) #1. Roulette selection
    paired_ind = pairing(selected_ind) #2. Perform the pairing
    new_cross_ind = mating(paired_ind, pop, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover) #3. Get the new chromosomes already crossovered
    new_gen.extend(new_cross_ind) #4. Add crossovered individuals
    #Mutation:
    num_mutated = int(len(pop)*mutation_rate) #Number of individuals chosen to be mutated
    list_IDs = list(range(num_elitism, num_cross + num_elitism)) #All the IDs of the individuals that may be mutated (not the elite)
    shuffle(list_IDs)
    mut_IDs = list_IDs[:num_mutated] #Get list of IDs of chromosomes to be mutated
    for i in mut_IDs:
        new_gen[i] = mutation(new_gen[i], mast_np, mutation_type, num_gen_changed_mutation, min_number_of_genes, max_number_of_genes)
    return new_gen