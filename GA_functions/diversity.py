#!/usr/bin/env python3
"""
In this file it is defined a function to keep the diversity.

Functions:
    keep_diversity: Function called to keep the diversity.
"""
from .population import create_individual

#KEEP_DIVERSITY
# Description: This function is called when it is wanted to do a great emphasis
#   in the diversity of the problem. When an individual is repeated in the
#   population, this repeated individual is substituted by a completely new 
#   and randomly generated individual.
#    
#   @Inputs:
#       ordered_pop: output of the function calculate_fitness_and_order.
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes of the chromosome.
#       max_number_of_genes: Maximum number of genes of the chromosome.
#   @Outputs:
#       None
def keep_diversity(ordered_pop, mast_np, min_number_of_genes, max_number_of_genes):
    for x in ordered_pop: #Sort the individuals
        x.chromosome.sort() 
    list_chrom = [x.chromosome for x in ordered_pop]
    for i in range(1,len(list_chrom)):
        if ordered_pop[i].chromosome in list_chrom[:i]:
            ordered_pop[i].chromosome = create_individual(min_number_of_genes, max_number_of_genes, mast_np).chromosome.copy()    