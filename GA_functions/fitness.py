#!/usr/bin/env python3
"""
In this file it is defined the fitness function

Functions:
    fitness: It calculates the fitness for a given chromosome
    calculate_fitness_and_order: This function loops over the whole population 
    	for calculating the fitness of each individual, it normalizes it and
    	then it orders the whole population from best to worst fitness.
  
IF NEEDED OTHER FITNESS, THIS FILE SHOULD BE CHANGED!

"""
from .aux_functions import calculate_min_row, calculate_max_row

#FITNESS_CALCULATION
# Description: This function calculates the fitness of an individual.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the weights of each item.
#       min_item_per_row: If 1, it is calculated the minimum selection of items 
#       (one item per row) for a given individual. If it is 0 the maximum is 
#       calculated.
#       MAX_LENGTH_CHROM: Maximum length of the chromosome up to which there is 
#           no penalization.
#       PENALIZATION_LENGTH: (fitness) Penalization value per each new 
#           transportist over the limit MAX_NUM_TRANS.
#       PERCENT: (fitness) Percentage of the total cost that penalizes each  
#           extra new trasnportists.
#       PENALIZATION_RATING: (fitness) Penalization value per each bad 
#           transportist in the rating.
#       RATING_TRANS: List with the ordered rating of each transportist.
#   @Outputs:
#       fitness_val: Value of the fitness  
def fitness(chrom, mast_np, min_item_per_row = 1, MAX_LENGTH_CHROM = 3, PENALIZATION_LENGTH = 0, PERCENT = 0, PENALIZATION_RATING = 0, RATING_TRANS = []):
    if min_item_per_row:
        price = calculate_min_row(chrom, mast_np)
    else:
        price = calculate_max_row(chrom, mast_np)
    if RATING_TRANS == []: #If empty list -> the array of 0s
        RATING_TRANS = [0]*mast_np.shape[1]
    #pen_num_trans = PENALIZATION_LENGTH*((max(MAX_NUM_TRANS,len(chrom))-MAX_NUM_TRANS)**2) #Penalization for the number of transportists, squared
    pen_num_trans = PENALIZATION_LENGTH*(max(MAX_LENGTH_CHROM,len(chrom))-MAX_LENGTH_CHROM) #Penalization for the number of transportists, linear
    pen_percent = price*(max(MAX_LENGTH_CHROM,len(chrom))-MAX_LENGTH_CHROM)*PERCENT #Penalization as a percentage of the total cost for each extra transportist
    pen_rating = PENALIZATION_RATING*sum([RATING_TRANS[i] for i in chrom])/len(chrom) #Penalization for the rating of the transportists
    fitness_val = price + pen_num_trans + pen_percent + pen_rating 
    return fitness_val


#CALCULATE_FITNESS_AND_ORDER
# Description: This function calculates the fitness of the Individuals, it 
#   normalizes it, and then it order the population from best to worst fitness.
#   Note that it order the input value pop, and it doesn´t return anything.
#
#   @Inputs:
#       pop: population (list with Individuals).
#       mast_np: numpy array with the weights of each item.
#       min_item_per_row: If 1, it is calculated the minimum selection of items 
#       (one item per row) for a given individual. If it is 0 the maximum is 
#       calculated.
#       minimize: If 1 the fitness value is inverse normalized, i.e, higher
#           values mappped to 0 and lower to 1. Otherwise it is normaly
#           normalized.
#       MAX_LENGTH_CHROM: Maximum lenght of the chromosome up to which there is 
#           no penalization.
#       PENALIZATION_LENGTH: (fitness) Penalization value per each new 
#           transportist over the limit MAX_NUM_TRANS.
#       PERCENT: (fitness) Percentage of the total cost that penalizes each  
#           extra new trasnportists.
#       PENALIZATION_RATING: (fitness) Penalization value per each bad 
#           transportist in the rating.
#       RATING_TRANS: List with the ordered rating of each transportist.
#   @Outputs:
#       None
def calculate_fitness_and_order(pop, mast_np, min_item_per_row = 1, minimize = 1, MAX_LENGTH_CHROM = 3, PENALIZATION_LENGTH = 0, PERCENT = 0, PENALIZATION_RATING = 0, RATING_TRANS = []):
    for ind in pop:
        chrom = ind.chromosome
        ind.fitness = fitness(chrom, mast_np, min_item_per_row, MAX_LENGTH_CHROM, PENALIZATION_LENGTH, PERCENT, PENALIZATION_RATING, RATING_TRANS)
    max_v = max(list(map(lambda x: x.fitness, pop)))
    min_v = min(list(map(lambda x: x.fitness, pop)))
    if max_v != min_v: #There is no 0-division
        if minimize: 
            for ind in pop:
                ind.normalized_fitness = 1-(ind.fitness-min_v)/(max_v-min_v)
        else:
            for ind in pop:
                ind.normalized_fitness = (ind.fitness-min_v)/(max_v-min_v) 
    else: #There would be 0-division
        print ('In this generation the whole population has converged to the same fitness value.')
        for ind in pop:
            ind.normalized_fitness = 1
    pop.sort(key=lambda x: x.normalized_fitness, reverse=True) #Sort by normalized fitness