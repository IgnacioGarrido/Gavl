#!/usr/bin/env python3
"""
In this file it is defined the function to perform a roulette wheel selection.
    
Functions: 
    roulette_selection: Given a population, this fucntion calculates the a 
        roulette wheel selection based in the normalized fitness.
"""
from random import random as rnd

#ROULETTE_SELECTION
# Description: This function returns a list with all the individuals selected 
#   by the roulette wheel selection.
#    
#   @Inputs:
#       pop: Population - output of the function calculate_fitness_and_order.
#       num_selected_ind: number of selected individuals.
#   @Outputs:
#       list_selected_ind: List with the IDs of the selected individuals.
def roulette_selection(pop, num_selected_ind):
    sumFit = sum(map(lambda x: x.normalized_fitness, pop)) #Sum of the normalized fitness
    list_selected_ind = []
    for i in range(num_selected_ind):
        cumSumFit = 0 #Cumulative sum
        rndFit = rnd()*sumFit #Sum fitness of the selected individual
        for i in range(len(pop)):
            cumSumFit += pop[i].normalized_fitness
            if rndFit <= cumSumFit:
                list_selected_ind.append(i)
                break
    return list_selected_ind