#!/usr/bin/env python3
"""
In this file it is defined the functions to perform crossover.

Functions: 
    cross_ind: this functio calculated the crossover of two individuals.
    mating: this function calculates the crossover of ALL the paired 
        individuals with the help of the function cross_ind.

IF A DIFFERENT CROSSOVER IS NEEDED, IT SHOULD BE PROGRAMMED IN THE FUNCTION cross_ind
"""
from .aux_functions import check_valid_chromosome
from random import shuffle
from itertools import combinations

#CROSS_IND
# Description: This function calculates a crossover that gives a valid
#   chromosome between two individuals. It selects a random number (between 1 
#   and max_num_gen_changed) of genes of the individuals a and b (note that the
#   number of genes from individual a to individual b may be different than the
#   number of genes changed from individual b to individual a) that are 
#   changed, and calculates if there is any possible crossover of that size. 
#   If so, one of the possible crossovers (selected randomly) is performed. If 
#   there is no possible crossover of that size, it is tested another different
#   combination of sizes. 
#
#   @Inputs:
#       chrom_a: Individual A's chromosome.
#       chrom_b: Individual B's chromosome.
#       mast_np: numpy array with the weights of each item.
#       min_number_of_genes: Minimum number of genes of the chromosome.
#       max_number_of_genes: Maximum number of genes of the chromosome.
#       max_num_gen_changed_crossover: Maximum number of genes that can be 
#           taken of an individual in order to make the crossover.
#   @Outputs:
#       crossed_1, crossed_2: The two already crossed individuals.
def cross_ind(chrom_a, chrom_b, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover):
    num_gen_a_interchanged = list(map(lambda x: x+1, list(range(max_num_gen_changed_crossover))))
    shuffle(num_gen_a_interchanged) #random number of genes of a that are interchanged
    num_gen_b_interchanged = list(map(lambda x: x+1, list(range(max_num_gen_changed_crossover)))) 
    shuffle(num_gen_b_interchanged) #random number of genes of b that are interchanged
    for num_a in num_gen_a_interchanged: #Take a random number of genes of individual a to interchange
        for num_b in num_gen_b_interchanged: #Take a random number of genes of individual b to interchange
            pairs_a = list(combinations(chrom_a, num_a)) #Get all the possible combinations of num_a elements in ind_a
            pairs_b = list(combinations(chrom_b, num_b)) #Get all the possible combinations of num_b elements in ind_b
            if not ((len(pairs_a) == 0) or (len(pairs_b) == 0)): #If there are no possible combination (more elements to interchange than the length of the list) -> THIS IS NOT SUPPOSED TO HAPPEN ANY TIME IF THE CROSSOVER IS DONE IN AN INTELLIGENT WAY...
                if not(((len(chrom_a)-num_a+num_b) < min_number_of_genes) or ((len(chrom_b)-num_b+num_a) < min_number_of_genes) or ((len(chrom_a)-num_a+num_b) > max_number_of_genes) or ((len(chrom_b)-num_b+num_a) > max_number_of_genes)): #Too big/small chromosomes would be formed   
                    shuffle(pairs_a)
                    shuffle(pairs_b)
                    for cross_a in pairs_a: #Take a random crossover possibility for individual a
                        if not (True in [i in chrom_b for i in cross_a]): #Some of the genes of individual a that is going to be interchanged are NOT already in individual b (wouldn't be a crossover)           
                            for cross_b in pairs_b: #Take a random crossover possibility for individual b
                                if not (True in [i in chrom_a for i in cross_b]): #Some of the genes of individual b that is going to be interchanged are NOT already in individual a (wouldn't be a crossover)
                                    crossed_1 = chrom_a.copy() #Crossed individual 1
                                    crossed_2 = chrom_b.copy() #Crossed individual 2
                                    for gen in cross_a: #Eliminate the genes of that were selected for crossover from ind a from crossed_1 and add them to crossed_2
                                        crossed_1.remove(gen) 
                                        crossed_2.append(gen)        
                                    for gen in cross_b: #Eliminate the genes of that were selected for crossover from ind b from crossed_2 and add them to crossed_1
                                        crossed_2.remove(gen) 
                                        crossed_1.append(gen)              
                                    if check_valid_chromosome(crossed_1, mast_np) and check_valid_chromosome(crossed_2, mast_np): #If the crossovers result in valid individuals, return them
                                        return crossed_1, crossed_2
    return chrom_a, chrom_b #If there are not possible crossovers, just return the two individuals
  
    
#MATING
# Description: This function ruturns the mated couples of individuals when it 
#   is possible to make this mating. It makes a random selection of one of the 
#   possible matings given by calc_possible_crossovers. If it is not posible to
#   make the mating because all the combinations give invalid individuals, 
#   then the two intended to be paired individuals are returned.
#
#   @Inputs:
#       paired_ind: List of tuples with the IDs of the mated of the 
#           individuals. It is the output of the function pairing.
#       ordered_pop: Ordered population. It is like the output of the function
#           calculate_fitness_and_order.
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes of the chromosome.
#       max_number_of_genes: Maximum number of genes of the chromosome.
#       max_num_gen_changed_crossover: Maximum number of genes that can be 
#           taken of an individual in order to make the crossover.
#   @Outputs:
#       list_new_ind: List with the new chromosomes. 
def mating(paired_ind, ordered_pop, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover):
    list_new_ind = []
    for (ind_a, ind_b) in paired_ind:
        chrom_a = ordered_pop[ind_a].chromosome.copy()
        chrom_b = ordered_pop[ind_b].chromosome.copy()
        new_chrom_a, new_chrom_b = cross_ind(chrom_a, chrom_b, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover)
        list_new_ind.append(new_chrom_a) #New individuals with temporary fitness value of 0
        list_new_ind.append(new_chrom_b)         
    return list_new_ind