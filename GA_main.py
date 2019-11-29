#!/usr/bin/env python3
"""
In this file it is defined the main function to execute the GA with chromosomes
    with variable length.

Functions:
    GA_vl: Main function to execute the GA with chromosomes with variable length.
        
Example call:
    best_chrom, historic_fitness, pop = GA_vl(num_individuals = 100, df = mast_np, min_number_of_genes = 3, max_number_of_genes = 5, PENALIZATION = 7000)
"""
from GA_functions.aux_functions import get_master_np
from GA_functions.population import population
from GA_functions.fitness import calculate_fitness_and_order
from GA_functions.termination_criteria import check_termination_criteria
from GA_functions.next_generation import next_generation
from GA_functions.diversity import keep_diversity

#GA_VL
# Description: This function is the main function for executing the GA and the
#   only one that should be execution.
#
#   @Inputs:
#       num_individuals: Number of individuals.
#       df: numpy/pandas array in which each column represents an entity to be 
#           chosen (eg. an enterprise) and each row an item in which that  
#           entity participates. If it is a pandas df, then the colnames are   
#           taken as the entities (for presentation of the results purposes).
#       min_number_of_genes: Minimum number of genes of a chromosome.
#       max_number_of_genes: Maximum number of genes of a chromosome.
#       max_num_gen_changed_crossover: Maximum number of genes that can be 
#           taken of an individual in order to make the crossover.
#       termination_criteria: Termination criteria's function to call.
#            'goal_fitness_reached' -> Minimum fitness reached
#            'max_num_generation_reached' -> Max num generations reached
#       elitism_rate: percentage of elitism. 
#       mutation_rate: Mutation rate per individual.
#       mutation_type: It can be either:
#           'mut_gen': One of the genes of an individual is randomly changed.
#           'addsub_gen': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed_mutation: It is the number of genes that is changed, 
#           ie, the number of genes that are added, substracted or changed by 
#           new ones, when an individual is selected for mutation.
#       max_gen: (if termination criteria = max_num_generation_reached) -> 
#           -> Maximum number of generations.
#       goal_fitness: (if termination criteria = goal_fitness_reached) -> 
#           -> goal fitness
#       diversity_5_gen: If 1, then the function keep_diversity is called every
#           5 generations. In this function it is checked if there are any  
#           repeated individuals, andif so, the repetitions are substituted by 
#           a completely and randomly generated new individual.
#       min_item_per_row: If 1, it is calculated the minimum selection of items 
#           (one item per row) for a given individual. If it is 0 the maximum is 
#           calculated.
#       minimize: If 1 the fitness value is inverse normalized, ie, higher
#           value mappped to 0 and lower to 1.
#       MAX_NUM_TRANS: (fitness) Maximum number of transportists up to which
#           there is no penalization.
#       PENALIZATION: (fitness) Penalization value per each new transportist 
#           (or for their rating).
#       PERCENT: (fitness) Percentage of the total cost that penalizes each  
#           extra new trasnportists.
#       RATING_TRANS: (fitness) List with the ordered rating of each 
#           transportist.
#   @Outputs:
#       best_individual_chromosome: Best Individual's chromosome.
#       best_fit_per_gen: Array with the best individual per generation.
#       pop: Whole population.
def GA_vl(num_individuals, df, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover = 2, termination_criteria  = 'max_num_generation_reached', elitism_rate = 0.1, mutation_rate = 0.2, mutation_type = 'both', num_gen_changed_mutation = 1, max_gen = 100, goal_fitness = None, diversity_5_gen = 1, min_item_per_row = 1, minimize = 1, MAX_NUM_TRANS = 3, PENALIZATION = 0, PERCENT = 0, RATING_TRANS = []):           
    mast_np, colnames = get_master_np(df)
    if num_individuals*elitism_rate < 1 and elitism_rate != 0: #Error if the number of individuals is not big enough...
        raise ValueError('With this elitism rate, it is needed, at least, ' + str(int(1/elitism_rate)) + ' individuals per generation')
    generation_num  = 1 #Initialize the generation counter to 0
    print('Generation: ' + str(generation_num))
    best_fit_per_gen = []
    pop = population(num_individuals, mast_np, min_number_of_genes, max_number_of_genes) #Creation of the initial population
    calculate_fitness_and_order(pop, mast_np, min_item_per_row, minimize, MAX_NUM_TRANS, PENALIZATION, PERCENT, RATING_TRANS) #Calculate and order the initial population by fitness
    best_fit_per_gen.append(pop[0].fitness)
    if check_termination_criteria(termination_criteria, generation_num, max_gen, pop[0].fitness, goal_fitness, minimize):
        return pop[0].chromosome, best_fit_per_gen, pop
    while(not check_termination_criteria(termination_criteria, generation_num, max_gen, pop[0].fitness, goal_fitness, minimize)): #Run GA while the termination criteria is not met
        generation_num += 1 # Add 1 to the generation number
        print('Generation: ' + str(generation_num))
        new_chromosomes = next_generation(pop, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed_crossover, elitism_rate, mutation_rate, mutation_type, num_gen_changed_mutation) #Get the new generation chromosomes
        for i in range(len(pop)):
            pop[i].chromosome = new_chromosomes[i] #Copy new generation
        calculate_fitness_and_order(pop, mast_np, min_item_per_row, minimize, MAX_NUM_TRANS, PENALIZATION, PERCENT, RATING_TRANS) #Calculate and order the initial population by fitness
        best_fit_per_gen.append(pop[0].fitness)
        if diversity_5_gen and generation_num%5 == 0: #Every 5 generations check diversity criteria
            keep_diversity(pop, mast_np, min_number_of_genes, max_number_of_genes)
    return pop[0].chromosome, best_fit_per_gen, pop