# Ignacio Garrido Botella
# This file includes all the functions needed to perform the GA operations:
#   Selection
#   Pairing
#   Mating
#   Mutation
#   Fitness calculation
#   Check termination criteria
# For the completion of this code it has been very useful the work of 
#   Cahit Bartu Yazici, showed in the post: https://towardsdatascience.com/continuous-genetic-algorithm-from-scratch-with-python-ff29deedd099
import numpy as np
import pandas as pd
from random import random as rnd
from random import randrange, seed, shuffle, choice
from itertools import combinations

seed(42)

df_prices = pd.read_csv('/Users/Ignacio/Documents/Universidad/MUIT/Tercero/TFM/Datasets/Optimization/final/prices_20172819.csv')
colnames = df_prices.columns
cantidad = df_prices.iloc[:,4]
for colname in colnames:
    df_prices[colname] = df_prices[colname]*cantidad
df_prices = df_prices.iloc[:200,5:] #Eliminate route columns
mast_np = df_prices.to_numpy()

del df_prices, colnames, colname, cantidad

#%% CLASS INDIVIDUAL:

# This is the class of the individuals. They have the attributes chromosome,
#   fitness and normalized fitness. 
class Individual:
  def __init__(self, chromosome, fitness, normalized_fitness):
    self.chromosome = chromosome
    self.fitness = fitness
    self.normalized_fitness = normalized_fitness

#%% AUX FUNCTIONS 

#CHECK_VALID_CHROMOSOME
# Description: This function checks if the given individual is valid. Note that
#   an individual is not valid if no trips are offered by any of transportists 
#   for a route.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the prices.
#   @Outputs:
#       Boolean: True if valid and false otherwise.
def check_valid_chromosome(chrom, mast_np):
    np_reduced = mast_np[:,chrom] #Get valid values
    return [True]*len(chrom) not in np.isnan(np_reduced).tolist()


#CALCULATE_MIN_ROW
# Description: This function calculates the minimum selection (one item per 
#   row) for a given individual.
#
#   @Inputs:
#       chrom: individual.
#       mast_np: numpy array with the prices.
#   @Outputs:
#       Minimum price.  
def calculate_min_row(chrom, mast_np):
    np_reduced = mast_np[:,chrom] #Get valid values
    price = 0
    for i in range(len(np_reduced)):
        x = np_reduced[i,:]
        x = x[~np.isnan(x)]
        price += min(x)
    return price


#CALCULATE_MAX_ROW
# Description: This function calculates the maximum selection of items (one 
#   item per row) for a given individual.
#
#   @Inputs:
#       chrom: individual.
#       mast_np: numpy array with the prices.
#   @Outputs:
#       Minimum price.  
def calculate_max_row(chrom, mast_np):
    np_reduced = mast_np[:,chrom] #Get valid values
    price = 0
    for i in range(len(np_reduced)):
        x = np_reduced[i,:]
        x = x[~np.isnan(x)]
        price += max(x)
    return price


#GET_MASTER_NP
# Description: This is a VERY auxiliar function. It receives either a pandas df
#   or a numpy df with the input data of the entities, and it transforms it to 
#   the structure needed to perform the analysis.
#
#   @Inputs:
#       df: pandas or numpy array.
#   @Outputs:
#       mast_np: Numpy array with the data.
#       colnames: Column names (names of the entities). If df is a numpy array
#           then these column names are taken as their indices.
def get_master_np(df):
    if type(df) == pd.core.frame.DataFrame: #Pandas 
        colnames = df.columns
        mast_np = df.to_numpy()
        return mast_np, colnames
    elif type(df) == np.ndarray: #Numpy
        colnames = list(map(lambda x: str(x),list(range(df.shape[1]))))
        return df, colnames
    else:
        raise ValueError('The dataframe (df) has to be either a numpy array or a pandas dataframe')


#%% FUNCTIONS GA

######################
# CREATE POPULATION: #
######################
    
#CREATE_INDIVIDUAL
# Description: This function creates individuals with a random and variable 
#   length of the genes between the specified limits. Each individual is a set 
#   of ints.
#
#   @Inputs:
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes
#       mast_np: numpy array with the prices.
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
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes
#   @Outputs:
#       list_individuals: List with the individuals.
def population(num_individuals, mast_np, min_number_of_genes, max_number_of_genes):
    return [create_individual(min_number_of_genes, max_number_of_genes, mast_np) for i in range(num_individuals)]

   
######################
# CALCULATE FITNESS: #
######################

#FITNESS_CALCULATION
# Description: This function calculates the fitness of an individual.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the prices.
#       min_item_per_row: If 1, it is calculated the minimum selection of items 
#       (one item per row) for a given individual. If it is 0 the maximum is 
#       calculated.
#       MAX_NUM_TRANS: Maximum number of transportists up to which there is no 
#           penalization.
#       PENALIZATION: Penalization value per each new transportist (or for 
#           their rating).
#       PERCENT: Percentage of the total cost that penalizes each extra new 
#           trasnportists.
#       RATING_TRANS: List with the ordered rating of each transportist.
#   @Outputs:
#       fitness_val: Value of the fitness  
def fitness_calculation(chrom, mast_np, min_item_per_row = 1, MAX_NUM_TRANS = 3, PENALIZATION = 0, PERCENT = 0, RATING_TRANS = []):
    if min_item_per_row:
        price = calculate_min_row(chrom, mast_np)
    else:
        price = calculate_max_row(chrom, mast_np)
    pen_num_trans_squ = PENALIZATION*((max(MAX_NUM_TRANS,len(chrom))-MAX_NUM_TRANS)**2) #Penalization for the number of transportists, squared
#    pen_num_trans_lin = PENALIZATION*(max(MAX_NUM_TRANS,len(chrom))-MAX_NUM_TRANS) #Penalization for the number of transportists, linear
#    pen_percent = price*(max(MAX_NUM_TRANS,len(chrom))-MAX_NUM_TRANS)*PERCENT #Penalization as a percentage of the total cost for each extra transportist
#    pen_rating = PENALIZATION*sum(RATING_TRANS[chrom])/len(chrom) #Penalization for the rating of the transportists
    fitness_val = price + pen_num_trans_squ #+ pen_num_trans_lin + pen_percent + pen_rating
    return fitness_val
    

#CALCULATE_FITNESS_AND_ORDER
# Description: This function calculates the fitness of the Individuals, it 
#   normalizes it, and then it order the population from best to worst fitness.
#   Note that it order the input value pop, and it doesn´t return anything.
#
#   @Inputs:
#       pop: population (list with Individuals).
#       mast_np: numpy array with the prices.
#       min_item_per_row: If 1, it is calculated the minimum selection of items 
#       (one item per row) for a given individual. If it is 0 the maximum is 
#       calculated.
#       minimize: If 1 the fitness value is inverse normalized, i.e, higher
#           values mappped to 0 and lower to 1. Otherwise it is normaly
#           normalized.
#       MAX_NUM_TRANS: Maximum number of transportists up to which there is no 
#           penalization.
#       PENALIZATION: Penalization value per each new transportist (or for 
#           their rating).
#       PERCENT: Percentage of the total cost that penalizes each extra new 
#           trasnportists.
#       RATING_TRANS: List with the ordered rating of each transportist.
#   @Outputs:
#       None
def calculate_fitness_and_order(pop, mast_np, min_item_per_row = 1, minimize = 1, MAX_NUM_TRANS = 3, PENALIZATION = 0, PERCENT = 0, RATING_TRANS = []):
    for ind in pop:
        chrom = ind.chromosome
        ind.fitness = fitness_calculation(chrom, mast_np, min_item_per_row, MAX_NUM_TRANS, PENALIZATION, PERCENT, RATING_TRANS)
    max_v = max(list(map(lambda x: x.fitness, pop)))
    min_v = min(list(map(lambda x: x.fitness, pop)))
    if minimize:
        for ind in pop:
            ind.normalized_fitness = 1-(ind.fitness-min_v)/(max_v-min_v)
    else:
        for ind in pop:
            ind.normalized_fitness = (ind.fitness-min_v)/(max_v-min_v)       
    pop.sort(key=lambda x: x.normalized_fitness, reverse=True) #Sort by normalized fitness
#pop = population(100, mast_np, 3, 5)
#calculate_fitness_and_order(pop, mast_np, minimize = 1, MAX_NUM_TRANS = 3, PENALIZATION = 0, PERCENT = 0, RATING_TRANS = [])

######################
# PERFORM SELECTION: #
######################

#ROULETTE_SELECTION
# Description: This function returns a list with all the individuals selected 
#   by the roulette wheel selection.
#    
#   @Inputs:
#       ordered_pop: output of the function calculate_fitness_and_order.
#       num_selected_ind: number of selected individuals.
#   @Outputs:
#       list_selected_ind: List with the IDs of the selected individuals.
def roulette_selection(ordered_pop, num_selected_ind):
    sumFit = sum(map(lambda x: x.normalized_fitness, ordered_pop)) #Sum of the normalized fitness
    list_selected_ind = []
    for i in range(num_selected_ind):
        cumSumFit = 0 #Cumulative sum
        rndFit = rnd()*sumFit #Sum fitness of the selected individual
        for i in range(len(ordered_pop)):
            cumSumFit += ordered_pop[i].normalized_fitness
            if rndFit <= cumSumFit:
                list_selected_ind.append(i)
                break
    return list_selected_ind

####################
# PERFORM PAIRING: #
####################
    
#PAIRING
# Description: This function performs random pairing.
#    
#   @Inputs:
#       list_selected_ind: List with the IDs of the selected individuals. It is
#           the output of the function roulette_selection.
#   @Outputs:
#       paired_ind: List of tuples with the IDs of the mated of the individuals.
def pairing(list_selected_ind):
    list_sel = list_selected_ind.copy()
    if len(list_sel) % 2 == 1:
        list_sel.pop() #If odd, pop last individual of the list
    shuffle(list_sel)
    paired_ind = []
    while len(list_sel) > 0:
        paired_ind.append((list_sel.pop(), list_sel.pop()))
    return paired_ind

######################
# PERFORM CROSSOVER: #
######################
    
#CALC_POSSIBLE_CROSSOVERS
# Description: This function calculate the possible crossovers that give a 
#   valid chromosome between two individuals. There are four possible types of 
#   crossover, taking 1 element from individual A and 1 from individual B; 
#   taking 2 elements from individual A and 1 from individual B; taking 1 
#   element from individual A and 2 from individual B; taking 2 elements from
#   individual A and 2 from individual B.
#
#   @Inputs:
#       chrom_a: Individual A's chromosome.
#       chrom_b: Individual B's chromosome.
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#   @Outputs:
#       pos_crossovers: List of tuples in which the first element has the 
#           elements of ind_a that give a possible crossover and the second the
#           elements of ind_b. Ie, the first element of the tuple belongs to 
#           ind_a and should be put in ind_b, and viceversa.
def calc_possible_crossovers(chrom_a, chrom_b, mast_np, min_number_of_genes, max_number_of_genes):
    pairs_a = list(combinations(chrom_a, 2)) #Get all the possible combinations of two elements in ind_a
    pairs_b = list(combinations(chrom_b, 2)) #Get all the possible combinations of two elements in ind_b
    pos_crossovers = []
    # 1 element from individual A and 1 from individual B
    for el_a in chrom_a:
        for el_b in chrom_b:
            if (el_a not in chrom_b) and (el_b not in chrom_a): #Check that the element is not repeated
                crossed_1 = [el_b if x==el_a else x for x in chrom_a]
                crossed_2 = [el_a if x==el_b else x for x in chrom_b]
                if check_valid_chromosome(crossed_1, mast_np) and check_valid_chromosome(crossed_2, mast_np):
                    pos_crossovers.append(([el_a], [el_b])) #If valid chromosomes, append in pos_crossovers
    # 2 elements from individual A and 1 from individual B
    if (len(el_a) > min_number_of_genes) and (len(el_b) < max_number_of_genes) and (len(el_a) >= 2): #If the length of both, a and b are not in the limits of length
        for el_a in pairs_a:
            el_a_1 = el_a[0]
            el_a_2 = el_a[1]
            for el_b in chrom_b:
                if (el_a_1 not in chrom_b) and (el_a_2 not in chrom_b) and (el_b not in chrom_a): #Check that the element is not repeated       
                    crossed_1 = chrom_a.copy()
                    crossed_2 = chrom_b.copy()
                    #Eliminate the two elements of chrom_a and append the new one
                    crossed_1.remove(el_a_1) 
                    crossed_1.remove(el_a_2)
                    crossed_1.append(el_b)
                    #Eliminate the element of chrom_b and append the new ones
                    crossed_2.remove(el_b)
                    crossed_2.append(el_a_1)
                    crossed_2.append(el_a_2)
                    if check_valid_chromosome(crossed_1, mast_np) and check_valid_chromosome(crossed_2, mast_np):
                        pos_crossovers.append(([el_a_1, el_a_2], [el_b])) #If valid chromosomes, append in pos_crossovers
    # 1 element from individual A and 2 from individual B
    if (len(el_b) > min_number_of_genes) and (len(el_a) < max_number_of_genes) and (len(el_b) >= 2): #If the length of both, a and b are not in the limits of length
        for el_b in pairs_b:
            el_b_1 = el_b[0]
            el_b_2 = el_b[1]
            for el_a in chrom_a:
                if (el_b_1 not in chrom_a) and (el_b_2 not in chrom_a) and (el_a not in chrom_b) and (len(chrom_b) > 2): #Check that the element is not repeated and the length of chrom_b is bigger than 2
                    crossed_1 = chrom_a.copy()
                    crossed_2 = chrom_b.copy()
                    #Eliminate the element of chrom_a and append the new ones              
                    crossed_1.remove(el_a) 
                    crossed_1.append(el_b_1)
                    crossed_1.append(el_b_2)
                    #Eliminate the two elements of chrom_b and append the new one
                    crossed_2.remove(el_b_1)
                    crossed_2.remove(el_b_2)
                    crossed_2.append(el_a)
                    if check_valid_chromosome(crossed_1, mast_np) and check_valid_chromosome(crossed_2, mast_np):
                        pos_crossovers.append(([el_a], [el_b_1, el_b_2])) #If valid chromosomes, append in pos_crossovers
    # 2 element from individual A and 2 from individual B
    if (len(el_a) >= 2) and (len(el_b) >= 2): #If the length of both, a and b are not in the limits of length
        for el_a in pairs_a:
            el_a_1 = el_a[0]
            el_a_2 = el_a[1]
            for el_b in pairs_b:
                el_b_1 = el_b[0]
                el_b_2 = el_b[1]
                if (el_b_1 not in chrom_a) and (el_b_2 not in chrom_a) and (el_a_1 not in chrom_b) and (el_a_2 not in chrom_b): #Check that the element is not repeated
                    crossed_1 = chrom_a.copy()
                    crossed_2 = chrom_b.copy()
                    #Eliminate the elements of chrom_a and append the new ones              
                    crossed_1.remove(el_a_1) 
                    crossed_1.remove(el_a_2)
                    crossed_1.append(el_b_1)
                    crossed_1.append(el_b_2)
                    #Eliminate the elements of chrom_b and append the new ones
                    crossed_2.remove(el_b_1)
                    crossed_2.remove(el_b_2)
                    crossed_2.append(el_a_1)
                    crossed_2.append(el_a_2)
                    if check_valid_chromosome(crossed_1, mast_np) and check_valid_chromosome(crossed_2, mast_np):
                        pos_crossovers.append(([el_a_1, el_a_2], [el_b_1, el_b_2])) #If valid chromosomes, append in pos_crossovers
    return pos_crossovers


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
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#   @Outputs:
#       list_new_ind: List with the new chromosomes. 
def mating(paired_ind, ordered_pop, mast_np, min_number_of_genes, max_number_of_genes):
    list_new_ind = []
    for (ind_a, ind_b) in paired_ind:
        chrom_a = ordered_pop[ind_a].chromosome.copy()
        chrom_b = ordered_pop[ind_b].chromosome.copy()
        pos_crossovers = calc_possible_crossovers(chrom_a, chrom_b, mast_np, min_number_of_genes, max_number_of_genes)
        if len(pos_crossovers) == 0: #There are no possible crossovers
            list_new_ind.append(ordered_pop[ind_a].chromosome.copy())
            list_new_ind.append(ordered_pop[ind_b].chromosome.copy())
        else: #Else, choose a random crossover, perform it, and append the new individuals to list_new_ind
            cross = choice(pos_crossovers) #Get a random crossover from the list of possible crossovers
            elim_a = cross[0]
            elim_b = cross[1]
            new_chrom_a = [e for e in chrom_a if e not in elim_a]
            new_chrom_a.extend(elim_b)
            new_chrom_b = [e for e in chrom_b if e not in elim_b]
            new_chrom_b.extend(elim_a)
            list_new_ind.append(new_chrom_a) #New individuals with temporary fitness value of 0
            list_new_ind.append(new_chrom_b)         
    return list_new_ind

#####################
# PERFORM MUTATION: #
#####################
    
#MUTATION
# Description: This function receives a chromosome and performs a random 
#   mutation over one of its elements. It iterates over all the possible 
#   mutations until one is found, moment in which the execution stops. If no 
#   mutation is found, then it returns theinput chromosome.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the prices.
#       mutation_type: It can be either:
#           'mut_gen': One of the genes of an individual is randomly changed.
#           'addsub_gen': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed: It is the number of genes that is changed, ie, the 
#           number of genes that are added, substracted or changed by new ones.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#   @Outputs:
#       new_chrom: mutated chromosome
def mutation(chrom, mast_np, mutation_type = 'both', num_gen_changed = 1, min_number_of_genes = None, max_number_of_genes = None):
    num_elem_to_choose = mast_np.shape[1] #Number of elements from which to choose
    list_elem = range(num_elem_to_choose)
    list_elem = [e for e in list_elem if e not in chrom] #Get the elements that can be chosen as new part of the chromosome (new elements).
    shuffle(list_elem) #randomnize
    both_tp = int(rnd() > 0.5) #If mutation_type = 'both' it is chosen randomly one of both mutation types.
    new_chrom = chrom.copy()
    shuffle(new_chrom) #randomnize
    FLAG_MUT = 0 #Number of genes changed.
    if (mutation_type == 'mut_gen') or ((mutation_type == 'both') and (both_tp == 0)):
        for el in chrom:
            for i in list_elem:
                aux_chrom = [i if x==el else x for x in new_chrom]
                if check_valid_chromosome(aux_chrom, mast_np): #Break when the first valid mutation is found
                    FLAG_MUT += 1
                    list_elem.remove(i)
                    new_chrom = aux_chrom
                    break
            if FLAG_MUT >= num_gen_changed:
                break
    if (mutation_type == 'addsub_gen') or ((mutation_type == 'both') and (both_tp == 1)):
        if len(chrom) > max_number_of_genes-num_gen_changed:
            add = 0 #sub
        elif len(chrom) < min_number_of_genes+num_gen_changed:
            add = 1 #add
        else: 
            add = int(rnd() > 0.5) #Randomly choose to add or sub a gen
        if add: #Add a new item
            for i in range(num_gen_changed):
                it = list_elem.pop()
                new_chrom.append(it)
        if not add: #Eliminate an item
            aux_chrom = new_chrom.copy()
            for i in range(num_gen_changed):
                for el in aux_chrom:
                    aux_chrom_test = aux_chrom.copy()
                    aux_chrom_test.remove(el)
                    if check_valid_chromosome(aux_chrom_test, mast_np):
                        aux_chrom = aux_chrom_test.copy()
                        break
                new_chrom = aux_chrom.copy()
    return new_chrom

########################
# GET NEXT GENERATION: #
########################

#NEXT_GENERATION
# Description: This function returns the chromosomes of the next generation, 
#   specified the elitism, and mutation rate.
#    
#   @Inputs:
#       ordered_pop: output of the function calculate_fitness_and_order.
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#       elitism_rate: percentage of elitism. 
#       mutation_rate: Mutation rate per individual.
#       mutation_type: It can be either:
#           'mut_gen': One of the genes of an individual is randomly changed.
#           'addsub_gen': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed: It is the number of genes that is changed, ie, the 
#           number of genes that are added, substracted or changed by new ones.
#   @Outputs:
#       new_gen: List with the chromosomes of the new generation.
def next_generation(ordered_pop, mast_np, min_number_of_genes, max_number_of_genes, elitism_rate = 0.05, mutation_rate = 0.4, mutation_type = 'both', num_gen_changed = 1):
    new_gen = [] #List with the chromosomes of the next generation
    num_elitism = int(len(ordered_pop)*elitism_rate) #Number of individuals chosen as elite
    num_cross = len(ordered_pop) - num_elitism #Number of individuals chosen for crossover
    if (num_cross % 2) == 1: #Odd number of Individuals for crossover
        num_elitism += 1 #Add one individual to elitism
        num_cross -= 1 #Remove one from crossover
    #Elite:
    new_gen.extend(x.chromosome for x in ordered_pop[:num_elitism]) #Add elite
    #Crosover:
    selected_ind = roulette_selection(ordered_pop, num_cross) #1. Roulette selection
    paired_ind = pairing(selected_ind) #2. Perform the pairing
    new_cross_ind = mating(paired_ind, ordered_pop, mast_np, min_number_of_genes, max_number_of_genes) #3. Get the new chromosomes already crossovered
    new_gen.extend(new_cross_ind) #4. Add crossovered individuals
    #Mutation:
    num_mutated = int(len(ordered_pop)*mutation_rate) #Number of individuals chosen to be mutated
    list_IDs = list(range(num_elitism, len(mast_np))) #All the IDs of the individuals that may be mutated (not the elite)
    shuffle(list_IDs)
    mut_IDs = list_IDs[:num_mutated] #Get list of IDs of chromosomes to be mutated
    for i in mut_IDs:
        new_gen[i] = mutation(new_gen[i], mast_np, mutation_type, num_gen_changed, min_number_of_genes, max_number_of_genes)
    return new_gen

#########################
# TERMINATION CRITERIA: #
#########################
    
#CHECK_TERMINATION_CRITERIA
# Description: This checks the termination criteria
#    
#   @Inputs:
#       termination_criteria: Termination criteria's function to call.
#       gen_num: (max_num_generation_reached) Number of the current generation.
#       max_gen: (max_num_generation_reached) Maximum number of generations.
#   @Outputs:
#       Boolean. True if the maximum number of generations is reached
def check_termination_criteria(termination_criteria, gen_num = None, max_gen = None):
    if termination_criteria == 'max_num_generation_reached':
        return max_num_generation_reached(gen_num, max_gen)
    

#MAX_NUM_GENERATION_REACHED
# Description: This function returns True if the maximum number of generations 
#   is reached.
#    
#   @Inputs:
#       gen_num: Number of the current generation.
#       max_gen: Maximum number of generations.
#   @Outputs:
#       Boolean. True if the maximum number of generations is reached
def max_num_generation_reached(gen_num, max_gen):
    return (gen_num >= max_gen)

################################
# GENETIC ALGORITHM EXECUTION: #
################################

#GENETIC_ALGORITHM
# Description: This function is the main function for executing the GA and the
#   only one that should be execution.
#
#   @Inputs:
#       num_individuals: Number of individuals.
#       df: numpy/pandas array in which each column represents an entity to be 
#           chosen (eg. an enterprise) and each row an item in which that  
#           entity participates. If it is a pandas df, then the colnames are   
#           taken as the entities (for presentation of the results purposes).
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#       termination_criteria: Termination criteria's function to call.
#       elitism_rate: percentage of elitism. 
#       mutation_rate: Mutation rate per individual.
#       mutation_type: It can be either:
#           'mut_gen': One of the genes of an individual is randomly changed.
#           'addsub_gen': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed: It is the number of genes that is changed, ie, the 
#           number of genes that are added, substracted or changed by new ones, 
#           when an individual is selected for mutation.
#       max_gen: (max_num_generation_reached) Maximum number of generations.
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
#       best_individual: Best Individual (of the class individual).
#       best_fit_per_gen: Array with the best individual per generation.
#       pop: Whole population.
def GA_oipr(num_individuals, df, min_number_of_genes, max_number_of_genes, termination_criteria  = 'max_num_generation_reached', elitism_rate = 0.1, mutation_rate = 0.2, mutation_type = 'both', num_gen_changed = 1, max_gen = 100, min_item_per_row = 1, minimize = 1, MAX_NUM_TRANS = 3, PENALIZATION = 0, PERCENT = 0, RATING_TRANS = None):
    mast_np, colnames = get_master_np(df)
    if num_individuals*elitism_rate < 1: #Error if the number of individuals is not big enough...
        raise ValueError('With this elitism rate, it is needed, at least, ' + str(int(1/elitism_rate)) + ' individuals per generation')
    generation_num  = 0 #Initialize the generation counter to 0
    best_fit_per_gen = []
    pop = population(num_individuals, mast_np, min_number_of_genes, max_number_of_genes) #Creation of the initial population
    calculate_fitness_and_order(pop, mast_np, min_item_per_row, minimize, MAX_NUM_TRANS, PENALIZATION, PERCENT, RATING_TRANS) #Calculate and order the initial population by fitness
    best_fit_per_gen.append(pop[0].fitness)
    if check_termination_criteria(termination_criteria, generation_num, max_gen):
        return pop[0], best_fit_per_gen, pop
    while(not check_termination_criteria(termination_criteria, generation_num, max_gen)): #Run GA while the termination criteria is not met
        generation_num += 1 # Add 1 to the generation number
        print('Generation: ' + str(generation_num))
        new_chromosomes = next_generation(pop, mast_np, min_number_of_genes, max_number_of_genes, elitism_rate, mutation_rate, mutation_type, num_gen_changed) #Get the new generation chromosomes
        for i in range(len(pop)):
            pop[i].chromosome = new_chromosomes[i] #Copy new generation
        calculate_fitness_and_order(pop, mast_np, minimize, MAX_NUM_TRANS, PENALIZATION, PERCENT, RATING_TRANS) #Calculate and order the population by fitness
        best_fit_per_gen.append(pop[0].fitness)
    return pop[0], best_fit_per_gen, pop

#%%
#a, b, c = GA_oipr(50, mast_np, min_number_of_genes = 3, max_number_of_genes = 5, PENALIZATION = 12111)

#%% TODO
    
# TODO_0: Hacer algo para mantener la diversidad -> Convergencia muy muy muy rápida
# TODO_1: Controlar con algún heurístico que no haya ninguna ruta con todo nans (por ejemplo que haya al menos 3 transportistas por ruta o si no quitarla).
# TODO_2: Crear una población inicial con individuos únicos (sin repeticiones).
# TODO_3: Test other selection/pairing/mating methods
# TODO_4: Revisar IMPORTANCIA_DIFF_TRANS
# TODO_5: Fitness de MINLP y de GA no coinciden...




#CALC_POSSIBLE_CROSSOVERS
# Description: This function calculate the possible crossovers that give a 
#   valid chromosome between two individuals. There are four possible types of 
#   crossover, taking 1 element from individual A and 1 from individual B; 
#   taking 2 elements from individual A and 1 from individual B; taking 1 
#   element from individual A and 2 from individual B; taking 2 elements from
#   individual A and 2 from individual B.
#
#   @Inputs:
#       chrom_a: Individual A's chromosome.
#       chrom_b: Individual B's chromosome.
#       mast_np: numpy array with the prices.
#       min_number_of_genes: Minimum number of genes.
#       max_number_of_genes: Maximum number of genes.
#       max_num_gen_changed: Maximum number of genes that can be taken 
#           of an individual in order to make the crossover.
#   @Outputs:
#       pos_crossovers: List of tuples in which the first element has the 
#           elements of ind_a that give a possible crossover and the second the
#           elements of ind_b. Ie, the first element of the tuple belongs to 
#           ind_a and should be put in ind_b, and viceversa.
def calc_crossover(chrom_a, chrom_b, mast_np, min_number_of_genes, max_number_of_genes, max_num_gen_changed):
    num_gen_a_interchanged = list(map(lambda x: x+1, list(range(max_num_gen_changed))))
    shuffle(num_gen_a_interchanged) #random number of genes of a that are interchanged
    num_gen_b_interchanged = list(map(lambda x: x+1, list(range(max_num_gen_changed)))) 
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
  





