#!/usr/bin/env python3
"""
In this file it is defined the function to perform mutation.

Function: 
    mutation: function that performs mutation.

IF A DIFFERENT MUTATION IS NEEDED, IT SHOULD BE PROGRAMMED IN THE FUNCTION mutation
"""
from .aux_functions import check_valid_chromosome
from random import random as rnd
from random import shuffle
from itertools import combinations

#MUTATION
# Description: This function receives a chromosome and performs a random 
#   mutation over one of its elements. It iterates over all the possible 
#   mutations until one is found, moment in which the execution stops. If no 
#   mutation is found, then it returns the input chromosome.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the weights of each item.
#       mutation_type: It can be either:
#           'mut_gene': One of the genes of an individual is randomly changed.
#           'addsub_gene': It is added or eliminated (randomly) a new element to
#               an Individual.
#           'both': 'mut_gen' and 'addsub_gen' are randomly applied.
#       num_gen_changed_mutation: It is the number of genes that is changed,
#           ie, the number of genes that are added, substracted or changed by 
#           new ones.
#       min_number_of_genes: Minimum number of genes of the chromosome.
#       max_number_of_genes: Maximum number of genes of the chromosome.
#   @Outputs:
#       new_chrom: mutated chromosome
def mutation(chrom, mast_np, mutation_type = 'both', num_gen_changed_mutation = 1, min_number_of_genes = None, max_number_of_genes = None):
    num_elem_to_choose = mast_np.shape[1] #Number of elements from which to choose
    list_elem = range(num_elem_to_choose)
    list_elem = [e for e in list_elem if e not in chrom] #Get the elements that can be chosen as new part of the chromosome (new elements).
    shuffle(list_elem) #randomnize
    both_tp = int(rnd() > 0.5) #If mutation_type = 'both' it is chosen randomly one of both mutation types.
    new_chrom = chrom.copy()
    shuffle(new_chrom) #randomnize
    FLAG_MUT = 0 #Number of genes changed.
    if (mutation_type == 'mut_gene') or ((mutation_type == 'both') and (both_tp == 0)):
        for el in chrom:
            for i in list_elem:
                aux_chrom = [i if x==el else x for x in new_chrom]
                if check_valid_chromosome(aux_chrom, mast_np): #Break when the first valid mutation is found
                    FLAG_MUT += 1
                    list_elem.remove(i)
                    new_chrom = aux_chrom
                    break
            if FLAG_MUT >= num_gen_changed_mutation:
                break
    if (mutation_type == 'addsub_gene') or ((mutation_type == 'both') and (both_tp == 1)):
        if len(chrom) > max_number_of_genes-num_gen_changed_mutation:
            add = 0 #sub
        elif len(chrom) < min_number_of_genes+num_gen_changed_mutation:
            add = 1 #add
        else: 
            add = int(rnd() > 0.5) #Randomly choose to add or sub a gen
        if add: #Add a new item
            if len(list_elem) >= num_gen_changed_mutation:
                add_combs = list(combinations(list_elem, num_gen_changed_mutation))
            else:
                add_combs = list_elem
            if len(list_elem) == 0: #If there are no possibilities to mutate (all the elements are being used)
                return chrom
            aux_chrom = new_chrom.copy()
            for comb in add_combs:
                for el in comb:
                    aux_chrom.append(el)
                if check_valid_chromosome(aux_chrom, mast_np):
                    new_chrom = aux_chrom.copy() 
                    break
        if not add: #Eliminate an item
            aux_chrom = new_chrom.copy()
            for i in range(num_gen_changed_mutation):
                for el in aux_chrom:
                    aux_chrom_test = aux_chrom.copy()
                    if len(aux_chrom) > 1: #Check that there are at least one element
                        aux_chrom_test.remove(el)
                    if check_valid_chromosome(aux_chrom_test, mast_np):
                        aux_chrom = aux_chrom_test.copy()
                        break
                new_chrom = aux_chrom.copy()
    return new_chrom