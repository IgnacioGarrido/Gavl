#!/usr/bin/env python3
"""
In this file it is defined auxiliar functions that are useful for the main 
processes.

Functions:
    check_valid_chromosome: Returns true if the chromosome is valid.
    calculate_min_row: Given a chromosome, it calculates the minimum selection 
        of items (one per row).
    calculate_max_row: Given a chromosome, it calculates the maximum selection 
        of items (one per row).
    get_master_np: Auxiliar function to get the master array with all the 
        weights of the items, given a pandas df.
        
IF A NEW DEFINITION OF HOW IT IS A VALID CHROMOSOME IS NEEDED, IT SHOULD BE CHANGED THE FUNCTION check_valid_chromosome
"""
import numpy as np
import pandas as pd

#CHECK_VALID_CHROMOSOME
# Description: This function checks if the given individual is valid. Note that
#   an individual is not valid if no trips are offered by any of transportists 
#   for a route.
#
#   @Inputs:
#       chrom: chromosome.
#       mast_np: numpy array with the weights of each item.
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
#       mast_np: numpy array with the weights of each item.
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
#       mast_np: numpy array with the weights of each item.
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
#       mast_np: numpy array with the weights of each item.
#       colnames: Column names (names of the entities). If df is a numpy array
#           then these column names are taken as their indices.
def get_master_np(df):
    if type(df) == pd.core.frame.DataFrame: #Pandas 
        colnames = list(df.columns)
        mast_np = df.to_numpy()
        return mast_np, colnames
    elif type(df) == np.ndarray: #Numpy
        #colnames = list(map(lambda x: str(x),list(range(df.shape[1]))))
        colnames = list(range(df.shape[1]))
        return df, colnames
    else:
        raise ValueError('The passed argument to df has to be either a numpy array or a pandas dataframe')
  
#INDEX_TO_COLNAME
# Description: This is a VERY auxiliar function. It receives the population, 
#   being its chromosome  the values of the columns that are chosen and it 
#   changes it to its real names.
#
#   @Inputs:
#       pop: Population - output of the function calculate_fitness_and_order.
#       colnames: Column names (names of the entities). If df is a numpy array
#           then these column names are taken as their indices.
#   @Outputs:
#       None  
def index_to_colname(pop, colnames):
    for ind in pop:
        ind.chromosome = [colnames[i] for i in ind.chromosome]