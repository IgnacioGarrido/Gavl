#!/usr/bin/env python3
"""
In this file it is defined the functions to check the termination criteria.
    
Functions:
    max_num_generation_reached: Function that checks if the max number of 
        generations is reached.
    goal_fitness_reached: Function that checks if the goal fitness is reached.
    check_termination_criteria: this function selects the termination criteria
"""  

#MAX_NUM_GENERATION_REACHED
# Description: This function returns True if the maximum number of generations 
#   is reached.
#    
#   @Inputs:
#       gen_num: Number of the current generation.
#       max_gen: Maximum number of generations.
#   @Outputs:
#       Boolean. True if the maximum number of generations is reached.
def max_num_generation_reached(gen_num, max_gen):
    return (gen_num >= max_gen)


#GOAL_FITNESS_REACHED
# Description: This function returns True if the goal fitness is reached. 
#    
#   @Inputs:
#       gen_fitness: Current generation best fitness.
#       goal_fitness:  min_fitness_reached) Fitness goal.
#       minimize: If 1, the goal is minimize the fitness function and 
#           viceversa.
#   @Outputs:
#       Boolean. True if the goal fitness is reached.
def goal_fitness_reached(gen_fitness, goal_fitness, minimize):
    if minimize:
        return (gen_fitness <= goal_fitness)
    else:
        return (gen_fitness >= goal_fitness)


#CHECK_TERMINATION_CRITERIA
# Description: This checks the termination criteria
#    
#   @Inputs:
#       termination_criteria: Termination criteria's function to call.
#       gen_num: (max_num_generation_reached) Number of the current generation.
#       max_gen: (max_num_generation_reached) Maximum number of generations.
#       gen_fitness: (goal_fitness_reached) Current generation best fitness.
#       goal_fitness:  (goal_fitness_reached) Fitness goal.
#       minimize: (goal_fitness_reached) If 1, the goal is minimize the fitness function and 
#           viceversa.
#   @Outputs:
#       Boolean. True if the maximum number of generations is reached
def check_termination_criteria(termination_criteria, gen_num = None, max_gen = None, gen_fitness = None, goal_fitness = None, minimize = None):
    if termination_criteria == 'max_num_generation_reached':
        return max_num_generation_reached(gen_num, max_gen)
    elif termination_criteria == 'goal_fitness_reached':
        return goal_fitness_reached(gen_fitness, goal_fitness, minimize)