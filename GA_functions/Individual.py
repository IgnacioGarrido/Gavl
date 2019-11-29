#!/usr/bin/env python3
"""
In this file it is defined the class Individual.
"""

# CLASS INDIVIDUAL:
#
# This is the class of the individuals. It has the attributes chromosome,
#   fitness and normalized fitness. 
class Individual:
  def __init__(self, chromosome, fitness, normalized_fitness):
    self.chromosome = chromosome
    self.fitness = fitness
    self.normalized_fitness = normalized_fitness