"""
In this file it is defined the class Individual.

Classes:
    :Individual: Main class.
"""
from inspect import signature
import uuid


class Individual:
    """ Class  of the individuals. It has the attributes chromosome, fitness and normalized fitness.
    """

    def __init__(self, chromosome):
        """ Constructor.

        :param chromosome: (list of genes) Chromosome of the individual.
        """
        if type(chromosome) != list:
            raise AttributeError('The chromosome must be a list of genes')
        else:
            self.chromosome = chromosome
            self._id = str(
                uuid.uuid4())  # Unique ID for each individual ---> It is unique because uuid uses the time component for creating this id
            self.fitness_value = None  # This attribute will be filled when the individual is evaluated
            self.normalized_fitness_value = None  # This attribute holds the normalized value of the fitness in comparison to the rest of the population
            self.inverse_normalized_fitness_value = None  # This attribute holds the inverse normalized value of the fitness in comparison to the rest of the population

    def set_new_chromosome(self, chromosome):
        """ Method to set a chromosome.

        :param chromosome: (list of genes) Chromosome of the individual.
        :return:
        """
        if type(chromosome) != list:
            raise AttributeError('The chromosome must be a list of genes')
        else:
            self.chromosome = chromosome

    def set_fitness_value(self, fitness_value):
        """ Method to set the fitness value.

        :param fitness_value: (float) Fitness value.
        :return:
        """
        if type(fitness_value) == int or type(fitness_value) == float:
            self.fitness_value = fitness_value
        else:
            raise ValueError('The fitness value must be a float.')

    def set_normalized_fitness_value(self, normalized_fitness_value):
        """ Method to set the normalized fitness value.

        :param normalized_fitness_value: (float) Fitness value.
        :return:
        """
        if type(normalized_fitness_value) == int or type(normalized_fitness_value) == float:
            self.normalized_fitness_value = normalized_fitness_value
        else:
            raise ValueError('The normalized fitness value must be a float.')

    def set_inverse_normalized_fitness_value(self, inverse_normalized_fitness_value):
        """ Method to set the inverse normalized fitness value.

        :param inverse_normalized_fitness_value: (float) Fitness value.
        :return:
        """
        if type(inverse_normalized_fitness_value) == int or type(inverse_normalized_fitness_value) == float:
            self.inverse_normalized_fitness_value = inverse_normalized_fitness_value
        else:
            raise ValueError('The inverse normalized fitness value must be a float.')

    def calculate_fitness(self, fitness):
        """ Method to calculate the fitness of the individual.

        :param fitness: (function) Function to evaluate the fitness. Its ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness.
        :return:
            * :fitness_value: (float) The fitness value
            * It is also set the value of the attribute Individual.fitness_value to the calculated value when this method is called.
        """
        try:
            value = fitness(self.chromosome)
            self.set_fitness_value(value)
            return value
        except Exception as e:
            raise (type(e))(str(e) + '\n' +
                            "Error when calculating the fitness of an individual")

    def kill_and_reset(self, chromosome):
        """ This is a bit of a "tricky" method. If a new individual of a new generation is going to be created, it will be computationally much faster to get the object of an individual that is going to be killed and reset its values to the ones of the new individual. ---> This will be done instead of directly creating a completely new individual with the new data and stop referencing the old individual.

        :param chromosome: (list of genes) Chromosome of the individual.
        """
        if type(chromosome) != list:
            raise AttributeError('The chromosome must be a list of genes')
        else:
            self.chromosome = chromosome
            self.fitness_value = None  # This attribute will be filled when the individual is evaluated
            self.normalized_fitness_value = None  # This attribute holds the normalized value of the fitness in comparison to the rest of the population
            self.inverse_normalized_fitness_value = None  # This attribute holds the inverse normalized value of the fitness in comparison to the rest of the population

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.chromosome == other.chromosome

    def __lt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.fitness < other.fitness

    def __gt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.fitness > other.fitness

    def __repr__(self):
        return '{self.__class__.__name__}({self.chromosome})'.format(self=self)

    def __str__(self):
        return "Chromosome: {0} \nFitness: {1}".format(self.chromosome, self.fitness_value)

    def __iter__(self):
        return iter(self.chromosome)

    def __len__(self):
        return len(self.chromosome)

    def __getitem__(self, pos):
        return self.chromosome[pos]
