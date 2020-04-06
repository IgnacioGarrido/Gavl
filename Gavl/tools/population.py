"""
In this file it is defined the class to hold the population.

Classes:
    :Population: Main class.
"""
from .individual import Individual
from abc import ABC, abstractmethod


class Population(ABC):
    """ Class of the population.
    """

    def __init__(self):
        """Constructor.

        """
        # Algorithm's attributes:
        self.population = []  # When the algorithm starts it will be filled with the population. It is a list with all the individuals.

    def __set_population(self, population):
        """ Method to set the whole passed population of individuals.
        WARNING: This method deletes all the previous individuals of the population. If wanted to maintain the olf individuals use the method add_individual(individual).

        :param population: (list of chromosomes, Individuals or Population) List of individuals or object of the class Population.
        :return:
        """
        try:
            if type(population) == list and False not in [type(ind) == list for ind in
                                                          population]:  # If the passed object is a list of chromosomes
                self.population = []  # Reset the population
                for individual in population:
                    self.add_individual(individual)
            elif type(population) == list and False not in [type(ind) == Individual for ind in
                                                            population]:  # If the passed object is a list of Individuals
                self.population = population
            elif type(population) == Population and False not in [type(ind) == Individual for ind in population]:  # If it is passed an object of the class Population.
                self.population = population.population
            else:
                raise ValueError()
        except ValueError:
            raise ValueError(
                'The population must be a list of individuals of the class Individual or a list of chromosomes.')  # There has been an error
        except Exception as e:
            raise (type(e))(
                str(
                    e) + '\n' + 'The population must be a list of individuals of the class Individual or a Population object.')  # There has been an error

    def __kill_and_reset_whole_population(self, new_population):
        """ This is a bit of a "tricky" method. If a new population of individuals of a new generation is going to be created, it will be computationally much faster to get the objects of the actual individuals that are already in the population and that are going to be killed and reset their values to the ones of the new population. ---> This will be done instead of directly creating a completely population of new individuals with the new data and stop referencing the individuals of the old population.

        :param new_population: (list of chromosomes, Individuals or Population) List of individuals or object of the class Population.
        :return:
        """
        if not hasattr(new_population, '__len__'):
            ValueError(
                'The passed argument must be a list of chromosomes, Individuals (class Individual) or a whole Population.')
        elif len(self.population) != len(new_population):
            raise ValueError(
                'The new population must have the same number of individuals than the existing one. Note that the current population has {} individuals.'.format(
                    len(self.population)))
        else:
            try:
                if type(new_population) == list and False not in [type(ind) == list for ind in
                                                                  new_population]:  # If the passed object is a list of chromosomes
                    for i in range(len(self.population)):
                        self.population[i].kill_and_reset(new_population[i])
                elif type(new_population) == list and False not in [type(ind) == Individual for ind in
                                                                    new_population]:  # If the passed object is a list of Individuals
                    for i in range(len(self.population)):
                        self.population[i].kill_and_reset(new_population[i].chromosome)
                elif type(new_population) == Population and False not in [type(ind) == Individual for ind in new_population]:  # If it is passed an object of the class Population.
                    for i in range(len(self.population)):
                        self.population[i].kill_and_reset(new_population[i].chromosome)
                else:
                    raise ValueError()
            except ValueError:
                raise ValueError(
                    'The population must be a list of individuals of the class Individual or a list of chromosomes.')  # There has been an error
            except Exception as e:
                raise (type(e))(
                    str(
                        e) + '\n' + 'The population must be a list of individuals of the class Individual or a Population object.')  # There has been an error

    def add_individual(self, individual):
        """ Method to add a new individual to the population. BEWARE: The population may have a bigger size than the specified if new individuals are added.

        :param individual: (Individual or list) individual of the class Individual or list representing the chromosome.
        :return:
        """
        if type(individual) == list:
            new_ind = Individual(individual)
            self.population.append(new_ind)
        elif type(individual) == Individual:
            self.population.append(individual)
        else:
            raise ValueError(
                'The argument Individual must be of the class Individual or a list representing the chromosome.')

    def get_individual_by_id(self, id_individual):
        """ This method returns an Individual given its ID. If there is no Individual with this ID in the population, it returns None.

        :param id_individual: (str) ID of the individual.
        :return:
            * :individual: (Individual) Individual whose id is id_individual.
        """
        if type(id_individual) != str:
            raise ValueError('This method MUST receive the id of the individual, which is a ')
        else:
            for individual in self.population:
                if individual._id == id_individual:
                    return individual
            return None  # No individual with that id has been found

    def __calculate_normalized_fitness(self):
        """ This method calculates the normalized fitness of the whole population (and it adds the parameter to each individual (to the attribute normalized_fitness)). The same is done with the inverse normalize fitness.
        """
        if True in [ind.fitness_value is None for ind in self]:
            self.__calculate_fitness_population()  # Calculate the fitness of the individuals
        max_v = max([ind.fitness_value for ind in self])
        min_v = min([ind.fitness_value for ind in self])
        if max_v != min_v:  # There is no 0-division
            for ind in self.population:
                normalized_fitness_value = (ind.fitness_value - min_v) / (max_v - min_v)
                ind.set_inverse_normalized_fitness_value(1 - normalized_fitness_value)
                ind.set_normalized_fitness_value(normalized_fitness_value)
        else:  # There would be 0-division ---> In this generation the whole population has converged to the same fitness value.
            for ind in self.population:
                ind.set_inverse_normalized_fitness_value(1)
                ind.set_normalized_fitness_value(1)

    @abstractmethod
    def __calculate_fitness_population(self):
        """ This method calculates the fitness of all the individuals of the population and then it sets this attribute to each individual.
        """
        pass

    @abstractmethod
    def __generate_population(self):
        """ This method should generate a new population and add it to the attribute population.
        """
        pass

    @abstractmethod
    def __sort_population(self):
        """ This method orders the population by its fitness.
        """
        pass

    def __calculate_fitness_and_sort(self):
        """ This method is called to calculate the fitness of the whole population, then calculate the normalized fitness and then sort the population according to this fitness.
        """
        self.__calculate_fitness_population()
        self.__calculate_normalized_fitness()
        self.__sort_population()

    @abstractmethod
    def __get_next_generation(self):
        """ This method is used to calculate the next generation.

        :return:
            * :population: (list of individuals) It returns the next generation.
        """
        pass

    @abstractmethod
    def best_individual(self):
        """ This method returns the best individual.

        :return:
            * :individual: (Individual) Best individual.
        """
        pass

    def get_population(self):
        """ Method to get the population

        :return:
        """
        self.__calculate_fitness_population()
        self.__calculate_normalized_fitness()
        self.__sort_population()
        return self.population

    def __len__(self):
        return len(self.population)

    def __iter__(self):
        return iter(self.population)

    def __getitem__(self, pos):
        return self.population[pos]
