"""
In this file it is defined the main class to execute the GA with chromosomes with variable length.

Classes:
    :Gavl: Main class.

Example call:
    ga = Gavl()
    ga.set_hyperparameter('size_population', 100)
    ga.set_hyperparameter('min_length_chromosome', 3)
    ga.set_hyperparameter('max_length_chromosome', 10)
    def fitness(chromosome):
        return sum(chromosome)
    ga.set_hyperparameter('fitness', fitness)
    ga.set_hyperparameter('possible_genes', list(range(20)))
    ga.set_hyperparameter('termination_criteria', {'max_num_generation_reached': 150})
    ga.set_hyperparameter('repeated_genes_allowed', 0)
    ga.set_hyperparameter('minimize', 1)
    ga.set_hyperparameter('elitism_rate', 0.3)
    ga.set_hyperparameter('mutation_rate', 0.1)
    ga.set_hyperparameter('mutation_type', 'both')
    ga.set_hyperparameter('keep_diversity',5)
    ga.set_hyperparameter('show_progress',1)    
    # Launch optimization:
    best_individual = ga.optimize()
    # Get the results:
	
"""
import random
from inspect import signature
from .tools.population import Population
from .tools.individual import Individual
from .tools.generate_chromosome import generate_chromosome
from .tools.termination_criteria import check_termination_criteria
from .tools.keep_diversity import keep_diversity
from .tools.selection import roulette_selection
from .tools.pairing import pairing
from .tools.crossover import mating
from .tools.mutation import mutation


class Gavl(Population):
    """ This class is the main class for executing the GA and the only one that should be called.
    """

    def __init__(self):
        """ Constructor of the class Gavl.
        """
        # Population's attributes:
        super().__init__()
        # Hyperparameters that must be defined with the method set_hyperparameter
        self.size_population = None
        self.min_length_chromosome = None
        self.max_length_chromosome = None
        # Methods ---> Functions of the GA
        self.fitness = None  # Function that calculates the fitness of an individual.
        self.generate_new_chromosome = generate_chromosome  # Function to generate new individuals.
        self.selection = roulette_selection  # Function to perform the selection of individuals for the crossover.
        self.pairing = pairing  # Pairing method given a selection of individuals.
        self.crossover = mating  # Function that performs the crossover
        self.mutation = mutation  # Function that performs the mutation
        # Hyperparameters of the individual
        self.possible_genes = None
        self.repeated_genes_allowed = 0
        self.check_valid_individual = lambda \
                chromosome: True  # By default all the individuals are valid ---> Recommended: The invalid individuals are allowed, but they are given a very bad fitness.
        # Default values of the other hyperparameters.
        self.minimize = 1
        self.elitism_rate = 0.05
        self.mutation_rate = 0.3
        self.mutation_type = 'both'
        self.max_num_gen_changed_mutation = None  # Once max_length_chromosome is set, it will be equal to max_length_chromosome unless it is specified other num_gen_changed_mutation.
        self.termination_criteria = {'max_num_generation_reached': 100}
        self._termination_criteria_args = {'termination_criteria': 'max_num_generation_reached', 'generation_goal': 100,
                                           'generation_count': 0}  # This is a dictionary with the arguments of the termination criteria (depending on the criteria, the function to check if the termination criteria is met will receive different arguments).
        self._check_termination_criteria_function = check_termination_criteria  # Termination criteria function ---> This cannot be chaged
        self.keep_diversity = -1
        self._keep_diversity_function = keep_diversity  # Keep diversity function ---> This cannot be chaged
        # GA info:
        self.best_fitness_per_generation = []  # Attribute that holds a list of the best fitness in each generation.
        # Execution:
        self.show_progress = 1  # If show_progress = 1 it will be printed the the progress of the genetic algorithm (the generation).
        self._generation_count = 0  # Number of generations of the current population

    def set_hyperparameter(self, id_hyperparameter, value):
        """ Method to set the hyperparameters (attributes of the class Gavl).

        :param id_hyperparameter: (str) Id (name) of the hyperparameter that is going to be set.
        :param value: (...) Value of the hyperparameter that is going to be set.
        :return:
        """
        # Definition of a check of how the hyperparameters should be and a ValueError message in case these conditions are not met.
        hyperparameter_conditions = {
            'size_population': ([lambda x: type(x) == int, lambda x: x > 0,
                                 lambda x: getattr(self, 'elitism_rate', None) == 0 or
                                           getattr(self, 'elitism_rate', None) * x >= 1
                                 ],
                                "The attribute 'size_population' must be an integer greater than 0. As well, it is needed an elitism rate and a population size such that Gavl.elitism_rate = 0 or Gavl.elitism_rate*Gavl.size_population >= 1."),
            'min_length_chromosome': ([lambda x: type(x) == int, lambda x: x >= 0,
                                       lambda x: True if getattr(self, 'max_length_chromosome', None) is None
                                       else x <= getattr(self, 'max_length_chromosome', None)
                                       ],
                                      "The attribute 'min_length_chromosome' must be an integer greater or equal to 0. As well it should be smaller or equal to the attribute 'max_length_chromosome'."),
            'max_length_chromosome': ([lambda x: type(x) == int, lambda x: x >= 1,
                                       lambda x: True if getattr(self, 'min_length_chromosome', None) is None
                                       else x >= getattr(self, 'min_length_chromosome', None),
                                       lambda x: True if getattr(self, 'max_num_gen_changed_mutation', None) is None
                                       else x > getattr(self, 'max_num_gen_changed_mutation', None),
                                       lambda x: True if getattr(self, 'possible_genes', None) is None or
                                                         getattr(self, 'repeated_genes_allowed', None) == 1
                                       else x < len(getattr(self, 'possible_genes', None))
                                       ],
                                      "The attribute 'max_length_chromosome' must be an integer greater or equal to 1. As well it should be greater or equal to the attribute 'min_length_chromosome'. Note that if the the hyperparameter 'max_num_gen_changed_mutation' has been defined before 'max_length_chromosome', then 'max_length_chromosome' must be greater or equal to 'max_num_gen_changed_mutation'. NOTE that if the attribute Gavl.repeated_genes_allowed = 0, then possible_genes must be greater than the attribute Gavl.max_num_gen_changed_mutation."),
            'fitness': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 1],
                        "The method 'fitness' must be a function whose ONLY argument is the individual's chromosome (ie. fitness(chromosome)) and returns the value of the fitness. NOTE that the chromosome is a list of genes."),
            'generate_new_chromosome': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 4],
                                        "The method 'generate_new_chromosome' must be a function with the arguments min_length_chromosome (minimum length of the chromosome), max_length_chromosome (maximum length of the chromosome), possible_genes (list with the possible values that the genes can take), repeated_genes_allowed (boolean that indicates if the genes can be repeated in the chromosome or not). Note that the arguments of the function must be specified in this order."),
            'selection': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 3],
                          "The method 'selection' must be a function that receives three arguments and returns a list of the selected individuals. It receives (in this order) a list with the population (list of objects of the class Individual), the attribute self.minimize (1 -> minimize; 0 -> maximize) and the number of individuals to be selected. It must return a list with the IDs of the selected individuals (individual._id). The default selection method is roulette wheel selection."),
            'pairing': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 1],
                        "The method 'pairing' must be a function that receives one argument and returns the list of the paired individuals. It receives a list with the IDs of the selected individuals (see selection method), and returns a list of tuples with the IDs of the paired individuals. The default pairing method is random pairing."),
            'crossover': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 5],
                          "The method 'crossover' must be a function that receives five arguments and returns a list with the chromosomes of the newly crossed individuals. It must receive (in this order) a list ([(Individual_a, Individual_b) , ...]) with tuples of the paired individuals (of the class Individual), the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed) and the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False). It must return a list with the chromosomes of the newly created individuals."),
            'mutation': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 8],
                         "The method 'mutation' must be a function that receives eight arguments and returns a list with the chromosomes of the newly mutated individuals. It must receive (in this order) a list with the chromosomes of the individuals that are going to be crossed (note, it is chromosomes what this function receives, ie, list of genes, NOT objects of the class individual), an string that represents the mutation type (if the mutation method is changed, this is useless), an int that represents the maximum number of genes that it is allowed to change in a mutation (it is the attribute .max_num_gen_changed_mutation), the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed), the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False) and a list of all the allowed values that the genes can take (it is the attribute .possible_genes). It must return a list with the chromosomes of the newly mutated individuals."),
            'possible_genes': ([lambda x: type(x) == list,
                                lambda x: True if getattr(self, 'max_length_chromosome', None) is None or
                                                  getattr(self, 'repeated_genes_allowed', None) == 1
                                else len(x) >= getattr(self, 'max_length_chromosome', None)
                                ],
                               "The attribute 'possible_genes' is a list with all the possible values that the genes can take. Eg. If the chromosome is a list of indices between 0 and 999, then possible_genes = list(range(1000)). NOTE that if the attribute Gavl.repeated_genes_allowed = 0, then possible_genes must be greater than the attribute Gavl.max_length_chromosome (please, if you allow repeated genes, call the method .set_hyperparameter('repeated_genes_allowed', 1) before setting the possible genes)."),
            'repeated_genes_allowed': ([lambda x: type(x) == int or type(x) == bool,
                                        lambda x: x == 0 or x == 1 or type(x) == bool],
                                       "The attribute 'repeated_genes_allowed' must be either 0 (repeated genes not allowed in an individual) or 1 (repeated genes allowed in an individual)."),
            'check_valid_individual': ([lambda x: callable(x), lambda x: len(signature(x).parameters) == 1],
                                       "The method 'check_valid_individual' must be a function whose ONLY argument is the individual's chromosome (ie. check_valid_individual(chromosome)) and returns a boolean (True if it is a valid solution and False otherwise). Note that it is recommended NOT changing this method, and giving the invalid individuals a penalization in the fitness function. NOTE that the chromosome is a list of genes."),
            'minimize': ([lambda x: type(x) == int or type(x) == bool, lambda x: x == 0 or x == 1 or type(x) == bool],
                         "The attribute 'minimize' must be an int that can take the values 0 (maximize) or 1 (minimize). Its default value is 1."),
            'elitism_rate': ([lambda x: type(x) == float or type(x) == int, lambda x: 0 <= x <= 1,
                              lambda x: True if getattr(self, 'size_population', None) is None
                              else x == 0 or getattr(self, 'size_population', None) * x >= 1
                              ],
                             "Elitism rate must be a number that only takes values between 0 and 1. As well, it is needed an elitism rate and a population size such that Gavl.elitism_rate = 0 or Gavl.elitism_rate*Gavl.size_population >= 1."),
            'mutation_rate': ([lambda x: type(x) == float or type(x) == int, lambda x: 0 <= x <= 1],
                              "Mutation rate must be a number between 0 and 1."),
            'mutation_type': ([lambda x: type(x) == str, lambda x: x in ['mut_gene', 'addsub_gene', 'both']],
                              "The attribute 'mutation_type' must be a string that takes the values 'mut_gene', 'addsub_gene' or 'both'."),
            'max_num_gen_changed_mutation': ([lambda x: type(x) == int,
                                              lambda x: True if getattr(self, 'max_length_chromosome', None) is None
                                              else x < getattr(self, 'max_length_chromosome', None)
                                              ],
                                             "The MAXIMUM number of genes changed in each mutation (attribute 'max_num_gen_changed_mutation') must be an integer between 1 and the maximum length of the chromosome. Note that by default 'max_num_gen_changed_mutation' will be int(max_length_chromosome/3 + 1)."),
            'termination_criteria': ([lambda x: type(x) == dict, lambda x: len(x) == 1,
                                      lambda x: list(x.keys())[0] in ['goal_fitness_reached',
                                                                      'max_num_generation_reached'],
                                      lambda x: type(list(x.values())[0]) == int or type(list(x.values())[0]) == float
                                      ],
                                     "The attribute 'termination_criteria' must be a dictionary with either the value '{'max_num_generation_reached': number of generations}' or '{'goal_fitness_reached': goal fitness}'. As well the number of generations or the goal fitness must be of type int or float."),
            'keep_diversity': ([lambda x: type(x) == int, lambda x: x != 0, lambda x: x >= -1],
                               "The attribute 'keep_diversity' must be an int that can either take the value -1 (default value) if it is not wanted use a technique to keep the diversity of the population or it can take an int value that represents the number of generations in which this method is applied (eg. if it takes the value 2, keep diversity will be applied every two generations). Note that 'keep_diversity' cannot be equal to 0."),
            'show_progress': (
                [lambda x: type(x) == int or type(x) == bool, lambda x: x == 0 or x == 1 or type(x) == bool],
                "The attribute 'show_progress' must be an int that can take the values 0 (do not show progress) or 1 (show progress). Its default value is 0.")
        }
        if id_hyperparameter not in list(hyperparameter_conditions.keys()):
            raise ValueError("""The parameter id_hyperparameter of the method set_hyperparameter() must be one of the parameters of the next list:
                * 'size_population': Int that represents the size of the population.
                * 'min_length_chromosome': Int that represents the minimum length of a chromosome.
                * 'max_length_chromosome': Int that represents the maximum length of the chromosome.
                * 'fitness': Function to evaluate the fitness. Its ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness. NOTE that the chromosome is a list of genes.
                * 'generate_new_chromosome': Function to create a new chromosome. It receives the four arguments (in this order) min_length_chromosome (minimum length of the chromosome), max_length_chromosome (maximum length of the chromosome), possible_genes (list with the possible values that the genes can take), repeated_genes_allowed (boolean that indicates if the genes can be repeated in the chromosome or not). This function must return a list of genes. Note that the arguments of the function must be specified in this order.
                * 'selection': Function to perform the selection method. It must be a function that receives three arguments and returns a list of the selected individuals. It receives (in this order) a list with the population (list of objects of the class Individual), the attribute self.minimize (1 -> minimize; 0 -> maximize) and the number of individuals to be selected. It must return a list with the IDs of the selected individuals (individual._id). The default selection method is roulette wheel selection.
                * 'pairing': Function to perform the pairing method. It must be a function that receives one argument and returns the list of the paired individuals. It receives a list with the IDs of the selected individuals (see selection method), and returns a list of tuples with the IDs of the paired individuals. The default pairing method is random pairing.
                * 'crossover': Function to perform the crossover method. It must be a function that receives five arguments and returns a list with the chromosomes of the newly crossed individuals. It must receive (in this order) a list ([(Individual_a, Individual_b) , ...]) with tuples of the paired individuals (of the class Individual), the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed) and the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False). It must return a list with the chromosomes of the newly created individuals.
                * 'mutation': Function to perform the mutation method. It must be a function that receives eight arguments and returns a list with the chromosomes of the newly mutated individuals. It must receive (in this order) a list with the chromosomes of the individuals that are going to be crossed (note, it is chromosomes what this function receives, ie, list of genes, NOT objects of the class individual), an string that represents the mutation type (if the mutation method is changed, this is useless), an int that represents the maximum number of genes that it is allowed to change in a mutation (it is the attribute .max_num_gen_changed_mutation), the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed), the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False) and a list of all the allowed values that the genes can take (it is the attribute .possible_genes). It must return a list with the chromosomes of the newly mutated individuals.
                * 'possible_genes': list with all the possible values that the genes can take.
                * 'repeated_genes_allowed': Int that represents if an individual can have repeated genes (repeated_genes_allowed = 1) or not (repeated_genes_allowed = 0). By default it is 0.
                * 'check_valid_individual': function whose ONLY argument is the individual's chromosome (ie. check_valid_individual(chromosome)) and returns a boolean (True if it is a valid solution and False otherwise). Note that it is recommended NOT changing this method, and giving the invalid individuals a penalization in the fitness function. NOTE that the chromosome is a list of genes.
                * 'minimize': Int that represents if the fitness will be minimized (minimize = 1) or maximized (minimize = 0). By default minimize = 1.
                * 'elitism_rate': Number between 0 and 1 that represents the elitism rate. By default elitism_rate = 0.05.
                * 'mutation_rate': Number between 0 and 1 that represents the mutation rate. By default mutation_rate = 0.3.
                * 'mutation_type': String that represents the mutation type. It can only take the values 'mut_gene', 'addsub_gene' or 'both'. By default mutation_type = 'both'.
                * 'max_num_gen_changed_mutation': Int that represents the MAXIMUM number of genes changed in each mutation. By default it is int(max_length_chromosome/3 + 1).
                * 'termination_criteria': The attribute 'termination_criteria' must be a dictionary that represents the termination criteria, with either the value '{'max_num_generation_reached': number of generations}' or '{'goal_fitness_reached': goal fitness}.
                * 'keep_diversity': Int that represents every how many generations the diversity techniques must be applied. Its default value is -1, which means that no diversity techniques will be applied.
                * 'show_progress': Int that indicates if it is wanted to show the progress. It can take the values 0 (do not show progress) or 1 (show progress). Its default value is 1.
                """)
        else:
            try:
                conditions = hyperparameter_conditions[id_hyperparameter]
                if False in list(map(lambda x: x(value), conditions[0])):  # Check all the conditions
                    raise ValueError(conditions[1])  # There has been an error
                else:
                    setattr(self, id_hyperparameter, value)  # Set the attribute
                    if id_hyperparameter == 'max_length_chromosome' and getattr(self, 'max_num_gen_changed_mutation',
                                                                                None) is None:  # If the hyperparameter that is being set is 'max_length_chromosome' and it has not been defined the hyperparameter 'max_num_gen_changed_mutation', then it is set to 'max_length_chromosome'.
                        setattr(self, 'max_num_gen_changed_mutation', int(value / 3 + 1))  # Set the attribute
                    elif id_hyperparameter == 'termination_criteria':  # Actualize the arguments of the function to check the termination criteria
                        if list(value.keys())[0] == 'goal_fitness_reached':
                            self._termination_criteria_args = {'termination_criteria': 'goal_fitness_reached',
                                                               'goal_fitness': list(value.values())[0],
                                                               'generation_fitness': 0,
                                                               'minimize': self.minimize}
                        elif list(value.keys())[0] == 'max_num_generation_reached':
                            self._termination_criteria_args = {'termination_criteria': 'max_num_generation_reached',
                                                               'generation_goal': list(value.values())[0],
                                                               'generation_count': 0}
                        else:  # Double check
                            ValueError(
                                "The termination criteria argument can have either the parameters 'goal_fitness_reached' or 'max_num_generation_reached'.")
                    elif id_hyperparameter == 'minimize' and self._termination_criteria_args[
                        'termination_criteria'] == 'goal_fitness_reached':
                        self._termination_criteria_args['minimize'] = value
            except ValueError:
                raise ValueError(conditions[1])  # There has been an error
            except Exception as e:
                print(str)
                raise (type(e))(str(e) + '\n' + conditions[1])  # There has been an error

    def optimize(self):
        """ This method is the one that is called to start the optimization (start the genetic algorithm). Before calling this method it must be defined Gavl.size_population, Gavl.min_length_chromosome, Gavl.max_length_chromosome, Gavl.fitness and Gavl.possible_genes.

        :return:
            (Individual) Best individual.
        """
        # First, check if the needed attributes are defined.
        if self.fitness is None:
            raise AttributeError(
                "The method fitness has to be defined before you can call this method. For doing so call the method Gavl.set_hyperparameter('fitness', value), where value is a function whose ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness.")
        elif self.size_population is None:
            raise AttributeError(
                "The attribute 'size_population' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('size_population', value).")
        elif self.min_length_chromosome is None:
            raise AttributeError(
                "The attribute 'min_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('min_length_chromosome', value).")
        elif self.max_length_chromosome is None:
            raise AttributeError(
                "The attribute 'max_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('max_length_chromosome', value).")
        elif self.possible_genes is None:
            raise AttributeError(
                "The attribute 'possible_genes' has to be defined before you can call this method. It must be a list with all the possible values that the genes can take.")
        else:
            # Start the algorithm
            self._generation_count = 0  # Start the generation counter
            self.best_fitness_per_generation = []  # Empty list
            # Create the population
            self._Population__generate_population()
            self._Population__calculate_fitness_and_sort()
            if self._check_termination_criteria_function(self._termination_criteria_args):
                return self.best_individual()
            while not self._check_termination_criteria_function(self._termination_criteria_args):
                self._generation_count += 1  # Start the generation counter
                if self.show_progress:
                    print('Generation: {}'.format(self._generation_count))
                new_population = self._Population__get_next_generation()  # Calculate the next generation.
                self._Population__kill_and_reset_whole_population(
                    new_population)  # Set the next generation.
                if self._generation_count % self.keep_diversity == 0 and self.keep_diversity > 0:
                    self._Population__calculate_fitness_and_sort()  # Calculate the fitness and sort the population
                    # Keep diversity protocol:
                    new_diverse_population = self._keep_diversity_function(self.population,
                                                                           self.generate_new_chromosome,
                                                                           self.min_length_chromosome,
                                                                           self.max_length_chromosome,
                                                                           self.possible_genes,
                                                                           self.repeated_genes_allowed,
                                                                           self.check_valid_individual)
                    self._Population__kill_and_reset_whole_population(
                        new_diverse_population)  # Set the next generation.
                self.__update_termination_criteria_args()  # Update the termination criteria arguments
                self.best_fitness_per_generation.append(
                    self.best_individual().fitness_value)  # Get the best fitness value per generation
            self._Population__calculate_fitness_and_sort()  # Calculate the fitness and sort the population
            return self.best_individual()

    def _Population__get_next_generation(self):
        """ This method is used to calculate the next generation.

        :return:
            * :new_generation: (list of chromosomes) List with the chromosomes of the next generation.
        """
        # First, sort the population by its fitness:
        self._Population__calculate_fitness_and_sort()  # Calculate the fitness and sort the population
        # List with the individuals of the next generation:
        new_generation = []
        # Get the size of the groups for the new population:
        size_elitism = int(len(self.population) * self.elitism_rate)  # Number of individuals chosen as elite
        size_crossover = int(len(self.population)) - size_elitism  # Number of individuals chosen for crossover
        if (size_crossover % 2) == 1:  # Odd number of Individuals for crossover
            size_elitism += 1  # Add one individual to elitism
            size_crossover -= 1  # Remove one from crossover
        # Elite:
        elite = [individual.chromosome for individual in self.population[:size_elitism]]
        new_generation.extend(
            elite)  # Add elite ---> Note that the population is sorted by fitness before calling this function.
        # Crosover:
        selected_individuals = self.selection(self.population, self.minimize, size_crossover)  # 1. Roulette selection
        paired_ids = self.pairing(selected_individuals)  # 2. Perform the pairing
        list_of_paired_ind = [(self.get_individual_by_id(id_a).chromosome, self.get_individual_by_id(id_b).chromosome)
                              for id_a, id_b in paired_ids]  # List with the chromosomes of the paired individuals
        new_crossed_ind = self.crossover(list_of_paired_ind,
                                         self.min_length_chromosome,
                                         self.max_length_chromosome,
                                         self.repeated_genes_allowed,
                                         self.check_valid_individual)  # 3. Get the new chromosomes already 'crossovered'
        for new_individual in new_crossed_ind:  # 4. Add crossovered individuals
            if type(new_individual) == Individual:
                new_generation.append(new_individual.chromosome)
            elif type(new_individual) == list:
                new_generation.append(new_individual)
            else:  # In case it is incorrectly redefined the crossover method
                raise ValueError(
                    'The method crossover must return a list with the chromosomes of the newly crossed individuals.')
        # Mutation:
        size_mutation = int(len(self.population) * self.mutation_rate)  # Size of the mutation
        if size_mutation >= len(new_generation) - size_elitism:
            size_mutation = int(len(new_generation) - size_elitism) - 1
        indices_mutation = random.sample(range(size_elitism, len(new_generation)),
                                         size_mutation)  # Get the indices of the individuals that are going to be mutated
        chromosomes_to_mutate = [new_generation[i] for i in
                                 indices_mutation]  # Get the chromosomes that are going to be mutated
        mutated_individuals = self.mutation(chromosomes_to_mutate, self.mutation_type,
                                            self.max_num_gen_changed_mutation, self.min_length_chromosome,
                                            self.max_length_chromosome, self.repeated_genes_allowed,
                                            self.check_valid_individual, self.possible_genes)
        for i in indices_mutation:  # Add the newly mutated individuals to the population
            m_ind = mutated_individuals.pop()
            if type(m_ind) == Individual:
                new_generation[i] = m_ind.chromosome
            elif type(m_ind) == list:
                new_generation[i] = m_ind
            else:  # In case it is incorrectly redefined the mutation method
                raise ValueError(
                    'The method mutation must return a list with the chromosomes of the newly mutated individuals.')
        return new_generation

    def _Population__calculate_fitness_population(self):
        """ This method calculates the fitness of all the individuals of the population and then it sets this attribute to each individual.
        """
        # First, check if the needed attributes are defined.
        if self.fitness is None:
            raise AttributeError(
                "The method fitness has to be defined before you can call this method. For doing so call the method Gavl.set_hyperparameter('fitness', value), where value is a function whose ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness.")
        elif not len(self.population):
            raise AttributeError('The population has not been generated yet.')
        else:
            for ind in self.population:
                ind.calculate_fitness(self.fitness)

    def _Population__generate_population(self):
        """ This method should generate a new population and add it to the attribute population.
        """
        if self.fitness is None:
            raise AttributeError(
                "The method fitness has to be defined before you can call this method. For doing so call the method Gavl.set_hyperparameter('fitness', value), where value is a function whose ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness.")
        elif self.size_population is None:
            raise AttributeError(
                "The attribute 'size_population' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('size_population', value).")
        elif self.min_length_chromosome is None:
            raise AttributeError(
                "The attribute 'min_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('min_length_chromosome', value).")
        elif self.max_length_chromosome is None:
            raise AttributeError(
                "The attribute 'max_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('max_length_chromosome', value).")
        elif self.possible_genes is None:
            raise AttributeError(
                "The attribute 'possible_genes' has to be defined before you can call this method. It must be a list with all the possible values that the genes can take.")
        else:
            while len(self.population) < self.size_population:
                # Create new individual:
                new_ind = self.generate_new_chromosome(self.min_length_chromosome, self.max_length_chromosome,
                                                       self.possible_genes, self.repeated_genes_allowed)
                if self.check_valid_individual(
                        new_ind):  # If check_valid_individual is very restrictive it may cause this step to take a VEEERY LOOONG time
                    self.add_individual(new_ind)

    def _Population__sort_population(self):
        """ This method sorts the population by its fitness.
        """
        if not len(self.population):
            raise AttributeError('The population has not been generated yet.')
        else:
            if any([ind.fitness_value is None for ind in self.population]):
                self._Population__calculate_fitness_population()  # Calculate the fitness of the individuals if it hasn't been done yet
            if self.minimize:  # If minimize == 1
                self.population.sort(key=lambda x: x.fitness_value,
                                     reverse=False)  # Sort from best fitness to worst fitness
            else:
                self.population.sort(key=lambda x: x.fitness_value,
                                     reverse=True)  # Sort from best fitness to worst fitness

    def __update_termination_criteria_args(self):
        """ This method is called to update the termination criteria arguments. It MUST be called in each new generation of the genetic algorithm.

        :return:
        """
        if list(self.termination_criteria.keys())[0] == 'max_num_generation_reached':
            self._termination_criteria_args['generation_count'] = self._generation_count
        elif list(self.termination_criteria.keys())[0] == 'goal_fitness_reached':
            self._termination_criteria_args['generation_fitness'] = self.best_individual().fitness_value

    def add_individual(self, individual):
        """ Method to add a new individual to the population. This method is overriding the method with the same name of the class Population ---> This is done to check that it is not added more individuals than self.size_population.

        :param individual: (Individual or list) individual of the class Individual or list representing the chromosome.
        :return:
        """
        if self.size_population is None:
            raise AttributeError(
                'The size of the population has not been generated yet. Please, define it by calling Gavl.set_hyperparameter("size_population", size).')
        else:
            if type(individual) == list:
                ind = individual.copy()
            elif type(individual) == Individual:
                ind = individual.chromosome.copy()
            else:
                raise ValueError(
                    'The given individual is not valid. It must be a list of genes or an object of the class Individual.')
            if self.check_valid_individual(ind):
                # Check that the number of individuals is lower than the maximum allowed size of the population
                if len(self.population) < self.size_population:
                    super().add_individual(ind)
                else:
                    raise AttributeError(
                        'The population is already filled (there are {} individuals, which is the specified maximum size of the population). If wanted more individuals change the attribute size_population by calling Gavl.set_hyperparameter("size_population", size).'.format(
                            self.size_population))
            else:
                raise ValueError(
                    'The given individual is not valid. Please, check the method check_valid_individual or change the individual.')

    def best_individual(self):
        """ This method returns the best individual.

        :return:
            * :individual: (Individual) Best individual.
        """
        if not len(self.population):
            raise AttributeError('The population has not been generated yet.')
        self._Population__sort_population()  # Order the population from best individual to worst.
        return self.population[0]

    def historic_fitness(self):
        """ This method returns the best fitness value in each generation.

        :return:
            * :bfv: (list of floats) This method returns a list with the best fitness value in each generation.
        """
        if not self.best_fitness_per_generation:
            raise ValueError(
                'You must make the optimization (call the method .optimize()) before you can call this method.')
        return self.best_fitness_per_generation

    def get_results(self):
        """ This method is used to get the results of the optimization process, once it has finished.

        :return:
            * :best_individual: (Individual) Best individual.
            * :population: (list of Individuals) List with all the individuals of the last generation.
            * historic_fitness: (list of floats) List with the best fitness value in each generation.
        """
        # First, check if the needed attributes are defined.
        if self.fitness is None:
            raise AttributeError(
                "The method fitness has to be defined before you can call this method. For doing so call the method Gavl.set_hyperparameter('fitness', value), where value is a function whose ONLY argument is the individual's chromosome (fitness(chromosome)) and returns the value of the fitness.")
        elif self.size_population is None:
            raise AttributeError(
                "The attribute 'size_population' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('size_population', value).")
        elif self.min_length_chromosome is None:
            raise AttributeError(
                "The attribute 'min_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('min_length_chromosome', value).")
        elif self.max_length_chromosome is None:
            raise AttributeError(
                "The attribute 'max_length_chromosome' has to be defined before you can call this method. It must be an integer greater than 0 and it can be set by calling the method Gavl.set_hyperparameter('max_length_chromosome', value).")
        elif self.possible_genes is None:
            raise AttributeError(
                "The attribute 'possible_genes' has to be defined before you can call this method. It must be a list with all the possible values that the genes can take.")
        best_individual = self.best_individual()
        population = self.population
        historic_fitness = self.historic_fitness()
        return best_individual, population, historic_fitness
