# Genetic algorithm with chromosomes with variable length

This repository contains a Python framework to launch a genetic algorithm based in individuals with chromosomes with variable length. The algorithm receives a fitness function and a set of configurable hyperparameters, and then it performs the optimization. The individuals are based in lists of variable length (the minimum and maximum length are configurable parameters) of non-ordered genes with or without repetition (this parameter is configurable as well). The crossover consists in taking a random number of genes of each of the parents to create the child (explained in the section of the crossover). As well, two different mutation algorithms are proposed.



## How to use:

  1. Download the files:
```shell
	git clone https://github.com/IgnacioGarrido/GA_chromosomeWithVariableLength.git
```

  2. Import the Gavl class from the file Gavl.py. Note that this is the main and only class needed to execute the genetic algorithm. All the hyperparameters should be passed to an object of this class.
```python
	from path.Gavl import Gavl
```

  3. Create an instance of Gavl():
```python
	ga = Gavl()
```

  4. Configure the hyperparameters that you want:
```python
    ga.set_hyperparameter('size_population', 100)  # Set the size of the population to 100 individuals
    ga.set_hyperparameter('min_length_chromosome', 3)  # Set the minimum length of the individual to 3 genes
    ga.set_hyperparameter('max_length_chromosome', 10)  # Set the maximum length of the individual to 10 genes

    # The fitness function receives the chromosome, which is a list with the genes. It must return the fitness value (number).
    def fitness(chromosome):
        return sum(chromosome)

    ga.set_hyperparameter('fitness', fitness)  # Set the fitness function
    ga.set_hyperparameter('possible_genes', [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])  # The possible values that the genes can take are the integer numbers between 0 and 15. Notice that this can be a list of integers, strings, ... 
    ga.set_hyperparameter('termination_criteria', {'max_num_generation_reached': 150})  # Set the termination criteria to maximum number of generations (in this case 150 generations).
    ga.set_hyperparameter('repeated_genes_allowed', 0)  # In this case the repeated genes are not allowed
    ga.set_hyperparameter('minimize', 1)  # In this case the goal is to minimize the fitness
    ga.set_hyperparameter('elitism_rate', 0.1)  # Set an elitism rate of 0.1
    ga.set_hyperparameter('mutation_rate', 0.2)  # Set a mutation rate of 0.2
    ga.set_hyperparameter('mutation_type', 'both')  # Set the mutation type (see documentation bellow)
    ga.set_hyperparameter('keep_diversity',5)  # Set a mechanism to keep the diversity (see documentation bellow)
    ga.set_hyperparameter('show_progress',1)  # Show progress of the genetic algorithm (this will print the generation number)
```

  5. Launch the optimization (notice that the method ```optimize()``` returns the best individual (best solution found for the current problem):
```python
    best_individual = ga.optimize()
```

  6. If needed, get more results:
```python
    best_individual, population, historic_fitness = ga.get_results()
```

  7. As well, if needed any changes, fork the repository and make the all the modifications you want for your project. This is open source software and any additional adaptations are welcome :).



## Set hyperparameters:

As shown in the example, all the hyperparameters can be set by calling the method ```set_hyperparameter()```. This method receives two arguments, the id of the hyperparameter (id_hyperparameter) and its value. Next it is shown all the possible configurable hyperparameters and an example of how they can be configured.



### Parameters that MUST be tunned:

Next it is shown a list with the hyperparameters that must be tunned before launching the genetic algorithm. These hyperparameters have not default value, and if it is tried to launch the program before tunning them an AttributeError will occur.

  * __'size_population'__: Int that represents the size of the population. It can be set by calling the method ```.set_hyperparameter('size_population', 100)``` ---> _Note that in this example the size of the population is being set to 100._
 
  * __'min_length_chromosome'__: Int that represents the minimum length of a chromosome. It can be set by calling the method ```.set_hyperparameter('min_length_chromosome', 3)``` ---> _Note that in this example the minimum length of the chromosome is being set to 3._
 
  * __'max_length_chromosome'__: Int that represents the maximum length of the chromosome. It can be set by calling the method ```.set_hyperparameter('max_length_chromosome', 10)``` ---> _Note that in this example the maximum length of the chromosome is being set to 10._
  
  * __'fitness'__: Function to evaluate the fitness. Its ONLY argument is the individual's chromosome (fitness(chromosome)), which is a list of genes, and returns the value of the fitness. It can be set by calling the method ```.set_hyperparameter('fitness', fun_fitness)``` ---> _Where fun_fitness is a function of the type ```fun_fitness(chromosome)```. For more information see the fitness chapter bellow._

  * __'possible_genes'__: list with all the possible values that the genes can take. Following the above example the possible genes could take the values ```pg = ['pen', 'pencil', 'food', 'rubber', 'book']```. It can be set by calling the method ```.set_hyperparameter('possible_genes', pg)``` ---> _Where pg is the previously defined list of items._



### Parameters that CAN be tunned:

  * __'generate_new_chromosome'__: Function to create a new chromosome. It receives the four arguments (in this order) _min_length_chromosome_ (minimum length of the chromosome), _max_length_chromosome_ (maximum length of the chromosome), _possible_genes_ (list with the possible values that the genes can take), _repeated_genes_allowed_ (boolean that indicates if the genes can be repeated in the chromosome or not). This function must return a list of genes. Note that the arguments of the function must be specified in this order. ---> _It can be set by calling the method ```.set_hyperparameter('generate_new_chromosome', generate_chromosome)```, where generate_chromosome is a function of the type ```generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes, repeated_genes_allowed)```. The default value is a function that takes a random length of the chromosome within the limits, and then it randomly selects the genes among the possibilities. The current default implementation of this function is shown bellow:_

```python
import random 


def generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes, repeated_genes_allowed):
	""" Function called to create a new individual (its chromosome). It randomly chooses its length (between min_length_chromosome and min_length_chromosome), and it randomly chooses genes among the list of possible_genes.

	:param min_length_chromosome: (int) Minimum allowed length of the chromosome.
	:param max_length_chromosome: (int) Maximum allowed length of the chromosome.
	:param possible_genes: (list of ...) List with the all the possible values that the genes can take.
	:param repeated_genes_allowed: (bool) It is a boolean that indicates whether the genes can be repeated in the chromosome (repeated_genes_allowed = 1) or they cannot be repeated (repeated_genes_allowed = 0).
	:return:
	    * (list of genes) List that represents the chromosome.
	"""
	# Choose a random number of genes
	number_of_genes = random.randrange(min_length_chromosome, max_length_chromosome + 1)
	# Create new chromosome:
	if repeated_genes_allowed:
	    chromosome = random.choices(possible_genes, weights=None, k=number_of_genes)
	    return chromosome
	else:
	    possible_genes_aux = possible_genes.copy()
	    random.shuffle(possible_genes_aux)
	    return possible_genes_aux[:number_of_genes]
```
  
  * __'selection'__: Function to perform the selection method. It must be a function that receives three arguments and returns a list of the selected individuals. It receives (in this order) a list with the population (list of objects of the class Individual), the attribute self.minimize (1 -> minimize; 0 -> maximize) and the number of individuals to be selected. It must return a list with the IDs of the selected individuals (individual._id). The default selection method is roulette wheel selection. ---> _It can be set by calling the method ```.set_hyperparameter('selection', new_selection_function)```. Its default value is roulette wheel selection (see selection chapter below)._
  
  * __'pairing'__: Function to perform the pairing method. It  must be a function that receives one argument and returns the list of the paired individuals. It receives a list with the IDs of the selected individuals (see selection method), and returns a list of tuples with the IDs of the paired individuals. The default pairing method is random pairing. ---> _It can be set by calling the method ```.set_hyperparameter('selection', new_pairing_function)```. Its default value is random pairing (see pairing chapter below)._

  * __'crossover'__: Function to perform the crossover method. It must be a function that receives five arguments and returns a list with the chromosomes of the newly crossed individuals. It must receive (in this order) a list ([(Individual_a, Individual_b) , ...]) with tuples of the paired individuals (of the class Individual) , the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed) and the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False). It must return a list with the chromosomes of the newly created individuals. ---> _It can be set by calling the method ```.set_hyperparameter('crossover', new_crossover_function)```. Its default value is explained in the crossover chapter below._

  * __'mutation'__: Function to perform the mutation method. It must be a function that receives eight arguments and returns a list with the chromosomes of the newly mutated individuals. It must receive (in this order) a list with the chromosomes of the individuals that are going to be crossed (note, it is chromosomes what this function receives, ie, list of genes, NOT objects of the class individual), an string that represents the mutation type (if the mutation method is changed, this is useless), an int that represents the maximum number of genes that it is allowed to change in a mutation (it is the attribute .max_num_gen_changed_mutation), the minimum allowed length of the chromosome, the maximum length of the chromosome, a boolean that indicates if the genes can be repeated (1 = repeated genes allowed), the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False) and a list of all the allowed values that the genes can take (it is the attribute .possible_genes). It must return a list with the chromosomes of the newly mutated individuals. ---> _It can be set by calling the method ```.set_hyperparameter('mutation', new_mutation_function)```. Its default value is explained in the mutation chapter below._

  * __'repeated_genes_allowed'__: Int that represents if an individual can have repeated genes (repeated_genes_allowed = 1) or not (repeated_genes_allowed = 0). ---> _It can be set by calling the method ```.set_hyperparameter('repeated_genes_allowed', 0)```. Its default value is 0._
  
  * __'check_valid_individual'__: function whose ONLY argument is the individual's chromosome (ie. check_valid_individual(chromosome)) and returns a boolean (True if it is a valid solution and False otherwise). Note that it is recommended NOT changing this method, and giving the invalid individuals a penalization in the fitness function (as shown in the fitness function example above). Note that the individual's chromosome is a list of genes. ---> _It can be set by calling the method ```.set_hyperparameter('check_valid_individual', check_valid_individual)```, where check_valid_individual is a function of the type ```check_valid_individual(chromosome)```. Its default value is a function that always returns True (ie. all the individuals are accepted)._
 
  * __'minimize'__: Int that represents if the fitness will be minimized (minimize = 1) or maximized (minimize = 0). ---> _It can be set by calling the method ```.set_hyperparameter('minimize', 1)```. Its default value is 1._
 
  * __'elitism_rate'__: Number between 0 and 1 that represents the elitism rate. ---> _It can be set by calling the method ```.set_hyperparameter('elitism_rate', 0.05)```. Its default value is 0.05._
 
  * __'mutation_rate'__ : Number between 0 and 1 that represents the mutation rate. ---> _It can be set by calling the method ```.set_hyperparameter('mutation_rate', 0.3)```. Its default value is 0.3._

  * __'mutation_type'__: String that represents the mutation type. It can ONLY take the values 'mut_gene', 'addsub_gene' or 'both'. Bellow it is explained what these mutation parameter mean (see the mutation chapter).  ---> _It can be set by calling the method ```.set_hyperparameter('mutation_type', 'both')```. Its default value is 'both'._

  * __'max_num_gen_changed_mutation'__: Int that represents the MAXIMUM number of genes changed in each mutation. ---> _It can be set by calling the method ```.set_hyperparameter('max_num_gen_changed_mutation', 5)```. Its default value is int(max_length_chromosome/3 + 1)._

  * __'termination_criteria'__: Dictionary that represents the termination criteria. The attribute 'termination_criteria' must be a dictionary with either the value '{'max_num_generation_reached': number of generations}' (maximum number of generations reached) or '{'goal_fitness_reached': goal fitness} (maximum goal value reached). ---> _It can be set by calling the method ```.set_hyperparameter('termination_criteria', {'max_num_generation_reached': 100})```. Its default value is ```{'max_num_generation_reached': 100}```, ie, the algorithm stops when computed 100 generations._

  * __'keep_diversity'__: Int that represents if it is wanted to apply the diversity techniques (keep_diversity = 1) or not (keep_diversity = 0). The diversity techniques are explained in the chapter "Keep diversity". ---> _It can be set by calling the method ```.set_hyperparameter('keep_diversity', x)``` which means that every x generations the diversity techniques will be applied. For example ```.set_hyperparameter('keep_diversity', 5)``` means that every 5 generations the diversity techniques are applied. Its default value is -1, which means that NO diversity techniques are applied._

  * __'show_progress'__: Int that indicates if it is wanted to show the progress. It can take the values 0 (do not show progress) or 1 (show progress). ---> _It can be set by calling the method ```.set_hyperparameter('show_progress', 1)```. Its default value is 1._



## The algorithm

Genetic algorithms (GA) are a well-known optimization tool used to solve problems in which the usage of the brute force may lead to inadmissible execution time. GA use nature-based heuristics that give a great approximation to the optimal solution and, luckily, converge to the best solution. 

Concretely, this project tries to solve problems in which the chromosome of the GA has variable length.



### Problems that can be solved with this GA

This optimization tool can be used to solve problems in which the chromosome has variable length and it can be represented as a list of genes. The nature of the genes and the fitness function must be defined before starting the algorithm, and it must be set a range of the possible lengths of the chromosome (the more fine is configured this parameter, the easier it is for the algorithm to converge).



### Individual

The individuals are represented by the class Individual with the following definition:

```python
import uuid

class Individual:
  def __init__(self, chromosome):
    """ Constructor.

    :param chromosome: (list of genes) Chromosome of the individual that is being created.
    """
    if type(chromosome) != list:
    	raise AttributeError('The chromosome must be a list of genes')
    else:
	  	self.chromosome = chromosome
	    self._id = str(uuid.uuid4())  # Unique ID for each individual ---> It is unique because uuid uses the time component for creating this id
	    self.fitness_value = None  # This attribute will be filled when the individual is evaluated
	    self.normalized_fitness_value = None  # This attribute holds the normalized value of the fitness in comparison to the rest of the population
	    self.inverse_normalized_fitness_value = None  # This attribute holds the inverse normalized value of the fitness in comparison to the rest of the population
```

Each individual has five attributes:

  * __chromosome__: (list) Is a variable length list (length between [_min_number_of_genes, max_number_of_genes_]). Examples of chromosomes are ```[1,10,6,3]```, ```['apple', 'banana', 'orange']``` and ```[{'a': 1, 'b': 2}, {'a': 2, 'b': 4}, {'a': 3, 'b': 2}]```.

  * __\_id__: (str) Unique string that identifies each individual. To create the ID it is used the module ```uuid```.

  * __fitness_value__: (float) Fitness value according to the defined fitness function (see fitness section below).

  * __normalized_fitness_value__: (float) Normalized fitness value with respect to the rest of the population (the individual with the highest fitness value will have a value of 1 in this attribute and the one with the lowest fitness value will have a value of 0).

  * __inverse_normalized_fitness_value__: (float) Inverse normalized fitness value with respect to the rest of the population (the individual with the lowest fitness will have a value of 1 in this attribute and the one with the highest fitness value will have a value of 0).

As well, each individual has the following methods:

  * __set_new_chromosome__:  Method to set a chromosome.

  * __set_fitness_value__: Method to set the fitness value.

  * __set_normalized_fitness_value__: Method to set the normalized fitness value.

  * __set_inverse_normalized_fitness_value__:  Method to set the inverse normalized fitness value.

  * __calculate_fitness__: Method to calculate the fitness of the individual.
  
  * __kill_and_reset__: This is a bit of a "tricky" method. If a new individual of a new generation is going to be created, it will be computationally much faster to get the object of an individual that is going to be killed and reset its values to the ones of the new individual. ---> This will be done instead of directly creating a completely new individual with the new data and stop referencing the old individual.



### Population

The population consists in a list of individuals. The population can be accessed with the method ```.get_population()```, which will return a list of ordered (from best to worst) individuals of the class Individual (shown above). 

There are three ways to access to the best individual's attributes:

```python
ga = Gavl()

# Configure the hyperparameters
best_individual = ga.optimize()

# You can get the best fitness and chromosome with the returned best individual of the method .optimize()
best_chromosome = best_individual.chromosome
best_fitness = best_individual.fitness_value

# You can also use the method best_individual(), which returns the best individual
bi = ga.best_individual()
best_chromosome = bi.chromosome
best_fitness = bi.fitness_value

# You can also get the whole ordered population (from best to worst individual) and get the first one
ordered_population = ga.get_population()
best_ind = ordered_population[0] # Get the first individual, ie, the best one
best_chromosome = best_ind.chromosome
best_fitness = best_ind.fitness_value
```



### Selection - Roulette wheel (Fitness proportionate)

The only available selection method in this project is roulette wheel selection. This method takes the whole population, and based in their normalized fitness (if the goal is maximize) or based in their inverse normalized fitness (if the goal is minimize) it performs roulette wheel selection.

This method is based in the selection of each individual in dependence of its fitness value. The better is the individual (better fitness), the more probable is to choose that individual for the mating process. Thus, the probability of choosing an individual for the mating process is:

<img src="https://render.githubusercontent.com/render/math?math=P_{i} = \frac{f_{i}}{\sum_{j=1}^{n} f_{j}}">

<span style="color:lightgray"> _Where f <sub>i</sub> is the normalized fitness (or inverse normalized if the goal is minimizing the fitness) of the individual i._</span>

In each generation a random selection process based on the fitness value is performed for the subsequent pairing and crossover.

_\* If other selection method is wanted, it can be set by calling the method ```.set_hyperparameter('selection', new_selection_function)```, where new_selection_function is the function that performs this new selection method. It must be a function that receives three arguments and returns a list of the selected individuals. It receives (in this order) a list with the population (list of objects of the class Individual), the attribute 'minimize' (1 -> minimize; 0 -> maximize) and the number of individuals to be selected. It must return a list with the IDs of the selected individuals (attribute individual.\_id)._

_\** If more information is wanted, go to the docs and the definition of the roulette wheel selection function in ga_variable_length/tools/selection.py._



### Pairing - Random pairing

The only available pairing method in this project is random pairing. Random pairing is based in, as its name says, given a selection of individuals, randomly pair them. In other words, given a list of the selected individuals (selection function) randomly pick the in pairs.

In each generation a random pairing process based in the selection is performed for the subsequent crossover.

_\* If other pairing method is wanted, it can be set by calling the method ```.set_hyperparameter('pairing', new_pairing_function)```, where new_pairing_function is the function that performs this new pairing method. It must be a function that receives one argument and returns the list of the paired individuals. Concretely, it receives a list with the IDs of the selected individuals (see selection method), and returns a list of tuples with the IDs of the paired individuals._

_\** If more information is wanted, go to the docs and the definition of the pairing default function in ga_variable_length/tools/pairing.py._

_\*** Note that with this pairing method there can be repeated individuals and one may be paired with itself. However, this rarely happens and it can be seen as an elitism mechanism. Even more, it can be used the keep_diversity functionality to make this more improvable._



### Crossover

Suppose there are two individuals, A (individual_a) and B (individual_b), that have already been selected and paired for the mating process. The mating process that has been defined consists in taking a random number of genes of individual A, removing them from individual A and then adding them to individual B and vice versa. Thus, two new individuals that result from the combination of A and B are created. The algorithm can be defined by the next steps:

  1. First, it is calculated which genes of individual A can be copied in individual B and vice versa. If the attribute 'repeated_genes_allowed' = 1, then all the genes of individual A can be copied in individual B, as there can be repetitions, and vice versa. However, if the attribute 'repeated_genes_allowed' = 0, then it could only be 'transfered' from A to B those genes of A that are not in B, and vice versa. Lets call *genes_a_to_b* and *genes_b_to_a* to the lists of the genes that can be transferred from A to B and vice versa.

  2. It is taken a random number between 1 and the maximum number of possible genes to transfer from each individual (```len(genes_a_to_b)``` and ```len(genes_b_to_a```)). So after this random process, it is decided to get *n_a* genes from individual A to transfer them to individual B and *n_b* genes from individual B to transfer them to individual A. Notice that *n_a* may be different from *n_b*.

  3. It is tested if taking *n_a* genes from *individual_a* and transferring (copying) them to *individual_b* or if taking *n_b* genes from *individual_b* and transferring (copying) them to *individual_a* would lead to an individual with a length over or under the limits *min_number_of_genes* and *max_number_of_genes* (attributes of ```Gavl()```). If it is not valid, it is repeated **step 2** for another different and random *n_a* or *n_b*. This condition is calculated in the next way:

```python
min_number_of_genes <= (len(individual_a) - n_a + n_b) <= max_number_of_genes
min_number_of_genes <= (len(individual_b) - n_b + n_a) <= max_number_of_genes
```

  4. If the previous condition is met, among all the possible combinations of length *n_a* of the array of possible genes to transfer from individual A to individual B (*genes_a_to_b*) it is taken one (*genes_change_a*). As well, among all the possible combinations of length *n_b* of the array of possible genes to transfer from individual B to individual A (*genes_b_to_a*) it is selected another one (*genes_change_b*).

  5. It is performed the transferring of genes from individual A to individual B and vice versa, creating the two new individuals *crossed_a* and *crossed_b*. Then, with the method ```check_valid_cromosome()``` it is calculated if the new individuals are valid. If the new individual is not valid, then it is repeated **step 4** until a valid combination is found. However, if no possible valid new individuals are found among the combinations, it is repeated **step 2** for a new *n_a* and *n_b*. Notice that, in case the method ```check_valid_cromosome()``` is very strict, this step can be computationally very expensive (this is why it is recommended to give as valid all the individuals, and then giving to the invalid individuals a big penalization in the fitness function). Furthermore, if 2000 unsuccessful crossovers are tried on the same pair of individuals, it is taken as impossible to couple those individuals and their original chromosomes are returned. This condition is calculated in the next way:

```python
crossed_a = individual_a.chromosome.copy()
crossed_b = individual_b.chromosome.copy()
for gen in genes_change_a:
    crossed_a.remove(gen)
    crossed_b.append(gen)
for gen in genes_change_b:
    crossed_b.remove(gen)
    crossed_a.append(gen)
if not check_valid_individual(crossed_b): # Two loops ---> Exit the closest one (the one that selects the genes of B to be copied in A)
    break
if check_valid_individual(crossed_a):  # else try with other combination of genes of B to copy in A
    return crossed_a, crossed_b
```

  6. When a possible crossover has been found, return the new crossed individuals. 

  7. If no possible crossover is found, return the two original individuals.

_\* If other crossover method is wanted, it can be set by calling the method ```.set_hyperparameter('crossover', new_crossover_function)```, where new_crossover_function is the function that performs this new crossover method. It must be a function that receives five arguments and returns a list with the chromosomes of the newly crossed individuals. It must receive (in this order) a list ([(Individual_a, Individual_b) , ...]) with tuples of the paired individuals (of the class Individual), the minimum allowed length of the chromosome, the maximum allowed length of the chromosome, a int (boolean) that indicates if the genes can be repeated (1 = repeated genes allowed) and the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False). It must return a list with the chromosomes (list of genes that represent the individuals, NOT objects of the class Individual) of the newly created individuals._

_\** If more information is wanted, go to the docs and the definition of the crossover default function in ga_variable_length/tools/mating.py._

_\*** Note that with this crossover it is not only being tested new combinations of genes in the individual, but also new lengths of the chromosomes of these new individuals (lengths that are dependent on the lengths of the two crossed individuals)._



### Mutation

In this project it has been defined two different mutation types, explained below. As well, there are three parameters that can be configured with the method ```.set_hyperparameter()```:

* __'mutation_type'__: Type of mutation (explained below). *(default value: 'both')*
* __'mutation_rate'__: Number between 0 and 1 that represents the mutation rate. *(default value: 0.3)*
* __'max_num_gen_changed_mutation'__: It is the maximum number of genes that can be changed, ie, the maximum number of genes that are mutated when an individual is selected for mutation. *(default value: int(max_length_chromosome/3 + 1))*

It has been defined two different mutation types:

* In gene mutation (```.set_hyperparameter('mutation_type', 'mut_gen')```): It is selected randomly a number of genes of the individual between 1 and *max_num_gen_changed_mutation* and they are changed by other random genes. Consequently, it is checked if this is a valid individual (with ```check_valid_individual(chromosome)```), and if so, the mutation is done. Note that the length of the individual does not change in this kind of mutation (unless there are no more possible genes to add to this individual). 

* In length mutation (```.set_hyperparameter('mutation_type', 'addsub_gen')```): The quality of the chromosome in this genetic algorithm is also dependent on the length. Because of that, it is also allowed mutations that affect the length of the chromosome. It is either added or eliminated to the individual (add or eliminate is randomly selected for each mutated individual) a random number of genes between 1 and *max_num_gen_changed_mutation*. Consequently, it is checked if this is a valid individual, and if so, the mutation is done. Note that the new individual would be valid if its length is between *min_number_of_genes* and *max_number_of_genes* and if ```check_valid_chromosome()``` returns True.

In the current project it has been also added the possibility of selecting both mutation methods. If mutation_type = 'both' (```.set_hyperparameter('mutation_type', 'both')```) it would be selected either length mutation or gene mutation randomly for each selected individual in each generation.

_\* If other mutation method is wanted, it can be set by calling the method ```.set_hyperparameter('mutation', new_muttaion_function)```, where new_mutation_function is the function that performs this new crossover method. It must be a function that receives five arguments and returns a list with the chromosomes of the newly crossed individuals. It must receive (in this order) a list ([(Individual_a, Individual_b) , ...]) with tuples of the paired individuals (of the class Individual), the minimum allowed length of the chromosome, the maximum allowed length of the chromosome, an int (boolean) that indicates if the genes can be repeated (1 = repeated genes allowed) and the function check_valid_individual(chromosome) that receives the chromosome of an individual and returns a boolean that indicates if the individual is valid (True) or not (False). It must return a list with the chromosomes (list of genes that represent the individuals, NOT objects of the class Individual) of the newly mutated individuals._

_\** If more information is wanted, go to the docs and the definition of the mutation default function in ga_variable_length/tools/mutation.py._



### Fitness

The fitness function is one of the parameters that must be defined before the algorithm starts. Its ONLY argument is the individual’s chromosome (fitness(chromosome)), which is a list of genes, and returns the value of the fitness. It can be set by calling the method ```.set_hyperparameter('fitness', fun_fitness)```.

For example, supposing the Knapsack problem in which a bag must be filled with a collection of items in order to get the greatest value, without overpassing a concrete weight, and with possible repetition of the items, there could be a chromosome like ```['food','pen','pen','rubber']```. An example of fitness for this problem could be:

```python
MAX_WEIGTH = 15  # Maximum allowed weight
PENALIZATION = 10  # Penalization for each kg over MAX_WEIGTH

# Dictionary ---> {'item': (price, weight)}
prices_weights = {
    'pen': (5, 2),
    'pencil': (4, 2),
    'food': (7, 6),
    'rubber': (3, 1),
    'book': (10, 9)
}

def fun_fitness(chromosome):
	""" Definition of the fitness function.

    :param chromosome: (list of genes) Chromosome of the individual that is being analyzed.
    """
    fitness = 0  # Cumulative fitness
    sum_weights = 0  # Cumulative weights
    for item in chromosome:
        fitness += prices_weights[item][0]
        sum_weights += prices_weights[item][1]
    fitness -= PENALIZATION * (max(sum_weights, MAX_WEIGTH) - MAX_WEIGTH)  # Add penalization
    return fitness
```

Thus, a chromosome like the one shown in the example above, ```['food','pen','pen','rubber']```, would have a fitness value of 20.

\* _Notice that, in case it is not wanted to give a negative penalization to individuals that surpass the maximum weight, and it is wanted to avoid them, it can be done by defining a function ```check_valid_individual(chromosome)```, and setting it by calling the method ```.set_hyperparameter('check_valid_individual', check_valid_individual)```. For more information about this show the example bellow._



### Keep diversity

It may be wanted to make a great emphasis in keeping the diversity (diversity is usually one of the biggest problems of this algorithm). Because of this, it has been defined the configurable hyperparameter 'keep_diversity'. 

The algorithm to keep the diversity of the population works as follows. When an individual is repeated in the population it is substituted by a completely new and randomly generated individual. As well the worst 25% of the population is substituted by completely new and randomly generated individuals. Thus, a greater space may be search thanks to this function.

\* _This functionality is activated by calling ```.set_hyperparameter('keep_diversity', x)``` which means that every x generations the diversity techniques will be applied. For example ```.set_hyperparameter('keep_diversity', 5)``` means that every 5 generations the diversity techniques are applied. Note that its default value is -1, which means that NO diversity techniques are applied. This parameter only accepts -1 (do not apply the keep diversity methodology) and any positive number different of zero that indicates every how many generation it is applied._

\** _If used this mechanism to keep the diversity, it is recommended to do it every a considerable number of generations (5 or more). This is because the newly created individuals will probably be much worst than the rest of the population, and it would be a good practice to let them 'mate' with the rest of the individuals for some generations._

_\** If more information is wanted, go to the docs and the definition of the keep diversity function in ga_variable_length/tools/keep_diversity.py._


### Termination criteria

In this project it has been defined two different termination criteria. Stop when a maximum number of generations is reached, or stop when a goal fitness is reached. They can be set as follows:

  * __'max_num_generation_reached'__: The maximum number of generations is reached. It can be set with the method ```.set_hyperparameter('termination_criteria', {'max_num_generation_reached': n})```, where *n* is the maximum number of generations.
 
  * __'goal_fitness_reached'__: The goal fitness is reached. It can be set with the method ```.set_hyperparameter('termination_criteria', {'goal_fitness_reached': m})```, where *m* is the goal fitness.



## Example:

Following the snippets shown above, it is going to be presented an example of how to use this tool with the knapsack problem without repetition.

_Problem: Imagine you have a knapsack that can hold at much 15 kg and there is a list of items with an associated weight and value. You want to fill the knapsack with items in order to get the maximum value, but without overpassing the 15 kg of maximum weight. Suppose there is the next list of items with weights and values associated:_

|           | Value | Weight |
|-----------|-------|--------|
| Pen       | 5     | 3      |
| Pencil    | 4     | 2      |
| Food      | 7     | 6      |
| Rubber    | 3     | 1      |
| Book      | 10    | 9      |
| Scissors  | 6     | 3      |
| Glasses   | 7     | 5      |
| Case      | 7     | 7      |
| Sharpener | 2     | 1      |

_Thus, in order to solve this problem we will write the next program:_

```python
from .ga_variable_length.ga_variable_length.Gavl import Gavl # Suppose that the file is in the same folder.

ga = Gavl()  # Initialize

# Compulsory hyperparameters

ga.set_hyperparameter('size_population', 50)  # Set the size of the population to 50 individuals
ga.set_hyperparameter('min_length_chromosome', 1)  # Set the minimum length of the individual to 1 genes
ga.set_hyperparameter('max_length_chromosome',
                      9)  # Set the maximum length of the individual to 9 genes (there are 9 items)

# Definition of the fitness:
# The fitness function receives the chromosome, which is a list with the genes. It must return the fitness value (number).

MAX_WEIGTH = 15  # Maximum allowed weight
PENALIZATION = 10  # Penalization for each kg over MAX_WEIGTH

# Dictionary ---> {'item': (price, weight)}
prices_weights = {
    'pen': (5, 3),
    'pencil': (4, 2),
    'food': (7, 6),
    'rubber': (3, 1),
    'book': (10, 9),
    'scissors': (6, 3),
    'glasses': (7, 5),
    'case': (7, 7),
    'sharpener': (2, 1)
}


def fun_fitness(chromosome):
    """ Definition of the fitness function.

    :param chromosome: (list of genes) Chromosome of the individual that is being analyzed.
    """
    fitness = 0  # Cumulative fitness
    sum_weights = 0  # Cumulative weights
    for item in chromosome:
        fitness += prices_weights[item][0]
        sum_weights += prices_weights[item][1]
    fitness -= PENALIZATION * (max(sum_weights, MAX_WEIGTH) - MAX_WEIGTH)  # Add penalization
    return fitness


ga.set_hyperparameter('fitness', fun_fitness)  # Set the fitness function

# Definition of the values that the genes can take

possible_genes = ['pen', 'pencil', 'food', 'rubber', 'book', 'scissors', 'glasses', 'case', 'sharpener']
ga.set_hyperparameter('possible_genes',
                      possible_genes)  # The possible values that the genes can take are the items' names

# # In case it is wanted to NOT allow chromosomes that overpass the maximum allowed weight it can be done the following (notice that it has preferred to define a penalization term instead):
#
# def check_valid_individual(chromosome):
#     """ Definition of a valid individual.
# 
#     :param chromosome: (list of genes) Chromosome of the individual that is being analyzed.
#     """
#     sum_weights = 0  # Cumulative weights
#     for item in chromosome:
#         sum_weights += prices_weights[item][1]
#     if sum_weights > MAX_WEIGTH:
#         return False
#     else:
#         return True
#
# ga.set_hyperparameter('check_valid_individual', check_valid_individual)

# Set other parameters

ga.set_hyperparameter('termination_criteria', {
    'max_num_generation_reached': 30})  # Set the termination criteria to maximum number of generations (in this case 30 generations).
ga.set_hyperparameter('repeated_genes_allowed',
                      0)  # In this case the repeated genes are NOT allowed ---> Each item can be taken only once
ga.set_hyperparameter('minimize', 0)  # In this case the goal is to maximize the fitness
ga.set_hyperparameter('elitism_rate', 0.1)  # Set an elitism rate of 0.1
ga.set_hyperparameter('mutation_rate', 0.2)  # Set a mutation rate of 0.2
ga.set_hyperparameter('mutation_type', 'both')  # Set the mutation type
ga.set_hyperparameter('keep_diversity', 5)  # Set a mechanism to keep the diversity every 5 generations
ga.set_hyperparameter('show_progress',
                      1)  # Show progress of the genetic algorithm (this will print the generation number)

best_individual = ga.optimize()  # Optimize

best_individual, population, historic_fitness = ga.get_results()  # Get other results

```

_The best individual returned would be ```['scissors', 'pencil', 'pen', 'rubber', 'sharpener', 'glasses']```._



<span style="color:lightgray"> _Please, if any bugs are found in this project, do not hesitate to contact me and I'll fix them as soon as I can._</span>




