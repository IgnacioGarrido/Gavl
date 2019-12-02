# GA_chromosomeWithVariableLength

This repository contains a python program that performs a genetic algorithm in which the chromosome has variable length. The algorithm receives a table in which each column represents an entity to be chosen (eg. an enterprise), and each row represents an item to be selected (eg. a product), with its corresponding weight.

## How to use:

  1. Download the files:
```shell
	git clone https://github.com/IgnacioGarrido/GA_chromosomeWithVariableLength.git
```
  2. Import the GA_vl function from the file GA_main.py. Note that this is the main and only function needed to execute the genetic algorithm. All the hyperparameters should be passed to this function.
```python
	from path.GA_main import GA_vl
```
  3. Call the GA:
```python
	best_chrom, historic_fitness, pop = GA_vl(num_individuals = 100, df = mast_np, min_number_of_genes = 3, max_number_of_genes = 5, PENALIZATION = 7000)
```
Where:
- best_chrom: Chromosome of the best individual.
- historic_fitness: Array containing the best fitness of each generation.
- pop: List containing the whole population of the last generation (Individuals with the attributes chromosome, fitness and normalized fitness).

4. If needed, fork the repository and make the all the modifications you need for your project. This is open source software and any additional adaptations are welcome :).

## Parameters that can be passed to the function ```GA_vl()```

* num_individuals: Number of individuals per generation.
* df: numpy/pandas array in which each column represents an entity to be chosen (eg. an enterprise) and each row an item in which that entity participates. If it is a pandas df, then the colnames are taken as the entities (for presenting in a nice way the results). Note that the chromosome consists in a selection of columns, ie, a good selection of enterprises (or whatever each column represents) according to the defined fitness. For example, the chromosome [1,5,8] represents that the enterprises (or whatever each column represents) in columns 1, 5 and 8 have been chosen.
* min_number_of_genes: Minimum number of genes of a chromosome (minimum allowed length of a chromosome). *(default value: 3)*
* max_number_of_genes: Maximum number of genes of a chromosome (maximum allowed length of a chromosome). *(default value: 5)*
* max_num_gen_changed_crossover: *(See crossover section below)* Maximum number of genes that can be taken of an individual in order to make the crossover. *(default value: 2)*
* termination_criteria: Termination criteria's function. It can take the next parameters:
	- 'max_num_generation_reached' -> Max num generations reached. *(default value)*
	- 'goal_fitness_reached' -> Minimum fitness reached. 
* elitism_rate: percentage of elitism. *(default value: 0.05)*
* mutation_rate: Mutation rate. *(default value: 0.4)*
* mutation_type: *(See mutation section below)* It represents the mutation type. It can take the next parameters:
	- 'mut_gen': One of the genes of an individual is randomly changed.
	- 'addsub_gen': It is added or eliminated (randomly) a new element to an Individual.
	- 'both': 'mut_gen' and 'addsub_gen' are randomly applied. *(default value)*
* num_gen_changed_mutation: *(See mutation section below)* It is the number of genes that is changed, ie, the number of genes that are added, substracted or changed by new ones, when an individual is selected for mutation. *(default value: 1)*
* max_gen: *(if termination criteria = 'max_num_generation_reached')* Specification of the maximum number of generations. *(default value: 100)*
* goal_fitness: *(if termination criteria = 'goal_fitness_reached')* Specification of the goal fitness. *(default value: None)*
* keep_diversity_5_gen: *(See keep diversity section below)* If 1, then the function keep_diversity is called every 5 generations. In this function it is checked if there are any repeated individuals, and if so, these repetitions are substituted by a completely and randomly generated new individual. *(default value: 1)*
* min_item_per_row: *(used to calculate the fitness -> See in the fitness section below)* If 1, it is calculated the minimum combination of items of the passed dataframe (one item per row) for a given individual. If it is 0 the maximum is calculated. For example, if *min_item_per_row = 1*, for the individual [1,5,8] it would be calculated the minimum combination of items (one item per row) for the combination of the columns 1, 5 and 8 of the passed dataframe *df*. *(default value: 1)*
* minimize: *(used to calculate the fitness -> See in the fitness section below)* If 1, the goal would be to minimize the fitness, and if 0, the goal would be maximizing the fitness. If 1 the fitness value is inverse normalized, ie, higher value mappped to 0 and lower to 1. *(default value: 1)*
* MAX_NUM_TRANS: *(used to calculate the fitness -> See in the fitness section below)* Maximum length of the chromosome up to which there is no penalization. *(default value: 3)*
* PENALIZATION: *(used to calculate the fitness -> See in the fitness section below)* Penalization value per each new extra item (extra with respect to *MAX_NUM_TRANS*) in the chromosome (or for their rating). *(default value: 0)*
* PERCENT: *(used to calculate the fitness -> See in the fitness section below)* Percentage of the total cost that penalizes each extra item (extra with respect to *MAX_NUM_TRANS*). *(default value: 0)*
* RATING_TRANS: *(used to calculate the fitness -> See in the fitness section below)* List with the ordered rating of each item represented in the columns. *(default value: [])*

## The algorithm

Genetic algorithms (GA) are a well-known optimization tool used to solve problems in which the usage of the brute force may lead to inadmisible execution time. GA use nature-based heuristics that give a great approximation to the optimal solution and, luckily, converge to the best solution. 

Concretely, this project tries to solve problems in which the chromosome of the GA has variable length.

### Problems that can be solved with this GA

This optimization tool can be used to solve problems in which the chromosome has variable length. Concretely, the tool has been optimized to solve problems in which it is given a database being its columns the elements that have to be selected, and it has to be selected one item of each row. Consecuently, the data in the given table, the data represents the weight of a item (represented in a row) for a concrete element that has to be selected (row). 

For example, suppose that it is given the next table in which it is represented the prices offered by six distributors to the food items you need. You would need to get all of the food items for the best price, but not all the distributors offer all the products (NaN values). As well, you probably would like to make a selection of distributors based in some heuristics you choose, such as the number of distributors (you may prefer a smaller number of distributors, even though the price may slightly increase) or the rating of each distributor (heuristics explained in the fitness section below). 


|         | Distributor 1 | Distributor 2 | Distributor 3 | Distributor 4 | Distributor 5 | Distributor 6 |
|---------|---------------|---------------|---------------|---------------|---------------|---------------|
| Apples  | 4.00          | 4.55          | 3.85          | NaN           | 4.05          | NaN           |
| Bread   | 2.05          | 1.99          | NaN           | 1.70          | 2.5           | 1.90          |
| Tomato  | 5.15          | 4.85          | NaN           | 5.00          | NaN           | 4.95          |
| Onions  | 2.50          | 2.45          | 2.75          | 2.60          | 2.60          | 2.80          |
| Bananas | 3.50          | NaN           | 3.55          | 3.40          | 3.70          | 3.75          |
| Pie     | NaN           | 10.50         | NaN           | 11.00         | 10.50         | NaN           |


The GA presented in this project optimizes this selection of individuals based in the heuristics you select. 

### Individual

The individuals are represented by the class Individual with the following definition:

```python
class Individual:
  def __init__(self, chromosome, fitness, normalized_fitness):
    self.chromosome = chromosome
    self.fitness = fitness
    self.normalized_fitness = normalized_fitness
```

Each individual has three attributes:
	- Chromosome: Is a variable length list (length between the limits [*min_number_of_genes, max_number_of_genes*] passed to the function ```GA_vl()```). Each element of the list represents an entity to be selected. For example, given the example of the previous table, the chromosome [1,5,6] would represent that it has been selected the distributors 1, 5 and 6.
	- Fitness: Fitness value according to the defined heuristics (see fitness section below).
	- Normalized_fitness: Normalized fitness. If the value of the parameter *minimize* passed to the function ```GA_vl()``` is 1, the fitness is inverse normalized (the minimum value of fitness of the population is normalized to 1, and the biggest is normalized to 0) and viceversa.

Of course, the individuals are randomly generated in the first generation.

### Selection - Roulette wheel

The only performed selection method in this project is roulette wheel selection. 

If it is wanted to program any other selection method, it should be changed the function ```roulette_selection()``` defined in the file *GA_functions/selection.py*. If other name is given to this new selection function, it should be changed its calling inside the function ```next_generation()``` in the file *GA_function/next_generation.py*. Note that if a new function of selection of individuals is made, it should return a list with the IDs of the new individuals. For example, it could return a list like [4,60,36,44,15,22] in which it is represented that the individuals 4, 60, 36, 44, 15 and 22 have been selected.

### Pairing - Random pairing

The only performed pairing method is random pairing given a selection. 

If it is wanted to program any other pairing method, it should be changed the function ```pairing()``` defined in the file *GA_functions/pairing.py*. If other name is given to this new pairing function, it should be changed its calling inside the function ```next_generation()``` in the file *GA_function/next_generation.py*. Note that if a new function of pairing of individuals is made, it should return a list with the paired individuals in tuples. For example, given a selection of individuals given by [4,60,36,44,15,22] this pairing method could return a pairing like [(4,36),(15,22),(60,44)].

### Fitness

### Crossover

### Mutation

### Keep diversity

It may be wanted to make a great emphasis in keeping the diversity. For so, it has been defined the function ```keep_diversity()```in the file *G*

keep_diversity_5_gen: *(See keep diversity section below)* If 1, then the function keep_diversity is called every 5 generations. In this function it is checked if there are any repeated individuals, and if so, these repetitions are substituted by a completely and randomly generated new individual. *(default value: 1)*


### Functions that may be needed to be defined again if some part of the project is changed

- check_valid_cromosome (in *GA_functions/aux_functions.py*): It receives a chromosome and the dataframe. This function defines what a valid chromosome is. If the passed chromosome is valid, then it returns *True* and if it is invalid it returns *False*. In the current definition of this function it is defined as a valid individual one that covers with at least one possible selection of an item for all of the rows. For example, passed the dataframe *df*, the individual [1,5,8] would be valid if there is NO row of ```df[:,[1,5,8]]``` that has all its items equal to NaN.


