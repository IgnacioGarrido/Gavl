"""
In this file it is defined the function to perform mutation.

Function: 
    mutation: Function that performs mutation.
    mutate_genes_manner: Auxiliary function to make the mutation by changing the genes.
    mutate_length_manner: Auxiliary function to make the mutation in length.
"""
import random
from .aux_functions.combinations import combinations


def mutation(chromosomes_to_mutate, mutation_type, max_num_gen_changed_mutation, min_length_chromosome,
             max_length_chromosome, repeated_genes_allowed, check_valid_individual, possible_genes):
    """ This function receives a chromosome and performs a random mutation over one of its elements. It iterates over all the possible mutations until one is found, moment in which the execution stops. If no mutation is found, then it returns the input chromosome. Notice that the function check_valid_individual is used to test the created individuals, and if it is done 1000 unsuccessful mutations on the same individual, it is taken as an impossible to mutate individual and its original chromosome is returned.

    :param chromosomes_to_mutate: (list of chromosomes) List of the chromosomes that are going to be mutated.
    :param mutation_type: (str) String that represents the mutation type. It can ONLY take the values ‘mut_gene’, ‘addsub_gene’ or ‘both’. It is the attribute .mutation_type of the class Gavl().
    :param max_num_gen_changed_mutation: (int) Maximum number of genes changed in the mutation. It is the attribute .max_num_gen_changed_mutation of the class Gavl().
    :param min_length_chromosome: (int) Minimum allowed number of genes of a chromosome. It is the attribute .min_length_chromosome of the class Gavl().
    :param max_length_chromosome: (int) Maximum allowed number of genes of a chromosome. It is the attribute .max_length_chromosome of the class Gavl().
    :param repeated_genes_allowed: (int) Boolean that indicates if an individual can have repeated genes. It can take the values 1 (repeated genes allowed) or 0 (repeated genes not allowed).
    :param check_valid_individual: (function) Function that receives a chromosome and returns a boolean that indicates whether the chromosome makes a valid individual (True) or not (False) according to some criteria.
    :param possible_genes: (list of genes) List with the possible genes.
    :return:
        * :list_new_mutated_chromosomes: (list of chromosomes) List with the chromosomes of the individuals that are mutated.
    """
    if mutation_type not in ['mut_gene', 'addsub_gene', 'both']:
        raise ValueError("The parameter 'mutation_type' can only take the values 'mut_gene', 'addsub_gene' or 'both'.")
    list_new_mutated_chromosomes = []  # List that will contain the chromosomes of the mutated individuals.
    for chromosome in chromosomes_to_mutate:
        if repeated_genes_allowed:
            mutation_genes = possible_genes  # List with the possible genes that can be mutated (genes that are not in the individual).
        else:
            mutation_genes = [e for e in possible_genes if
                              e not in chromosome]  # List with the possible genes that can be mutated (genes that are not in the individual).
        both_mutations_selection = int(
            random.random() > 0.5)  # If mutation_type = 'both' it is chosen randomly one of both mutation types.
        # Start the algorithm:
        if (mutation_type == 'mut_gene') or ((mutation_type == 'both') and (
                both_mutations_selection == 0)):  # Mutate the genes, ie, change some genes for others.
            new_mutated_chromosome = mutate_genes_manner(chromosome, max_num_gen_changed_mutation, mutation_genes,
                                                         check_valid_individual)  # Mutate the individual
            list_new_mutated_chromosomes.append(new_mutated_chromosome)  # Add the mutated individual to the population
        elif (mutation_type == 'addsub_gene') or ((mutation_type == 'both') and (both_mutations_selection == 1)):
            new_mutated_chromosome = mutate_length_manner(chromosome, max_num_gen_changed_mutation, min_length_chromosome,
                                                          max_length_chromosome, mutation_genes,
                                                          check_valid_individual)  # Mutate the individual
            list_new_mutated_chromosomes.append(new_mutated_chromosome)  # Add the mutated individual to the population
    return list_new_mutated_chromosomes


def mutate_genes_manner(chromosome, max_num_gen_changed_mutation, mutation_genes, check_valid_individual):
    """ This function performs the mutation of the genes (not the length).

    :param chromosome: (list of genes) Chromosome to mutate.
    :param max_num_gen_changed_mutation: (int) Maximum number of genes changed in the mutation. It is the attribute .max_num_gen_changed_mutation of the class Gavl().
    :param mutation_genes: (list of genes) List with the possible genes that can be mutated (genes that are not in the individual).
    :param check_valid_individual: (function) Function that receives a chromosome and returns a boolean that indicates whether the chromosome makes a valid individual (True) or not (False) according to some criteria.
    :return:
        * :new_chromosome: (list of genes) New mutated chromosome.
    """
    random.shuffle(chromosome)  # Randomize
    random.shuffle(mutation_genes)  # Randomize
    # Get randomly the number of genes to change
    num_genes_to_mutate = list(range(1, min(max_num_gen_changed_mutation, len(
        mutation_genes), len(chromosome)) + 1))  # List from which it is selected the number of genes to mutate
    random.shuffle(num_genes_to_mutate)  # Randomize
    count_mutations_tried = 0  # Count of the number of mutations tested ---> if count_mutations_tried == 1000 then break
    # Start the algorithm
    for num_gen in num_genes_to_mutate:  # Take a random number of genes to mutate
        genes_in_combinations = combinations(mutation_genes, num_gen)  # Combinations of the genes that will be added
        for gen_in_comb in genes_in_combinations:
            genes_out_combinations = combinations(chromosome,
                                                  num_gen)  # Combinations of the genes that will be taken from the individual
            for gen_out_comb in genes_out_combinations:
                count_mutations_tried += 1
                new_chromosome = chromosome.copy()  # New mutated chromosome ---> It is mutated below
                for gen in gen_out_comb:
                    new_chromosome.remove(gen)
                new_chromosome.extend(gen_in_comb)
                if check_valid_individual(new_chromosome):  # Else try other combination of genes
                    return new_chromosome
                if count_mutations_tried >= 1000:  # Try 1000 combinations or break
                    return chromosome
    return chromosome  # If not valid mutation is found, return the original chromosome


def mutate_length_manner(chromosome, max_num_gen_changed_mutation, min_length_chromosome, max_length_chromosome,
                         mutation_genes, check_valid_individual):
    """ This function performs the mutation of the genes (not the length).

    :param chromosome: (list of genes) Chromosome to mutate.
    :param max_num_gen_changed_mutation: (int) Maximum number of genes changed in the mutation. It is the attribute .max_num_gen_changed_mutation of the class Gavl().
    :param min_length_chromosome: (int) Minimum allowed number of genes of a chromosome. It is the attribute .min_length_chromosome of the class Gavl().
    :param max_length_chromosome: (int) Maximum allowed number of genes of a chromosome. It is the attribute .max_length_chromosome of the class Gavl().
    :param mutation_genes: (list of genes) List with the possible genes that can be mutated (genes that are not in the individual).
    :param check_valid_individual: (function) Function that receives a chromosome and returns a boolean that indicates whether the chromosome makes a valid individual (True) or not (False) according to some criteria.
    :return:
        * :new_chromosome: (list of genes) New mutated chromosome.
    """
    random.shuffle(chromosome)  # Randomize
    random.shuffle(mutation_genes)  # Randomize
    # Choose if we are adding genes or not
    if len(chromosome) == min_length_chromosome:  # Add
        add = 1
    elif len(chromosome) == max_length_chromosome:  # Sub
        add = 0
    else:
        add = int(random.random() > 0.5)  # Randomly choose to add or sub a gen
    # Get randomly the number of genes to add/sub
    if add:  # Add
        num_genes_to_mutate = list(range(1, min(max_num_gen_changed_mutation, len(
            mutation_genes), max_length_chromosome - len(
            chromosome)) + 1))  # List from which it is selected the number of genes to mutate
    else:  # Sub
        num_genes_to_mutate = list(range(1, min(max_num_gen_changed_mutation, len(
            chromosome) - min_length_chromosome) + 1))  # List from which it is selected the number of genes to mutate
    random.shuffle(num_genes_to_mutate)
    count_mutations_tried = 0  # Count of the number of mutations tested ---> if count_mutations_tried == 1000 then break
    # Start the algorithm
    if add:  # Add new genes
        for num_gen in num_genes_to_mutate:  # Take a random number of genes to mutate
            genes_in_combinations = combinations(mutation_genes,
                                                 num_gen)  # Combinations of the genes that will be added
            for gen_comb in genes_in_combinations:
                count_mutations_tried += 1
                new_chromosome = chromosome.copy()
                new_chromosome.extend(gen_comb)
                if check_valid_individual(new_chromosome):  # Else try other combination of genes
                    return new_chromosome
                if count_mutations_tried >= 1000:  # Try 1000 combinations or break
                    return chromosome
    else:  # Eliminate an item
        for num_gen in num_genes_to_mutate:  # Take a random number of genes to mutate
            genes_in_combinations = combinations(chromosome,
                                                 num_gen)  # Combinations of the genes that will be eliminated
            for gen_comb in genes_in_combinations:
                count_mutations_tried += 1
                new_chromosome = chromosome.copy()
                for gen in gen_comb:
                    new_chromosome.remove(gen)
                if check_valid_individual(new_chromosome):  # Else try other combination of genes
                    return new_chromosome
                if count_mutations_tried >= 1000:  # Try 1000 combinations or break
                    return chromosome
    return chromosome  # If not valid mutation is found, return the original chromosome
