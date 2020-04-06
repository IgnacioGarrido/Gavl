"""
In this file it is defined a function to keep the diversity.

Functions:
    keep_diversity: Function called to keep the diversity.
"""


def keep_diversity(population, generate_chromosome, min_length_chromosome, max_length_chromosome, possible_genes,
                   repeated_genes_allowed, check_valid_chromosome):
    """ This function is called when it is wanted to do a great emphasis in the diversity of the population. When an individual is repeated in the population it is substituted by a completely new and randomly generated individual. As well the worst 25% of the population is substituted by completely new and random individuals. Note that before calling this function the population MUST be sorted by fitness.

    :param population: (list of Individuals) Sorted population by fitness (from best to worst).
    :param generate_chromosome: (function) Function to generate a new chromosome.
    :param min_length_chromosome: (int) Minimum allowed length of the chromosome.
    :param max_length_chromosome: (int) Maximum allowed length of the chromosome.
    :param possible_genes: (list of ...) List with the all the possible values that the genes can take.
    :param repeated_genes_allowed: (bool) It is a boolean that indicates whether the genes can be repeated in the chromosome (repeated_genes_allowed = 1) or they cannot be repeated (repeated_genes_allowed = 0).
    :param check_valid_chromosome: (function) Function that receives the chromosome and returns True if it creates a valid individual and False otherwise.
    :return:
        (list of chromosomes) List of chromosomes that will represent the next generation.
    """
    new_population = [population[0].chromosome]  # List that will hold the new population
    list_chromosomes = [ind.chromosome for ind in population]  # List with all the chromosomes
    copy_list_chromosomes = list_chromosomes.copy()
    try:
        for i in range(len(list_chromosomes)):
            chrom = list_chromosomes[i]
            if len(chrom):
                if type(chrom[0]) == int or type(chrom[0]) == float or type(chrom[0]) == str or type(chrom[0]) == chr:
                    chrom.sort()  # So repeated individuals are correctly eliminated
                elif type(chrom[0]) == dict:
                    list_chromosomes[i] = sorted(chrom, key=lambda i: (list(i.keys())[0], list(i.values())[
                        0]))  # Not 100% sure it will work and sort uniquely the individuals, but it may do the trick
    except:
        print('Sorting error in keep_diversity.')
        list_chromosomes = copy_list_chromosomes.copy()  # If some error in the sorting, just copy the original list_chromosomes
    for i in range(1, int(len(list_chromosomes) / 4)):
        if population[i].chromosome in list_chromosomes[:i]:
            new_ind = generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes,
                                          repeated_genes_allowed)  # Generate a new individual
            while not check_valid_chromosome(new_ind):
                new_ind = generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes,
                                              repeated_genes_allowed)  # Check if the new individual is valid
            new_population.append(new_ind)
        else:
            new_population.append(population[i].chromosome)
    while len(new_population) < len(
            population):  # Add new randomly generated individuals in substitution to the worst 25%.
        new_ind = generate_chromosome(min_length_chromosome, max_length_chromosome, possible_genes,
                                      repeated_genes_allowed)  # Generate a new individual
        if check_valid_chromosome(new_ind):
            new_population.append(new_ind)
    return new_population
