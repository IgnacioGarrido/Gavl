# GA_chromosomeWithVariableLength

This repository contains the files needed to perform a genetic algorithm in which the chromosome has a variable length. The algorithm receives a table in which each column represents an entity to be chosen (eg. an enterprise), and each row represents an item to be selected, with its corresponding weight.

How to use:
  1. Download files
  2. Import GA_main
  3. Call:     
      best_chrom, historic_fitness, pop = GA_vl(num_individuals = 100, df = mast_np, min_number_of_genes = 3, max_number_of_genes = 5, PENALIZATION = 7000)
