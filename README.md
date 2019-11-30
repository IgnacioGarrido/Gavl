# GA_chromosomeWithVariableLength

This repository contains the files needed to perform a genetic algorithm in which the chromosome has variable length. The algorithm receives a table in which each column represents an entity to be chosen (eg. an enterprise), and each row represents an item to be selected (eg. a product), with its corresponding weight.

How to use:
    1. Download the files:
    
    ```shell
        git clone https://github.com/IgnacioGarrido/GA_chromosomeWithVariableLength.git
    ```
    
    2. Import the GA_vl function from the file GA_main.py. Note that this is the main and only function needed to execute the genetic algorithm. All the hyperparameters should be passed to this function.
    
    ```python
        from path.GA_main import GA_vl
    ```
    
    3. Call:
    
    ```python
      best_chrom, historic_fitness, pop = GA_vl(num_individuals = 100, df = mast_np, min_number_of_genes = 3, max_number_of_genes = 5, PENALIZATION = 7000)
    ```
    
    
      
## The algorithm:

### 
        
