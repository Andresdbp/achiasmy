import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd

def simulate_multitype_branching(num_simulations, s_driver, s_del, mu, max_generations=100000, threshold=1e5):
    """
    Simulate a multitype branching process for a beneficial allele that is 
    subject to accumulating deleterious mutations.
    
    Parameters:
      - num_simulations: number of independent lineage simulations.
      - s_driver: the beneficial effect.
      - s_del: per-mutation deleterious cost (assumed multiplicative).
      - mu: probability per generation that an offspring acquires an extra deleterious mutation.
      - max_generations: maximum generations to simulate.
      - threshold: population size threshold to consider the allele as established.
      
    Returns:
      - survival_prob: estimated probability that the allele escapes early stochastic loss.
    """
    survival_count = 0

    for sim in range(num_simulations):
        lineage = {0: 1}
        extinct = False
        generation = 0

        while generation < max_generations:
            new_lineage = {}
            total_count = 0
            for k, count in lineage.items():
                fitness = (1 + s_driver) * ((1 - s_del) ** k)
                for i in range(count):
                    num_offspring = np.random.poisson(fitness)
                    for _ in range(num_offspring):
                        new_k = k + 1 if np.random.rand() < mu else k
                        new_lineage[new_k] = new_lineage.get(new_k, 0) + 1
                        total_count += 1
            generation += 1
            
            if total_count == 0:
                extinct = True
                break
            
            if total_count >= threshold:
                break
            
            lineage = new_lineage
        
        if not extinct and total_count >= threshold:
            survival_count += 1

    survival_prob = survival_count / num_simulations
    return survival_prob

# Function to run one simulation for a given pair of parameters.
def run_simulation(params):
    Gs, s_driver = params
    U = mu * Gs  # Computed if needed; not used in this simulation function.
    survival_prob = simulate_multitype_branching(num_simulations, s_driver, s_del, mu)
    return (Gs, s_driver, survival_prob)

# Simulation grid parameters
num_simulations = 10000  # number of simulations for each parameter combination
s_del = 0.030           # deleterious cost per mutation
mu = 1e-8               # mutation probability per generation

# Parameter ranges
gs_min = 2e4
gs_max = 2e6
s_driver_min = 0.0
s_driver_max = 0.5

# Adjustable steps (change these as desired)
gs_steps = 100
s_driver_steps = 100

gs_values = np.linspace(gs_min, gs_max, gs_steps)
s_driver_values = np.linspace(s_driver_min, s_driver_max, s_driver_steps)

if __name__ == '__main__':
    # Prepare the list of parameter combinations.
    param_list = [(Gs, s_driver) for Gs in gs_values for s_driver in s_driver_values]

    print("Running simulations in parallel using multiprocessing (4 cores)...")
    # Create a pool with 4 explicitly defined processes.
    with mp.Pool(processes=4) as pool:
        results_list = pool.map(run_simulation, param_list)

    # Convert results_list into a matrix.
    results = np.empty((gs_steps, s_driver_steps))
    for res in results_list:
        Gs, s_driver, survival_prob = res
        i = np.argmin(np.abs(gs_values - Gs))
        j = np.argmin(np.abs(s_driver_values - s_driver))
        results[i, j] = survival_prob
        print(f"Gs: {Gs:.1e}, s_driver: {s_driver:.3f}, survival_prob: {survival_prob:.4f}")

    # Save the results matrix to a CSV file
    df_results = pd.DataFrame(results, index=gs_values, columns=s_driver_values)
    df_results.to_csv("results.csv")
    print("Results saved to results.csv")

    # Plot the heatmap.
    plt.figure(figsize=(8, 6))
    cax = plt.imshow(results, origin='lower', aspect='auto', 
                     extent=[s_driver_min, s_driver_max, gs_min, gs_max])
    plt.colorbar(cax, label='Survival Probability')
    plt.xlabel('s_driver')
    plt.ylabel('Gs')
    plt.title('Heatmap of Survival Probability across Parameter Space')
    plt.show()

