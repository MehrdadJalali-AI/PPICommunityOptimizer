"""
Lotus Effect Algorithm (LEA) - Preserved from MultiModelLEA.ipynb
This implementation maintains the original LEA mechanics as-is.
Only the fitness function and decision variables are adapted for community detection.
"""

import numpy as np
from scipy.special import gamma
import logging

logger = logging.getLogger(__name__)


class LotusEffectAlgorithm:
    """
    Lotus Effect Algorithm for optimization.
    Preserves original mechanics: mutation, self-cleaning, selection.
    """
    
    def __init__(self, population_size, dimensions, lower_bound, upper_bound, 
                 max_function_evaluations, fitness_function, random_seed=None):
        """
        Initialize LEA.
        
        Args:
            population_size: Number of individuals in population
            dimensions: Number of decision variables
            lower_bound: Lower bound for each dimension (scalar or array)
            upper_bound: Upper bound for each dimension (scalar or array)
            max_function_evaluations: Maximum number of fitness evaluations
            fitness_function: Function to minimize (takes solution vector, returns scalar)
            random_seed: Random seed for reproducibility
        """
        if random_seed is not None:
            np.random.seed(random_seed)
            
        self.population_size = population_size
        self.dimensions = dimensions
        self.lower_bound = np.array(lower_bound) if isinstance(lower_bound, (list, np.ndarray)) else lower_bound
        self.upper_bound = np.array(upper_bound) if isinstance(upper_bound, (list, np.ndarray)) else upper_bound
        self.max_function_evaluations = max_function_evaluations
        self.fitness_function = fitness_function
        
        # Initialize population
        if isinstance(self.lower_bound, np.ndarray) and len(self.lower_bound) == dimensions:
            self.population = np.random.uniform(
                self.lower_bound, self.upper_bound, (population_size, dimensions)
            )
        else:
            self.population = np.random.uniform(
                self.lower_bound, self.upper_bound, (population_size, dimensions)
            )
        
        self.best_solution = self.population[0].copy()
        self.best_fitness = float('inf')
        self.function_evaluations = 0
        self.alpha = 0.5  # Adaptive step-size
        self.beta = 1.5  # Levy flight parameter
        self.memory = []  # Store previous best solutions
        self.history = []  # Store best fitness per iteration
        
    def levy_flight(self):
        """
        Generate Levy flight step for exploration.
        """
        sigma = (gamma(1 + self.beta) * np.sin(np.pi * self.beta / 2) /
                 (gamma((1 + self.beta) / 2) * self.beta * 2 ** ((self.beta - 1) / 2))) ** (1 / self.beta)
        u = np.random.normal(0, sigma, size=self.dimensions)
        v = np.random.normal(0, 1, size=self.dimensions)
        step = u / np.abs(v) ** (1 / self.beta)
        return 0.01 * step
    
    def update_positions(self):
        """
        Update positions of all individuals using LEA mechanics.
        """
        for i in range(self.population_size):
            if self.function_evaluations >= self.max_function_evaluations:
                return
            step = self.alpha * (self.best_solution - self.population[i]) + self.levy_flight()
            self.population[i] += step
            
            # Clip to bounds
            if isinstance(self.lower_bound, np.ndarray):
                self.population[i] = np.clip(self.population[i], self.lower_bound, self.upper_bound)
            else:
                self.population[i] = np.clip(self.population[i], self.lower_bound, self.upper_bound)
            
            self.function_evaluations += 1
    
    def evaluate_fitness(self):
        """
        Evaluate fitness for all individuals and update best solution.
        """
        for i in range(self.population_size):
            if self.function_evaluations >= self.max_function_evaluations:
                return
            fitness = self.fitness_function(self.population[i])
            self.function_evaluations += 1
            
            if fitness < self.best_fitness:
                self.best_fitness = fitness
                self.best_solution = self.population[i].copy()
                self.memory.append(self.best_solution)
        
        self.history.append(self.best_fitness)
    
    def optimize(self):
        """
        Run LEA optimization.
        
        Returns:
            best_solution: Best solution found
            best_fitness: Best fitness value
            function_evaluations: Total function evaluations used
        """
        while self.function_evaluations < self.max_function_evaluations:
            self.update_positions()
            self.evaluate_fitness()
            
            if len(self.history) % 100 == 0:
                logger.debug(f"Evaluations: {self.function_evaluations}, Best Fitness: {self.best_fitness:.6f}")
        
        return self.best_solution, self.best_fitness, self.function_evaluations

