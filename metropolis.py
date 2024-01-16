# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 13:55:37 2024

@author: ljh3218
"""

import numpy as np
import pandas as pd
import matlab.engine as matlab
import uuid

#import multiprocessing
#import matplotlib.pyplot as plt

# Define the function that runs the IonMonger simulation and returns the outputs
def run_simulation(eng, inputs):
    # Convert numpy array into python list (matlab doesn't accept numpy arrays)
    params_list = []
    for i in range(len(inputs)):
        params_list.append(10**float(inputs[i]))

    # Call the master function in MATLAB using the MATLAB engine and calculate outputs
    #Set the bias voltages and pulse voltages here 
    #NB: remeber that theese need to be lists! Also, all need to be the same data type 
    #so use 0.0, not just 0
    Vbias = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    Vpulse = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
              0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2]
    
    sol = eng.SaP_test_params(params_list, Vbias, Vpulse)
    sol = np.asarray(sol)
    sol = sol.flatten()

    return sol

# Define the Metropolis Hastings algorithm
def metropolis_hastings(eng, initial_state, num_samples):
    # Initialize the current state and likelihood
    current_state = initial_state
    current_log_posterior = log_posterior(eng, current_state)
    
    # Initialize the samples and acceptance rate
    samples = [current_state]
    acceptance_rate = 0.0
    
    # Loop over the desired number of samples
    for i in range(num_samples):
        # Propose a new state using the proposal distribution
        proposed_state = proposal_distribution(current_state)
        while log_prior(proposed_state) == -np.inf:
            proposed_state = proposal_distribution(current_state)
        
        # Calculate the posterior of the proposed state
        proposed_log_posterior = log_posterior(eng, proposed_state)
        
        # Debugging print
        # print(f"Proposed state: {proposed_state}")
        print(f"Posterior: {np.exp(proposed_log_posterior)}")

        # Calculate the acceptance ratio
        acceptance_ratio = min(1, np.exp(proposed_log_posterior - current_log_posterior))
        
        # Accept or reject the proposed state
        if np.random.uniform() < acceptance_ratio:
            current_state = proposed_state
            current_log_posterior = proposed_log_posterior
            acceptance_rate += 1.0
        
        # Add the current state to the samples
        samples.append(current_state)
        print(i)
    
    # Return the samples and acceptance rate
    return samples, acceptance_rate / num_samples

# Scale the outputs so they are same order of magnitude
def scale_outputs(outputs):
    for i in range(len(y)):
        outputs[i] = outputs[i] / y[i]
    return outputs

# Randomly sample inputs from uniform priors
def initial_sample():
    size = len(prior_ranges)
    initial_inputs = np.zeros(size)
    for i in range(size):
        initial_inputs[i] = np.random.uniform(prior_ranges[i, 0], prior_ranges[i, 1])
    return initial_inputs

# Define the log-prior distribution, which is a uniform distribution over the prior ranges
def log_prior(inputs):
    for i in range(len(inputs)):
        if (inputs[i] < prior_ranges[i][0]) or (inputs[i] > prior_ranges[i][1]):
            return -np.inf
    return 0.0

# Define the log-likelihood function using the run_simulation and likelihood functions
def log_likelihood(eng, inputs):
    outputs = run_simulation(eng, inputs)
    outputs = scale_outputs(outputs)

    # Calculate the mean squared error between the logged outputs and experimental logged outputs y
    mse = np.mean((outputs - np.ones((len(y),)))**2)
    # Return the log_likelihood, which is proportional to -0.5 * mse for a normal distribution
    return (-0.5 * mse/0.05)

# Log posterior is just the log prior + log likelihood
def log_posterior(eng, inputs):
    return log_prior(inputs) + log_likelihood(eng, inputs)

# Define the proposal distribution, which is a normal distribution centered on the current state
def proposal_distribution(current_state):
    return np.random.normal(current_state, jump_dist_sigmas*np.abs(current_state))

def run_single_chain(n_iter):
    ranges = np.asarray([[1e16, 1e19],      # mobile ion vacancy density
                         [0.01, 0.5],       # HTL Energy Offset
                         [0.01, 0.5],       # ETL Energy Offset 
                         [0.075, 0.5],      # HTL Fermi Level Offset 
                         [0.075, 0.5],      # ETL Fermi Level Offset
                         [5e-9, 1e-6],      # electron bulk lifetime
                         [5e-9, 1e-6],      # hole bulk lifetime
                         [0.1, 100],        # v_surf at HTL/pero interface
                         [0.1, 100],        # v_surf at ETL/pero interface 
                         [0.1, 10],         # e- mobility in perovskite
                         [0.1, 10],         # h+ mobility in perovskite 
                         [1e-4, 0.1],       # e- mobility in ETL 
                         [1e-4, 0.1],       # h+ mobility in HTL
                         [10, 80]])         # permittivity of ETL 
    
    global prior_ranges, jump_dist_sigmas, y
    prior_ranges = np.log10(ranges)
    
    y = get_input_data()
    
    jump_dist_sigmas = 0.01*np.ones(len(prior_ranges[:,0]))
    
    np.random.seed()
    initial_inputs = initial_sample()
    # Start the MATLAB engine
    try:
        eng = matlab.start_matlab("-nosplash -nodisplay")
        eng.initialise_df(nargout=0)
        eng.SetUpParallelPool(nargout=0)

        # Run the Metropolis Hastings algorithm to sample from the likelihood distribution
        samples, acceptance_rate = metropolis_hastings(eng, initial_inputs, num_samples=n_iter)

        # Print the acceptance rate and the mean and standard deviation of the samples
        print(f"Acceptance rate: {acceptance_rate:.2f}")
        print(f"Mean: {np.mean(samples, axis=0)}")
        print(f"Standard deviation: {np.std(samples, axis=0)}")

        # Stop the MATLAB engine
        eng.quit()

        return np.asarray(samples)
    except:
        print('matlab engine error')
        
def get_input_data():
    data = pd.read_excel('input_data.xlsx', header = None)
    data_ar = data.to_numpy()
    flat_data_ar = data_ar.flatten()
    return flat_data_ar

if __name__ == "__main__":
    n_iter = 100
    result = run_single_chain(n_iter)   
    df = pd.DataFrame(result)
    df.to_csv('{}.csv'.format(uuid.uuid4().hex))
