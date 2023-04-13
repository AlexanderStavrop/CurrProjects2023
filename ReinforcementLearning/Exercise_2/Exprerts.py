import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


######################################################## Tester ########################################################
class Tester:    
    def __init__(self, numOfExperts=None, numOfRounds=None, dataPath=None, fname=None):
        # Assigning the number of bandits and the rounds to the class
        self.dataPath = dataPath
        self.fname = fname

        # Extracting the data from the csv file
        dataHandler = DataHandler(self.dataPath, self.fname)
        self.K, self.T, self.data  = dataHandler.extract_data()

        # If the user has given the number of experts we use the value extracted from the csv file
        if numOfExperts is not None and numOfExperts <= self.K:
            self.K = numOfExperts
        
        # If the user has given the number of rounds we use the value extracted from the csv file
        if numOfRounds is not None and numOfRounds <= self.T:
            self.T = numOfRounds


    # Function for running the experiment for the experts algorithm
    def run_experiment_BExperts(self, type=None):
        # Creating an instance of the Experts class
        bexperts = BExperts(self.K, self.T, self.data, type)
        
        # Running the experiment
        cumulative_regret = bexperts.run_experiment()

        # Returning the cumulative regret
        return cumulative_regret


    def plot_results(self, regrets):
        # Creating a vector for representing the rounds
        round_vector = np.arange(1, self.T+1)

        # Plotting the cumulative regret
        plt.figure(figsize=(10,5))

        plt.title("Cumulative regrets")
        # plt.plot(round_vector, regrets[1], label='Bandits')
        plt.plot(round_vector, regrets, label='Experts')
        plt.legend()
        plt.show()


















###################################################### Data Handler ####################################################
class DataHandler:
    def __init__(self, dataPath, fname):
        self.path = dataPath
        self.fname = fname

    def extract_data(self):
        # Extracting the data from the csv file
        data_frame = pd.read_csv(self.path + '/' + self.fname, sep=',', header=None, dtype=np.float64)

        # Converting the data frame to a numpy array
        np_data = data_frame.to_numpy()
        
        # Extracting the number of experts
        numOfExperts = np_data.shape[0]

        # Extracting the number of rounds
        numOfRounds = np_data.shape[1]

        # Creating an array for every time step
        timeSteps = np.array([np_data[:, x] for x in range(numOfRounds)])

        # Returning the number of experts, the number of rounds and the data
        return numOfExperts, numOfRounds, timeSteps



######################################### Multiplicate Weights Experts - Bandits ########################################
class BExperts:
    def __init__(self, bexperts, rounds, data, type):
        # Assigning the number of 
        self.K = bexperts
        self.T = rounds
        self.data = data
        self.type = type
        self.eta = np.sqrt(np.log(self.K) / (self.T))

        print(self.type)
        # Creating an array of dictionaries containing the needed weights for the experts
        self.bexperts = np.array([{'weight': 1} for x in range(self.K)])
        
        # Creating an array of dictionaries for measuring the performance of the algorithm 
        self.algo_stats = {'optimal': [0 for _ in range(self.T)],  'selected': [0 for _ in range(self.T)], 'round_regret': [0 for _ in range(self.T)], 'cumulative_regret': [0 for _ in range(self.T)]}


    # Function for updating the bexperts variables 
    def update_bexperts(self, curr_round, bexpert_index, optimal_server_value):
        # If the type of the algorithm is experts
        if self.type == 'experts':
            # Calculating the loss for each expert
            for i in range(self.K):
                self.bexperts[i]['loss'] = self.data[curr_round][i] - optimal_server_value

            # Updating the weight of each expert
            for i in range(self.K):
                self.bexperts[i]['weight'] = np.power(1 - self.eta, self.bexperts[i]['loss']) * self.bexperts[i]['weight']

        # If the type of the algorithm is bandits
        elif self.type == 'bandits':
            # Calculating the loss foe the chosen bandit 
            self.bexperts[bexpert_index]['loss'] = self.data[curr_round][bexpert_index] - optimal_server_value

            # Updating the weight of each expert
            self.bexperts[bexpert_index]['weight'] = np.power(1 - self.eta, self.bexperts[bexpert_index]['loss']) * self.bexperts[bexpert_index]['weight']


    # Function for updating the algo_stats variables
    def updating_algoStats(self, curr_round, optimal_server_value, chosen_expert_value):
        self.algo_stats['optimal'][curr_round] = self.algo_stats['optimal'][curr_round-1] + optimal_server_value
        self.algo_stats['selected'][curr_round] = self.algo_stats['selected'][curr_round-1] + chosen_expert_value
        self.algo_stats['round_regret'][curr_round] = (self.algo_stats['selected'][curr_round] - self.algo_stats['optimal'][curr_round]) 
        self.algo_stats['cumulative_regret'][curr_round] = self.algo_stats['cumulative_regret'][curr_round - 1] + self.algo_stats['round_regret'][curr_round]
    
    
    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for curr_round in range(self.T):
            # Selecting the server with the lowest traffic
            optimal_server_index = np.argmin([self.data[curr_round]])
            optimal_server_value = self.data[curr_round][optimal_server_index]

            # Calculating the total weight of the experts
            total_weight = np.sum([bexpert['weight'] for bexpert in self.bexperts])
            
            # Selecting the expert with the highest weight value
            chosen_expert_index = np.argmax([bexpert['weight'] for bexpert in self.bexperts])                       
            chosen_expert_value = self.data[curr_round][chosen_expert_index]

            # Updating the bexperts variables
            self.update_bexperts(curr_round, chosen_expert_index, optimal_server_value)                     
            
            # Updating the algorithm statistics
            self.updating_algoStats(curr_round, optimal_server_value, chosen_expert_value)


        return self.algo_stats['round_regret'], self.algo_stats['cumulative_regret']






















if __name__ == "__main__":

    tester = Tester(dataPath='ReinforcementLearning/Exercise_2', fname='Milano_timeseries.csv', numOfRounds=7000)

    exp_regret, cum_exp_regret = tester.run_experiment_BExperts(type='experts')
    # band_regret, cum_band_regret= tester.run_experiment_BExperts(type='bandits')
    
    # print("here")
    tester.plot_results(exp_regret)
    