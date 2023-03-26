import numpy as np
import matplotlib.pyplot as plt


######################################################## Tester ########################################################
class Tester:
    def __init__(self, bandits, rounds):
        # Assigning the number of bandits and the rounds to the class
        self.K = bandits
        self.T = rounds

        # Creating an array of individual dictionaries for each bandit
        self.bandits = np.array([{'prob': 0, 'reward': 0} for _ in range(self.K)])
        self.best_bandit = {'arm': 0,  "arm_score": 0, 'round_score': [0 for _ in range(self.T)]}

        # Creating the bandits
        self.create_bandits()

        # Choosing the best arm based on the data currently available
        self.find_best_arm()


    # Function for creating the bandits
    def create_bandits(self):
        # For each bandit
        for i in range(self.K):
            # Assigning a random uniform probability to each bandit
            self.bandits[i]['prob'] = np.random.uniform(0, 1)

            # Creating a and b values
            a = np.random.random()
            b = np.random.random()

            # Checking which one is bigger and assigning the uniform success probability accordingly
            if a < b:
                self.bandits[i]['reward'] = np.random.uniform(a, b)
            else:
                self.bandits[i]['reward'] = np.random.uniform(b, a)

    
    # Function for finding the best arm based on the data currently available
    def find_best_arm(self):
        # Calculating the product of the probability of success and the reward of each arm
        prob_reward = [bandit['prob'] * bandit['reward'] for bandit in self.bandits]

        # Choosing the arm with the highest product
        self.best_bandit['arm'] = np.argmax(prob_reward)

        # Updating the score of the arm based on the product
        self.best_bandit['arm_score'] = prob_reward[self.best_bandit['arm']]


    # Function for running the epsilonGreedy algorithm and returning the regret
    def run_EpsilonGreedy(self):
        # Creating the epsilonGreedy algorithm
        epsilon_greedy = EpsilonGreedy(self.K, self.T, self.bandits, self.best_bandit)

        # Running the algorithm and returning the regret
        return epsilon_greedy.run_experiment()

    
    # Function for running the epsilonGreedy algorithm and returning the regret
    def run_UCB(self):
        # Creating the epsilonGreedy algorithm
        ucb = UCB(self.K, self.T, self.bandits, self.best_bandit)

        # Running the algorithm and returning the regret
        return ucb.run_experiment()
    

    # Function for creating a vector containing the cumulative 
    def calculate_cum_regret(self, epsilon_regret, usb_regret):
        cum_epsilon = [0 for _ in range(T)]
        cum_ucb = [0 for _ in range(T)]

        for i in range(self.T):
            cum_epsilon[i] = cum_epsilon[i-1] + epsilon_regret[i]
            cum_ucb[i] = cum_ucb[i-1] + ucb_regret[i]
            
        return cum_epsilon, cum_ucb
    

    # Function for plotting the regret and the cumulative regret of the algorithms
    def plot_regret(self, epsilon_regret, label_epsilon, ucb_regret, label_ucb):
        # Creating a vector for representing the rounds
        round_vector = np.arange(1, self.T+1)
        cum_epsilon, cum_ucb = self.calculate_cum_regret(epsilon_regret, ucb_regret)
        
        EpsilonGreedy_complexity = (round_vector + 1)**(-1/3) * (self.K * np.log10(1 + round_vector))**(1/3)
        ucb_complexity = np.sqrt(self.K * self.T * np.log10(T))

        # Plotting the regret for the epsilonGreedy algorithm
        plt.figure(figsize=(10,5))
        plt.subplot(1, 2, 1)
        plt.plot(round_vector, epsilon_regret, label=label_epsilon)
        plt.plot(round_vector, ucb_regret, label=label_ucb)

        plt.title("Regret for epsilonGreedy and UCB")
        plt.xlabel('Rounds')
        plt.ylabel('Regret')
        plt.legend()
    

        # Plotting the cumulative regret
        plt.subplot(1, 2, 2)
        plt.plot(round_vector, EpsilonGreedy_complexity,label='O(t^2/3 * (k*log(t))^1/3)')
        plt.plot(round_vector, cum_epsilon,label='Cumulative regret of ε-Greedy')
        # plt.plot(round_vector, )
        plt.plot(round_vector, cum_ucb,label='Cumulative regret of ucb')
        plt.legend()

        plt.savefig('Reinforcement Learning/Exercise_1/Review/Images/Regret' + str(self.K) + '_' + str(self.T) + '.eps', format='eps')
        print("done")

#################################################### Epsilon Greedy ####################################################

# Class for implementing the epsilonGreedy algorithm
class EpsilonGreedy:
    def __init__(self, numOfBandits, rounds, bandits, best_bandit):
        # Assigning the number of bandits and the rounds to the class
        self.K = numOfBandits
        self.T = rounds
        self.best_bandit = best_bandit


        # Creating two dictionaries containing the needed values for measuring the performance of the algorithm
        self.bandits = np.array([{'prob': bandits[x]['prob'], 'reward': bandits[x]['reward'], 'score': 0, 'Q': 0, 'mu': 0} for x in range(self.K)])
        self.algo_stats = {'round_score': [0 for _ in range(self.T)], 'algo_score': [0 for _ in range(self.T)], 'regret': [0 for _ in range(self.T)]}


    # Function for choosing between exploration and exploitation
    def choose_action(self, epsilon):
        # If the random number is less than epsilon, explore
        if np.random.random() < epsilon:
            # Return a random arm
            return np.random.choice(self.K)
        # Else, exploit
        else:
            # Return the arm with the highest mu value
            return np.argmax([bandit['mu'] for bandit in self.bandits])


    # Function for updating the measurements
    def update_measurements(self, arm, round):
        # Incrementing the number the arm was chosen
        self.bandits[arm]['Q'] += 1

        # Calculating the reward based on the probability of success and the reward of the arm
        score = self.bandits[arm]['reward'] * np.random.binomial(1,p=self.bandits[arm]['prob'])

        # Updating the score of the arm and saving the score of the round accordingly
        self.bandits[arm]['score'] += score
        self.algo_stats['round_score'][round] = score

        # Updating the mu value of the arm accordingly
        self.bandits[arm]['mu'] = self.bandits[arm]['score'] / self.bandits[arm]['Q']


    # Function for calculating the regret of the algorithm
    def calculate_regret(self):
        for i in range(self.T):
            # Calculating the best score
            if i > 0:
                # Assigning the score of the best arm to the algorithm score
                self.best_bandit['round_score'][i] = self.best_bandit['round_score'][i-1] + self.best_bandit['arm_score']

                # Assigning the score of the epsilon greedy algorithm based on the score of each round
                self.algo_stats['algo_score'][i] = self.algo_stats['algo_score'][i-1] + self.algo_stats['round_score'][i]

            else:
                # Assigning the score of the best arm to the algorithm score for the first round
                self.best_bandit['round_score'][i] = self.best_bandit['round_score'][i]

                # Assigning the score of the epsilon greedy algorithm based on the score of each round
                self.algo_stats['algo_score'][i] = self.algo_stats['round_score'][i]

            # Calculating the regret
            self.algo_stats['regret'][i] = (self.best_bandit['round_score'][i] - self.algo_stats['algo_score'][i]) / (i + 1)


    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for curr_round in range(self.T):

            # Updating the epsilon value
            epsilon = (1 + curr_round)**(-1/3) * (self.K * np.log10(1 + curr_round))**(1/3)

            # Choosing the arm
            arm = self.choose_action(epsilon)

            # Updating the measurements using the chosen arm
            self.update_measurements(arm, curr_round)

        # Calculating the regret
        self.calculate_regret()

        return self.algo_stats['regret']

###################################################################### UCB ######################################################################
# Creating a class for UCB algorithm
class UCB:
    def __init__(self, numOfBandits, rounds, bandits, best_bandit):
        self.K = numOfBandits
        self.T = rounds
        self.best_bandit = best_bandit
        for i in range(self.T):
            self.best_bandit['round_score'][i] = 0

        # self.best_bandit['round_score'][x] for x in range(self.T)
        # self.best_bandit = np.array([{'prob': best_bandit[x]['prob'], 'reward': best_bandit[x]['reward'], 'score':0} for x in range(self.K)])

        # Creating two dictionaries containing the needed values for measuring the performance of the algorithm
        self.bandits = np.array([{'prob': bandits[x]['prob'], 'reward': bandits[x]['reward'], 'score':0, 'Q':0, 'mu':0, 'ucb':2**10} for x in range(self.K)])
        self.algo_stats = {'round_score': [0 for _ in range(self.T)], 'algo_score': [0 for _ in range(self.T)], 'regret': [0 for _ in range(self.T)]}

    # Function for updating the measurements
    def update_measurements(self, arm, round):
        # Incrementing the number the arm was chosen
        self.bandits[arm]['Q'] += 1

        # Calculating the reward based on the probability of success and the reward of the arm
        score = self.bandits[arm]['reward'] * np.random.binomial(1,p=self.bandits[arm]['prob'])

        # Updating the score of the arm and saving the score of the round accordingly
        self.bandits[arm]['score'] += score
        self.algo_stats['round_score'][round] = score

        # Updating the mu value of the arm accordingly
        self.bandits[arm]['mu'] = self.bandits[arm]['score'] / self.bandits[arm]['Q']

        # Updating the ucb value of the arm accordingly
        self.bandits[arm]['ucb'] = self.bandits[arm]['mu'] + np.sqrt(np.log(self.T) / self.bandits[arm]['Q'])


    # Function for calculating the regret of the algorithm
    def calculate_regret(self):
        for i in range(self.T):
            # Calculating the best score
            if i > 0:
                # Assigning the score of the best arm to the algorithm score
                self.best_bandit['round_score'][i] = self.best_bandit['round_score'][i-1] + self.best_bandit['arm_score']

                # Assigning the score of the epsilon greedy algorithm based on the score of each round
                self.algo_stats['algo_score'][i] = self.algo_stats['algo_score'][i-1] + self.algo_stats['round_score'][i]

            else:
                # Assigning the score of the best arm to the algorithm score for the first round
                self.best_bandit['round_score'][i] = self.best_bandit['round_score'][i]

                # Assigning the score of the epsilon greedy algorithm based on the score of each round
                self.algo_stats['algo_score'][i] = self.algo_stats['round_score'][i]

            # Calculating the regret
            self.algo_stats['regret'][i] = (self.best_bandit['round_score'][i] - self.algo_stats['algo_score'][i]) / (i + 1)


    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for curr_round in range(self.T):
            # Choosing the arm
            arm = np.argmax([bandit['ucb'] for bandit in self.bandits])

            # Updating the measurements using the chosen arm
            self.update_measurements(arm, curr_round)
            
        # Calculating the regret
        self.calculate_regret()

        return self.algo_stats['regret']

###################################################################### Main #####################################################################
if __name__ == '__main__':
    K = 10
    T = 10000

    ########################################################### Running the test once ###########################################################
    # Create bandits
    tester = Tester(K, T)

    # Running the epsilon greedy algorithm
    epsilon_regret = tester.run_EpsilonGreedy()

    # Running the UCB algorithm
    ucb_regret = tester.run_UCB()
    
    # Plotting the results
    tester.plot_regret(epsilon_regret, 'Epsilon Greedy', ucb_regret, 'UCB')


    # epsilon_regrets = [0 for _ in range (tests)] 
    # ucb_regrets = [0 for _ in range (tests)]

    # for i in range(tests):
    #     # Create bandits
    #     tester = Tester(K, T)

    #     # Running the epsilon greedy algorithm
    #     epsilon_regrets[i] = tester.run_EpsilonGreedy()
        
    #     # Running the UCB algorithm
    #     ucb_regrets[i] = tester.run_UCB()


    # final_epsilon_regret = [0 for _ in range (T)] 
    # final_ucb_regret = [0 for _ in range (T)] 

    # for round in range(T):
    #     print(round)
    #     e_reg = 0
    #     usb_reg = 0

    #     for test in range(tests):
    #         e_reg =+ epsilon_regrets[test][round]
    #         ucb_reg =+ ucb_regrets[test][round]

    #     final_epsilon_regret[round] = e_reg / tests
    #     final_ucb_regret[round] = ucb_reg / tests
         
    # # Plotting the results
    # tester.plot_regret(final_epsilon_regret, 'Epsilon Greedy', final_ucb_regret, 'UCB')