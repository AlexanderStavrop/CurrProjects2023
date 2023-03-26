import numpy as np
import matplotlib.pyplot as plt


######################################################## Tester ########################################################
class Tester:
    def __init__(self, bandits, rounds):
        # Assigning the number of bandits and the rounds to the class
        self.K = bandits
        self.T = rounds

        # Creating an array of individual dictionaries for each bandit
        self.bandits = np.array([{'prob': 0, 'reward': 0, 'score': 0, 'Q': 0, 'mu': 0} for _ in range(self.K)])
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

    def find_best_arm(self):
        # Calculating the product of the probability of success and the reward of each arm
        prob_reward = [bandit['prob'] * bandit['reward'] for bandit in self.bandits]

        # Choosing the arm with the highest product
        self.best_bandit['arm'] = np.argmax(prob_reward)

        # Updating the score of the arm based on the product
        self.best_bandit['arm_score'] = prob_reward[self.best_bandit['arm']]

    def run_EpsilonGreedy(self):
        # Creating the epsilonGreedy algorithm
        epsilon_greedy = EpsilonGreedy(self.K, self.T, self.bandits, self.best_bandit)

        # Running the algorithm and returning the regret
        return epsilon_greedy.run_experiment()


#################################################### Epsilon Greedy ####################################################

# Class for implementing the epsilonGreedy algorithm
class EpsilonGreedy:
    def __init__(self, numOfBandits, rounds, bandits, best_bandit):
        # Assigning the number of bandits and the rounds to the class
        self.K = numOfBandits
        self.T = rounds
        self.bandits = bandits
        self.best_bandit = best_bandit

        # Creating a dictionary containing the needed values for measuring the performance of the algorithm
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
            epsilon = (1 + curr_round)**(-1/3) * (np.log10(1 + curr_round))**(1/3)

            # Choosing the arm
            arm = self.choose_action(epsilon)

            # Updating the measurements using the chosen arm
            self.update_measurements(arm, curr_round)

            # Calculating the regret
            self.calculate_regret()

        return self.algo_stats['regret']


###################################################################### UCB ######################################################################
# Creating a class for espilonGreedy algorithm
class UBC:
    def __init__(self, k, t):
        self.K = k
        self.T = t

        # Creating the needed variables
        self.bandits = np.array([{'prob':0, 'reward':0, 'score':0, 'Q':0, 'mu':0, 'ubc':2**23} for _ in range(self.K)])

        # Creating the bandits
        self.create_bandits()


    # Function for creating the bandits
    def create_bandits(self):
        [self.bandits[x].update([('prob', np.random.uniform(0, 1))]) for x in range(self.K)]
        
        # For each bandit
        for i in range(self.K):
            # Creating a and b values
            a = np.random.random()
            b = np.random.random()

            # Checking which one is bigger and assigning the uniform success probability accordingly
            if a < b:
                self.bandits[i]['reward'] = np.random.uniform(a, b)
            else:
                self.bandits[i]['reward'] = np.random.uniform(b, a)
        

    # Function for choosing between exploration and exploitation
    def choose_action(self):
        return np.argmax([bandit['ubc'] for bandit in self.bandits])
    

    # Function for updating the measurements
    def update_measurements(self, arm):
        self.bandits[arm]['Q'] += 1
        self.bandits[arm]['score'] += self.bandits[arm].get('reward')


    # Function for updating the mu value
    def update_mu(self, arm):
        self.bandits[arm]['mu'] = self.bandits[arm].get('score') / self.bandits[arm].get('Q')


    # Function for updating the ucb value
    def update_ubc(self, arm):
        self.bandits[arm]['ubc'] = self.bandits[arm]['mu'] + np.sqrt(np.log(self.T) / self.bandits[arm]['Q'])


    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for i in range(self.T):
            # Choosing the arm
            arm = self.choose_action()

            # Updating the measurements using the chosen arm
            self.update_measurements(arm)

            # Update the mu value
            self.update_mu(arm)

            # Updating the ucb value
            self.update_ubc(arm)
            

        # Printing the results
        for i in range(self.K):
            print('prob = %.2f, reward = %.2f, score = %6.2f, mu = %.2f, Q = %3d, ubc = %.3f' %  (self.bandits[i].get('prob'), self.bandits[i].get('reward'), self.bandits[i].get('score'), self.bandits[i].get('mu'), self.bandits[i].get('Q'), self.bandits[i].get('ubc')))




def plot_regret(rounds, regret, label):
    plt.plot(np.arange(1, rounds+1), regret, label=label)
    plt.xlabel('Rounds')
    plt.ylabel('Regret')
    plt.legend()
    plt.savefig('regret_' + str(rounds) + '.png')
    print("done")

def print_stats(bandits_list):
        for i in range(10):
            prob = bandits_list[i]['prob'] * bandits_list[i]['reward']
            print('%d) prob = %.2f, reward = %.2f, score = %6.2f, mu = %.2f, Q = %3d , aaa = %.2f' %  ((i), bandits_list[i].get('prob'), bandits_list[i].get('reward'), bandits_list[i].get('score'), bandits_list[i].get('mu'), bandits_list[i].get('Q'), prob))

###################################################################### Main #####################################################################
if __name__ == '__main__':
    K = 10
    T = 1000

    # Create bandits
    tester = Tester(K, T)

    # Run the experiment for each algorithm
    reg = tester.run_EpsilonGreedy()
    plot_regret(T, reg, label='Epsilon Greedy')




    # epsilon = EpsilonGreedy(K, T, bandits_list, chosen_bandit)
    # epsilon_regret = epsilon.run_experiment()

    # plot_regret(T, epsilon_regret, 'Epsilon Greedy')

