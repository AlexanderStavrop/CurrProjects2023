import numpy as np

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
