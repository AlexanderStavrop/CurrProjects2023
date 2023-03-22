import numpy as np

# Creating a class for espilonGreedy algorithm
class epsilonGreedy:
    def __init__(self, k, t):
        self.K = k
        self.T = t

        # Creating the needed variables
        self.epsilon = 0
        self.bandits = np.array([{'prob':0, 'reward':0, 'score':0, 'Q':0, 'mu':0} for _ in range(self.K)])

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
        

    # Function for updating the epsilon value
    def update_epsilon(self, curr_round):
        self.epsilon = 5**(-1/3) * (np.log(1 + curr_round))**(1/3)


    # Function for choosing between exploration and exploitation
    def choose_action(self):
        # If the random number is less than epsilon, explore
        if np.random.random() < self.epsilon:
            return np.random.choice(self.K)
        # Else, exploit
        else:
            return np.argmax([bandit['mu'] for bandit in self.bandits])


    # Function for updating the measurements
    def update_measurements(self, arm):
        self.bandits[arm]['Q'] += 1
        self.bandits[arm]['score'] += self.bandits[arm].get('reward')


    # Function for updating the mu value
    def update_mu(self, arm):
        self.bandits[arm]['mu'] = self.bandits[arm].get('score') / self.bandits[arm].get('Q')


    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for i in range(self.T):
            # Updating the epsilon value
            self.update_epsilon(i)

            # Choosing the arm
            arm = self.choose_action()

            # Updating the measurements using the chosen arm
            self.update_measurements(arm)

            # Update the mu value
            self.update_mu(arm)


        # Printing the results
        for i in range(self.K):
            print('prob = %.2f, reward = %.2f, score = %6.2f, Q = %3d, mu = %.2f' %  (self.bandits[i].get('prob'), self.bandits[i].get('reward'), self.bandits[i].get('score'), self.bandits[i].get('Q'), self.bandits[i].get('mu')))


if __name__ == '__main__':
    N = 20
    K = 10
    T = 1000
    
    epsilon = epsilonGreedy(K, T)
    epsilon.run_experiment()



