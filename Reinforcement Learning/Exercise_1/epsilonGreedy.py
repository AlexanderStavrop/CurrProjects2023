import numpy as np

# Creating a class for espilonGreedy algorithm
class epsilonGreedy:
    def __init__(self, bandits, rounds):
        self.K = bandits
        self.T = rounds

        # Creating a variable for epsilon
        self.epsilon = 0

        # Creating an array of individual dictionaries for each bandit
        self.bandits = np.array([{'prob':0, 'reward':0, 'score':0, 'Q':0, 'mu':0} for _ in range(self.K)])

        # Creating the bandits
        self.create_bandits()
        

    # Function for creating the bandits
    def create_bandits(self):

        # [self.bandits[x].update([('prob', np.random.uniform(0, 1))]) for x in range(self.K)]

        # For each bandit
        for i in range(self.K):
            # Assigning a random uniform probability to each bandit
            self.bandits[i].update([('prob', np.random.uniform(0, 1))])

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
        # If the random number is less than epsilon, explore
        if np.random.random() < self.epsilon:
            # Return a random arm
            return np.random.choice(self.K)
        # Else, exploit
        else:
            # Return the arm with the highest mu value
            return np.argmax([bandit['mu'] for bandit in self.bandits])


    # Function for updating the measurements
    def update_measurements(self, arm):
        self.bandits[arm]['Q'] += 1
        self.bandits[arm]['score'] += self.bandits[arm].get('reward') * np.random.binomial(1,p=self.bandits[arm].get('prob'))
        self.bandits[arm]['mu'] = self.bandits[arm].get('score') / self.bandits[arm].get('Q')


    # Function for running the algorithm and extracting the results
    def run_experiment(self):
        # For each round
        for curr_round in range(self.T):

            # Updating the epsilon value
            self.epsilon = (1 + curr_round)**(-1/3) * (np.log10(1 + curr_round))**(1/3)
            
            # Choosing the arm
            arm = self.choose_action()

            # Updating the measurements using the chosen arm
            self.update_measurements(arm)
            

        # Printing the results
        for i in range(self.K):
            print('prob = %.2f, reward = %.2f, score = %6.2f, Q = %3d, mu = %.2f' %  (self.bandits[i].get('prob'), self.bandits[i].get('reward'), self.bandits[i].get('score'), self.bandits[i].get('Q'), self.bandits[i].get('mu')))


if __name__ == '__main__':
    K = 10
    T = 1000

    # epsilon = epsilonGreedy(K, T)
    # epsilon.run_experiment()



