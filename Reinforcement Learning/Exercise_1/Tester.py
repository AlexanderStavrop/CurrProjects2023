import epsilonGreedy
import UBC

if __name__ == '__main__':
    N = 20
    K = 10
    T = 1000
    
    epsilon_greedy = epsilonGreedy.epsilonGreedy(K, T)
    epsilon_greedy.run_experiment()

    ucb = UBC.UBC(K, T)
    ucb.run_experiment()