import numpy as np

from gaussian_dynamics.benchmark_suite import spa_order_comparison

out=spa_order_comparison(dt=2e-4,final_time=0.006)

print('SPA0 final populations:',out['spa0_populations'])
print('SPA1 final populations:',out['spa1_populations'])
print('L2 difference:',out['difference_l2'])
print('\nThe difference is a controlled approximation diagnostic for this short analytic CI run.')
