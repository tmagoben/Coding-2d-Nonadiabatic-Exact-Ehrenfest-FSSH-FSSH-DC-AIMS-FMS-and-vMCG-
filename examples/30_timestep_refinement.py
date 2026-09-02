from gaussian_dynamics.benchmark_suite import managed_timestep_refinement

out=managed_timestep_refinement(
    dts=(8e-4,4e-4,2e-4),
    final_time=0.006,
    spa_order=0,
)

print('dt values:')
print(out['dts'])
print('\nfinal populations:')
print(out['populations'])
print('\nsuccessive population-vector errors:')
print(out['successive_errors'])
print('\nobserved refinement orders:')
print(out['observed_orders'])
