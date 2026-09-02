from gaussian_dynamics.adaptive_spawning import CouplingExposureTracker


def trigger_time(dt,rate=2.0,threshold=0.1):
    tracker=CouplingExposureTracker(action_threshold=threshold,coupling_floor=0.0)
    t=0.0
    for _ in range(10000):
        t+=dt
        ready,action=tracker.update(('parent',1),rate,dt)
        if ready:
            return t,action
    raise RuntimeError('no trigger')

for dt in [0.01,0.005,0.0025]:
    t,action=trigger_time(dt)
    print(f'dt={dt:8.4g}  trigger_time={t:8.4g}  accumulated_action={action:8.4g}')

print('\nThe integrated |v.d| dt criterion is substantially less tied to an arbitrary per-step threshold.')
