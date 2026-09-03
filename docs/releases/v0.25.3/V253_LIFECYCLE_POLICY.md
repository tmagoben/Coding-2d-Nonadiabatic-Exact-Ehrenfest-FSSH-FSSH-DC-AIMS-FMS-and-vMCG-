# v0.25.3 lifecycle policy and frozen gates

| Quantity | Default |
|---|---:|
| Adaptation interval | 1 step |
| Minimum / maximum packets | 1 / 8 |
| Minimum removal age | 2 steps |
| Maximum dormant activation age | 64 steps |
| Position / momentum displacement | 2 local standard deviations |
| Spawn residual-capture threshold | `1e-5` Hartree |
| Minimum spawn novelty | `1e-4` |
| Projection relative / absolute SVD cutoff | `1e-11` / `1e-13` |
| Maximum basis condition number | `1e8` |
| Projection linear residual tolerance | `2e-10` |
| Maximum spawn projection loss | `2e-10` |
| Shape-activation population | `1e-6` |
| Maximum ordinary prune population | `1e-8` |
| Maximum prune projection loss | `1e-8` |
| Minimum merge overlap | `0.999` |
| Maximum merge projection loss | `1e-8` |
| Maximum event energy jump | `1e-6` Hartree |

The wider merge limits used by the deterministic merge oracle are explicit settings
in that receipt; production defaults remain the conservative values above.

No diagonal metric loading, coefficient seeding, averaged merge geometry, silent
renormalization without a receipt, or multiple topology events at one checkpoint is
permitted.
