from dataclasses import dataclass, field


@dataclass
class CouplingExposureTracker:
    """Accumulate dimensionless nonadiabatic coupling action over a coupling region.

    The first-order electronic amplitude change scales as

        |Delta c| ~ |v.d| dt.

    Accumulating |v.d| dt therefore gives a timestep-aware spawn trigger without
    choosing an unrelated fixed instantaneous threshold for every dt.
    """
    action_threshold: float = 2.0e-4
    coupling_floor: float = 1.0e-8
    exposures: dict = field(default_factory=dict)

    def reset(self, key):
        self.exposures[key] = 0.0

    def update(self, key, coupling_rate, dt):
        rate = abs(float(coupling_rate))
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        if rate < self.coupling_floor:
            self.exposures[key] = 0.0
            return False, 0.0

        value = self.exposures.get(key, 0.0) + rate * dt
        self.exposures[key] = value
        return value >= self.action_threshold, value

    def consume(self, key):
        value = self.exposures.get(key, 0.0)
        self.exposures[key] = 0.0
        return value


def first_order_transfer_bound(coupling_rate, dt):
    return abs(float(coupling_rate)) * float(dt)
