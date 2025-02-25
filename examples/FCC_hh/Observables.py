"""Collects the observables for the analysis."""

import numpy as np
import pylhe


class InvariantMass:
    """
    Evaluates the invariant mass of the final particles in the event.
    If a list of PID is given, only the particles with the PIDs in the list are taken into account.
    """

    def __init__(self, part_pids: list = None):
        self._pids = part_pids      # Stores the final particle PIDs

    def __call__(self, event: pylhe.LHEEvent):
        """Evaluates the invariant mass of the event."""
        # Four-momentum the particles in the event
        momentum = np.array([
            [getattr(part, comp) for comp in "e px py pz".split()]
            for part in event.particles if part.status == 1 and self._include_particle(part.id)
        ])
        # Total four-momentum
        event_momentum = np.sum(momentum, axis=0)
        # Invariant mass
        invariant_mass = np.sum(
            [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(event_momentum)]
        )
        return np.sqrt(invariant_mass)

    def _include_particle(self, part_pid) -> bool:
        """Checks if the particle must be included in the computation."""
        if self._pids is not None:
            return abs(part_pid) in self._pids
        return True


def missing_energy(event: pylhe.LHEEvent):
    """Computes the missing energy of the event"""
    # Transverse momentum of all neutrinos
    transv_momentum_nu = np.array([
        [getattr(part, comp) for comp in "px py".split()]
        for part in event.particles if part.status == 1 and abs(part.id) in [12, 14, 16]
    ])
    # Missing energy of the event
    met = np.sqrt(np.sum([np.power(comp, 2) for comp in np.sum(transv_momentum_nu, axis=0)]))

    return met


def met_mll_ratio(event: pylhe.LHEEvent):
    """Evaluates the MET/mll ratio of the event."""
    inv_mass = InvariantMass(part_pids=[11, 13])(event)    # Invariant mass
    met = missing_energy(event)   # missing energy
    # Returns the ratio
    return met / inv_mass
