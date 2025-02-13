"""Defines the important functions for the analysis"""

import pylhe
import numpy as np


def invariant_mass_taus(event: pylhe.LHEEvent):
    """Computes the invariant mass for a pair of taus"""
    # Four-momentum of each tau in the event
    momentum = np.array([
        [getattr(tau, comp) for comp in "e px py pz".split()]
        for tau in event.particles if abs(tau.id) == 15
    ])
    # Total four-momentum of the event
    event_momentum = np.sum(momentum, axis=0)
    # Invariant mass of the event
    invariant_mass = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(event_momentum)]
    )
    return np.sqrt(invariant_mass)


def taus_pt_cut(event: pylhe.LHEEvent) -> bool:
    """Selects only events where we have two taus with pT > 250 GeV."""
    # Transverse momentum of each tau in the event
    taus_pt = np.array([np.sqrt(tau.px**2 + tau.py**2) for tau in event.particles if abs(tau.id) == 15])
    # Returns True if both taus have pT > 250 GeV
    return sum(taus_pt > 250) > 1
