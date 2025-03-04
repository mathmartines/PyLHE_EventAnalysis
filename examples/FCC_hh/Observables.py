"""Collects the observables for the analysis."""

import pylhe
import numpy as np


def evaluate_total_momentum(event: pylhe.LHEEvent, part_pids: list):
    """Calculates the total momentum taking into account only the particles with PIDs in the part_pid list."""
    # Total momentum vector
    # Four-momentum of each tau in the event
    momentum = np.array([
        [getattr(part, comp) for comp in "e px py pz".split()]
        for part in event.particles if abs(part.id) in part_pids
    ])
    # Total four-momentum of the event
    return np.sum(momentum, axis=0)


def invariant_mass_emu(event: pylhe.LHEEvent):
    """Evaluates the invariant mass of the e-mu pair"""
    charged_leptons_momentum = evaluate_total_momentum(event, [11, 13])
    # Invariant mass of the event
    invariant_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(charged_leptons_momentum)]
    )
    return np.sqrt(invariant_mass_squared)


def missing_energy(event: pylhe.LHEEvent):
    """Calculates the missing energy of the event."""
    nu_total = evaluate_total_momentum(event, [12, 14, 16])
    # Transverse momentum
    return np.sqrt(np.power(nu_total[1], 2) + np.power(nu_total[2], 2))


def met_mll_ratio(event: pylhe.LHEEvent):
    """Evaluates the MET/mll ratio of the event."""
    inv_mass = invariant_mass_emu(event)
    met = missing_energy(event)  # missing energy
    # Returns the ratio
    return met / inv_mass


def rapidity(event: pylhe.LHEEvent):
    """Rapidity"""
    charged_lep_m = evaluate_total_momentum(event, [11, 13])
    rap = 0.5 * np.log((charged_lep_m[0] + charged_lep_m[3]) / (charged_lep_m[0] - charged_lep_m[3]))
    return rap
