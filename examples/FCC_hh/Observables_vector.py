"""Collects the observables for the analysis."""

import pylhe
import vector


def evaluate_total_momentum(event: pylhe.LHEEvent, part_pids: list):
    """Calculates the total momentum taking into account only the particles with PIDs in the part_pid list."""
    # Total momentum vector
    total_momentum = vector.MomentumObject4D(px=0, py=0, pz=0, e=0)
    for part in event.particles:
        if abs(part.id) in part_pids and part.status == 1:
            total_momentum += vector.MomentumObject4D(**{comp: getattr(part, comp) for comp in "e px py pz".split()})
    return total_momentum


def evaluate_total_momentum_pids(event: pylhe.LHEEvent, part_pids: list):
    """Calculates the total momentum taking into account only the particles with PIDs in the part_pid list."""
    # Total momentum vector
    total_momentum = vector.MomentumObject4D(px=0, py=0, pz=0, e=0)
    for part in event.particles:
        if part.id in part_pids and part.status == 1:
            total_momentum += vector.MomentumObject4D(**{comp: getattr(part, comp) for comp in "e px py pz".split()})
    return total_momentum


def invariant_mass_emu(event: pylhe.LHEEvent):
    """Evaluates the invariant mass of the e-mu pair"""
    charged_leptons_momentum = evaluate_total_momentum(event, [11, 13])
    return charged_leptons_momentum.m


def missing_energy(event: pylhe.LHEEvent):
    """Calculates the missing energy of the event."""
    nu_total = evaluate_total_momentum(event, [12, 14, 16])
    return nu_total.pt


def met_mll_ratio(event: pylhe.LHEEvent):
    """Evaluates the MET/mll ratio of the event."""
    inv_mass = invariant_mass_emu(event)
    met = missing_energy(event)   # missing energy
    # Returns the ratio
    return met / inv_mass


def pseudo_rapidity(event: pylhe.LHEEvent):
    """Rapidity"""
    charged_lep_m = evaluate_total_momentum(event, [11, 13])
    return charged_lep_m.eta


def rapidity(event: pylhe.LHEEvent):
    """Rapidity"""
    charged_lep_m = evaluate_total_momentum(event, [11, 13])
    return charged_lep_m.rapidity
