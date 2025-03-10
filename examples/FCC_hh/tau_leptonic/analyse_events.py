"""Reads the events from .lhe files and constructs the parton-level m_{tautau} distribution"""

import numpy as np
from PyLHE_EventAnalysis.src.Analysis import EventAnalysis, EventLoop
from PyLHE_EventAnalysis.src.Histogram import ObservableHistogram, HistogramCompound, CorrelatedHist
from PyLHE_EventAnalysis.examples.FCC_hh import Observables
from PyLHE_EventAnalysis.src.Utilities import read_xsection
import pylhe
import copy
import json


class CutRatioMETMLL:
    """Applies the cut on the met/mll ratio."""

    def __init__(self, ratio_cut):
        self.cut = ratio_cut

    def __call__(self, event: pylhe.LHEEvent) -> bool:
        """Applies the cut on the met/mll ratio."""
        # Computes the invariant mass
        mll = Observables.invariant_mass_emu(event)
        # Missing energy
        met = Observables.missing_energy(event)
        # Applies the cut
        return met > self.cut * mll


def rapidity_cut(event: pylhe.LHEEvent):
    """Rapidity of the charged lepton pair"""
    charged_lep_rap = Observables.rapidity(event)
    return abs(charged_lep_rap) < 2.5


# def taus_cut(event: pylhe.LHEEvent):
#     """Reconstructs the taus momentum and applies the original cuts."""
#     tau = Observables.evaluate_total_momentum_pids(event, [11, -12, 16])
#     antitau = Observables.evaluate_total_momentum_pids(event, [-11, 12, -16])
#
#     if tau.pt < 20 or antitau.pt < 20:
#         return False
#
#     return abs(tau.eta) < 4 and abs(antitau.eta) < 4


if __name__ == "__main__":
    # Folder where the files are stored
    folder_path = "/Users/martines/Desktop/PhD/Data/FCC-hh/ditau-leptonic/mU1_10TeV/mU1_10TeV"

    # Histograms for the analysis
    bin_edges = list(range(0, 16400, 400)) + [1000000000000]

    # 1. Invariant mass of the charged leptons
    inv_mass_observable = Observables.invariant_mass_emu
    inv_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=inv_mass_observable)

    # 2. For the met/mll distribution
    met_mll_observable = Observables.met_mll_ratio
    met_mll_dist = CorrelatedHist(xobservable=inv_mass_observable, yobservable=met_mll_observable, bin_edges=bin_edges)

    # All histograms that need to be booked for the analysis
    compound_hist = HistogramCompound(histograms={"MLL": inv_mass_hist})

    # Creates the EventAnalysis object - handles the selection of the event
    event_analyses = {
        # "no_cut": EventAnalysis(selection_cuts=[taus_cut]),
        # "MET_MLL_1": EventAnalysis(selection_cuts=[rapidity_cut, CutRatioMETMLL(0.1)]),
        "MET_MLL_2": EventAnalysis(selection_cuts=[rapidity_cut, CutRatioMETMLL(0.2)])
    }

    # Creates the EventLoop object - iterates over all the events in a .lhe file and books the histogram
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram_template=compound_hist)

    # Dictionary to store the constructed histograms
    mll_hists = {analysis: copy.copy(inv_mass_hist) for analysis in event_analyses}

    # MET/Mll distribution
    met_mll_sum = np.zeros(len(bin_edges) - 1)
    event_count = copy.copy(inv_mass_hist)

    # Iterates over all the different phase-space simulated regions
    for bin_index in range(1, 42):
        # Path to the .lhe file
        lhe_filename = f"{folder_path}/x1L-x1L-x1L-x1L-bin-{bin_index}.lhe"
        # Reads the cross-section for the current simulated bin
        cross_section = read_xsection(path_to_file=lhe_filename)
        # Reads the total number of events in the file
        num_events = pylhe.read_num_events(filepath=lhe_filename)

        # Runs the analysis in each event - returns the disttributions constructed out of the
        # selected events in the current file
        bin_hist = event_loop.analyse_events(
            filename=lhe_filename, event_analyses=event_analyses
        )

        # Extracts the info from the histograms
        for analysis in event_analyses:
            mll_hists[analysis] += (2 * cross_section / num_events) * bin_hist[analysis].get_hist("MLL")

        # For the met/mll plot
        # met_mll_sum += bin_hist["no_cut"].get_hist("MET_MLL").bin_sum
        # event_count += bin_hist["no_cut"].get_hist("MLL")

        print(cross_section)
        print(bin_hist["no_cut"].get_hist("MLL"))
        # print(bin_hist["no_cut"].get_hist("MET_MLL").bin_sum)

    # Invariant mass dists
    with open(f"{folder_path}/LQ_ML_tau_cuts.json", "w") as file_:
        simulations_hists = {term: dist.tolist() for term, dist in mll_hists.items()}
        json.dump(simulations_hists, file_, indent=4)

    # with open(f"{folder_path}/LQ_MET_ML_tau_cuts.json", "w") as file_:
    #     json.dump((met_mll_sum / event_count).tolist(), file_, indent=4)
