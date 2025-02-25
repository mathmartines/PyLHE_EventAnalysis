"""Reads the events from .lhe files and constructs the parton-level m_{tautau} distribution"""

from PyLHE_EventAnalysis.src.Analysis import EventAnalysis, EventLoop
from PyLHE_EventAnalysis.src.Histogram import ObservableHistogram, HistogramCompound, CorrelatedHist
from PyLHE_EventAnalysis.examples.FCC_hh import Observables
from PyLHE_EventAnalysis.src.Utilities import read_xsection
import pylhe
import copy
import json


def met_mll_ratio_cut(event: pylhe.LHEEvent) -> bool:
    """Applies the cut on the met/mll ratio."""
    # Computes the invariant mass
    mll = Observables.InvariantMass(part_pids=[11, 13])(event)
    # Missing energy
    met = Observables.missing_energy(event)
    # Applies the cut
    return met > 0.1 * mll


if __name__ == "__main__":
    # Folder where the files are stored
    folder_path = "/home/martines/work/MG5_aMC_v3_1_1/PhD/FCC-hh/ditau-prod-leptonic-decays/lhe_files/mU1_8TeV"

    # Cross-section terms we simulated
    xsec_terms = ["x1L-x1L-x1L-x1L"]

    # Histograms for the analysis
    bin_edges = list(range(0, 16400, 400))
    # 1. Invariant mass of the charged leptons
    inv_mass_observable = Observables.InvariantMass(part_pids=[11, 13])
    inv_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=inv_mass_observable)

    # Creates the EventAnalysis object - handles the selection of the event
    event_analysis = EventAnalysis(selection_cuts=[met_mll_ratio_cut])  # No cuts being applied
    # Creates the EventLoop object - iterates over all the events in a .lhe file and books the histogram
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram_template=inv_mass_hist)

    # Dictionary to store the constructed histograms
    simulations_hists = {}

    # Iterates over all the terms
    for term in xsec_terms:
        # Histograms for the current term
        simulations_hists[term] = copy.copy(inv_mass_hist)

        # Iterates over all the different phase-space simulated regions
        for bin_index in range(1, 41):
            # Path to the .lhe file
            lhe_filename = f"{folder_path}/{term}-bin-{bin_index}.lhe"
            # Reads the cross-section for the current simulated bin
            cross_section = read_xsection(path_to_file=lhe_filename)
            # Reads the total number of events in the file
            num_events = pylhe.read_num_events(filepath=lhe_filename)

            # Runs the analysis in each event - returns the disttributions constructed out of the
            # selected events in the current file
            bin_hist = event_loop.analyse_events(
                filename=lhe_filename, event_analysis=event_analysis
            )

            print(cross_section)
            print(bin_hist)

            # Updates the full histogram
            simulations_hists[term] += (2 * cross_section / num_events) * bin_hist

    # Saves the distribution in a .json file
    with open(f"{folder_path}/mU1_8TeV_MET_MLL_1.json", "w") as file_:
        simulations_hists = {term: dist.tolist() for term, dist in simulations_hists.items()}
        json.dump(simulations_hists, file_, indent=4)
