"""Reads the events from .lhe files and constructs the parton-level m_{tautau} distribution"""

import numpy as np
from PyLHE_EventAnalysis.src.Analysis import EventAnalysis, EventLoop
from PyLHE_EventAnalysis.src.Histogram import ObservableHistogram
from PyLHE_EventAnalysis.examples.FCC_hh import Observables
from PyLHE_EventAnalysis.src.Utilities import read_xsection
import pylhe
import copy
import json
from PyLHE_EventAnalysis.examples.FCC_hh.tau_leptonic.analyse_events import CutRatioMETMLL, rapidity_cut

if __name__ == "__main__":
    # Folder where the files are stored
    folder_path = "/home/martines/work/MG5_aMC_v3_1_1/PhD/FCC-hh/ditau-prod-leptonic-decays/lhe_files"

    # Folders for different masses
    leptoquark_masses = [2, 4, 6, 8, 10]

    # Histograms for the analysis
    bin_edges = list(range(0, 16400, 400)) + [1000000000000]

    # 1. Invariant mass of the charged leptons
    inv_mass_observable = Observables.invariant_mass_emu
    inv_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=inv_mass_observable)

    # Creates the EventAnalysis object - handles the selection of the event
    event_analyses = {
        # Only cuts on the rapidity and on the ratio met/mll
        "MLL": EventAnalysis(selection_cuts=[rapidity_cut, CutRatioMETMLL(0.2)])
    }

    # Creates the EventLoop object - iterates over all the events in a .lhe file and books the histogram
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram_template=inv_mass_hist)

    # Dictionary to store the constructed histograms for each mass
    mll_hists = {leptoquark_mass: copy.copy(inv_mass_hist) for leptoquark_mass in leptoquark_masses}

    # Iterates over all the masses
    for leptoquark_mass in leptoquark_masses:
        # Iterates over all the different phase-space simulated regions
        for bin_index in range(1, 42):
            # Path to the .lhe file
            lhe_filename = f"{folder_path}/mU1_{leptoquark_mass}TeV/x1L-x1L-x1L-x1L-bin-{bin_index}.lhe"
            # Reads the cross-section for the current simulated bin
            cross_section = read_xsection(path_to_file=lhe_filename)
            # Reads the total number of events in the file
            num_events = pylhe.read_num_events(filepath=lhe_filename)

            # Runs the analysis in each event - returns the disttributions constructed out of the
            # selected events in the current file
            mll_hist = event_loop.analyse_events(filename=lhe_filename, event_analyses=event_analyses)

            # Extracts the info from the histograms
            mll_hists[leptoquark_mass] += (2 * cross_section / num_events) * mll_hist["MLL"]

            print(cross_section)
            print(mll_hist["MLL"])

    # Invariant mass dists
    with open("/home/martines/work/MG5_aMC_v3_1_1/PhD/FCC-hh/betaL_mll.json", "w") as file_:
        simulations_hists = {term: dist.tolist() for term, dist in mll_hists.items()}
        json.dump(simulations_hists, file_, indent=4)
        print(".json file saved")
