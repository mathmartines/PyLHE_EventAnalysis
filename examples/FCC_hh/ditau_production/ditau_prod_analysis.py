"""Reads the events from .lhe files and constructs the parton-level m_{tautau} distribution"""

from PyLHE_EventAnalysis.src.Analysis import EventAnalysis, EventLoop
from PyLHE_EventAnalysis.src.Histogram import ObservableHistogram
from PyLHE_EventAnalysis.examples.FCC_hh.ditau_production import analysis_funcs
from PyLHE_EventAnalysis.src.Utilities import read_xsection
import pylhe
import copy
import json

if __name__ == "__main__":
    # Folder where the files are stored
    folder_path = "/home/martines/work/MG5_aMC_v3_1_1/PhD/FCC-hh/ditau-production/lhe_files/mU1_4TeV"

    # Cross-section terms we simulated
    xsec_terms = [
        "x1L-x1L", "x1L-x1L-x1L-x1L"
    ]

    # Definition of the histogram - by default the last bin contains overflow
    # bin_edges = list(range(0, 15400, 400))
    bin_edges = [0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 40000]
    mtautau_distribution = ObservableHistogram(bin_edges=bin_edges, observable=analysis_funcs.invariant_mass_taus)

    # Creates the EventAnalysis object - handles the selection of the event
    event_analysis = EventAnalysis(selection_cuts=[])  # No cuts being applied
    # Creates the EventLoop object - iterates over all the events in a .lhe file and books the histogram
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram_template=mtautau_distribution)

    # Dictionary to store the constructed m_tautau distribution of each term
    terms_dists = {}

    # Iterates over all the terms
    for term in xsec_terms:
        # Initializes an empty histogram for the current term
        terms_dists[term] = copy.copy(mtautau_distribution)
        # Iterates over all the different phase-space simulated regions
        for bin_index in range(1, 13):
            # Path to the .lhe file
            lhe_filename = f"{folder_path}/{term}-bin-{bin_index}.lhe"
            # Reads the cross-section for the current simulated bin
            cross_section = read_xsection(path_to_file=lhe_filename)
            # Reads the total number of events in the file
            num_events = pylhe.read_num_events(filepath=lhe_filename)
            # Runs the analysis in each event - returns the mtautau dist constructed out of the
            # selected events in the current file
            bin_mtautau_dist = event_loop.analyse_events(
                filename=lhe_filename, event_analysis=event_analysis
            )
            # Prints the number of events in each bin
            print(bin_mtautau_dist)
            # Adds the bin distribution to the term distribution
            # Note: the distribution is corrected by the weight, which means that
            #       the final distribution is in pb
            terms_dists[term] += (cross_section / num_events) * bin_mtautau_dist

    # Saves the distribution in a .json file
    with open(f"{folder_path}/mU1_4TeV.json", "w") as file_:
        terms_dists = {term: dist.tolist() for term, dist in terms_dists.items()}
        json.dump(terms_dists, file_, indent=4)
