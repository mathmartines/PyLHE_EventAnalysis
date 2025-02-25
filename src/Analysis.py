"""Simple classes responsible for iterating over events and performing event-by-event analysis."""

from typing import List, Callable, Dict
from PyLHE_EventAnalysis.src.Histogram import Histogram
import copy


class EventAnalysis:
    """
    Performs the analysis of a single event.
    Holds information about particle selections and event selection cuts.
    """

    def __init__(self, selection_cuts: List[Callable]):
        """
        :param selection_cuts:
            A list of functions representing the selection cuts an event must satisfy.
            Each function must return True if the event passes the cut and False otherwise.
        """
        self._cuts = selection_cuts

    def launch_analysis(self, event) -> bool:
        """
        Launches the analysis on the event.
        Returns True if the event is selected for analysis, and False otherwise.
        """
        # Apply the event selection cuts
        passed_cuts = all(cut(event) for cut in self._cuts)
        # Return a boolean indicating whether the event was selected
        return passed_cuts


class EventLoop:
    """
    Iterates over all events in an .lhe file
    and manages histogram booking with the selected events.
    """

    def __init__(self, file_reader: Callable, histogram_template: Histogram):
        # Function responsible for reading events
        self._file_reader = file_reader
        # Store the histogram template to be used for constructing histograms
        self._hist_template = histogram_template

    def analyse_events(self, filename: str, event_analyses: Dict[str, EventAnalysis]):
        """
        Runs the analysis on events from the .lhe file and returns a histogram
        constructed from the selected events.

        :param filename: Path to the .lhe file storing the events.
        :param event_analyses: Dictionary with all the diferent analysis that must be applied to the
                               list of events.

        :return: Dict with the booked histogram for each analysis.
        """
        print(f"Reading events from file: {filename}")

        # Count the number of processed events
        evt_number = 0
        # Initialize an empty for each of the analysis
        analyses_hist = {analysis_name: copy.copy(self._hist_template) for analysis_name in event_analyses}

        # Iterate over events in the file
        for event in self._file_reader(filename):
            if evt_number > 0 and evt_number % 10000 == 0:
                print(f"INFO: Processed {evt_number} events")
            # Iterates over all the analyses
            for analysis_name, event_analysis in event_analyses.items():
                # Launch the analysis on the event
                passed_cuts = event_analysis.launch_analysis(event=event)
                # Update the histogram if the event passes selection cuts
                if passed_cuts:
                    analyses_hist[analysis_name].update_hist(event=event)
            # Increment event counter
            evt_number += 1

        # Returns the dictionary with booked histogram for each analysis
        return analyses_hist
