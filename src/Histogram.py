"""Interface to construct a binning histogram."""

from abc import ABC, abstractmethod
import numpy as np
from typing import List, Callable, Dict
import copy


class Histogram(ABC):
    """Defines the interface for a histogram class."""

    @abstractmethod
    def update_hist(self, event):
        """Updates the histogram with a given event."""
        raise RuntimeError("Trying to use a method from an abstract class.")

    @abstractmethod
    def __copy__(self):
        """Clones an empty histogram."""
        pass


class BinIndexFinder:
    def __init__(self, bin_edges):
        self.bin_edges = bin_edges

    def find_bin_index(self, observable_value: float) -> int:
        """Finds the respective bin index for the given value of the observable."""
        # Look for the right bin number
        for bin_index in range(len(self.bin_edges) - 1):
            if self.bin_edges[bin_index] <= observable_value < self.bin_edges[bin_index + 1]:
                return bin_index
        return -1


class ObservableHistogram(Histogram, np.ndarray, BinIndexFinder):
    """
    One-dimensional histogram for a given observable.
    It behaves like a numpy array.
    The histogram is updated using the method:

    def update_hist(event)

    where it takes a single event as the argument.
    """

    def __new__(cls, bin_edges: List[float], observable: Callable):
        """
        The two necessary parameters for construction are:

        :param bin_edges: The respective bin edges for the histogram.
        :param observable: A function or callable object that computes
                           the observable for a single Event object.
        """
        # Create an empty histogram
        hist = np.zeros(shape=len(bin_edges) - 1).view(cls)
        # Store the bin_edges and observable as attributes
        hist.bin_edges = bin_edges
        hist.observable = observable
        # Return the histogram
        return hist

    def __init__(self, bin_edges, observable):
        BinIndexFinder.__init__(self, bin_edges=bin_edges)

    def __array_finalize__(self, hist):
        if hist is None:
            return
        # Add the attributes
        self.observable = getattr(hist, "observable", None)
        self.bin_edges = getattr(hist, "bin_edges", None)

    def update_hist(self, event):
        """Updates the histogram using the Event object."""
        # Calculate the observable
        obs_value = self.observable(event)
        # Find the bin index
        bin_index = self.find_bin_index(observable_value=obs_value)
        # Update the histogram if the observable is inside the histogram limits
        if 0 <= bin_index < len(self):
            self[bin_index] += 1

    def __copy__(self):
        """Shallow copy of the current histogram."""
        return self.__new__(self.__class__, bin_edges=self.bin_edges, observable=self.observable)


class CorrelatedHist(Histogram, BinIndexFinder):
    """..."""

    def __init__(self, xobservable, yobservable, bin_edges):
        super().__init__(bin_edges=bin_edges)
        self.xobs = xobservable
        self.yobs = yobservable
        self.bin_sum = np.zeros(len(bin_edges) - 1)

    def update_hist(self, event):
        xobs_value = self.xobs(event)
        yobs_value = self.yobs(event)
        bin_index = self.find_bin_index(xobs_value)
        if 0 <= bin_index < len(self.bin_sum):
            self.bin_sum[bin_index] += yobs_value

    def __copy__(self):
        return self.__class__(xobservable=self.xobs, yobservable=self.yobs, bin_edges=self.bin_edges)

    def merge_hist(self, hist):
        self.bin_sum += hist.bin_sum


class HistogramCompound(Histogram):
    """Stores a set of Histogram objects that must be updated."""

    def __init__(self, histograms: Dict[str, Histogram]):
        # Stores a dictionary whose values represent histogram, and the key is a name used to identify them
        self._hist_dict = histograms

    def update_hist(self, event):
        """Updates all the histograms with the given event."""
        for hist_name in self._hist_dict:
            self._hist_dict[hist_name].update_hist(event=event)

    def get_hist(self, hist_name: str):
        """Returns the Histogram object associated with the key 'hist_name'"""
        if hist_name in self._hist_dict:
            return self._hist_dict[hist_name]

    def __copy__(self):
        """Returns a shallow clone of all histograms."""
        clone_dict = {hist_name: copy.copy(hist) for hist_name, hist in self._hist_dict.items()}
        # Creates a container with now the cloned hists
        return self.__class__(histograms=clone_dict)

