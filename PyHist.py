class PyHist:
    def __init__(self, histo):
        """
        Initializes the class PyHist object.
        Parameters:
        - histo (TH1): The input ROOT histogram
        """
        self.histo_name = histo.GetName()
        self.bin_values = [histo.GetBinContent(i) for i in range(1, histo.GetNbinsX() + 1)]
        self.bin_error_low = [histo.GetBinErrorLow(i) for i in range(1, histo.GetNbinsX() + 1)]
        self.bin_error_hi = [histo.GetBinErrorUp(i) for i in range(1, histo.GetNbinsX() + 1)]
        self.bin_edges = [histo.GetBinLowEdge(i) for i in range(1, histo.GetNbinsX() + 2)]
        self.bin_widths = [histo.GetBinWidth(i) for i in range(1, histo.GetNbinsX() + 1)]
        self.is_normalized_by_width = False

    def divide_by_bin_width(self):
        """
        Divides the bin values and errors by the bin widths.
        """
        if not self.is_normalized_by_width:
            self.bin_values = [val / width for val, width in zip(self.bin_values, self.bin_widths)]
            self.bin_error_low = [err / width for err, width in zip(self.bin_error_low, self.bin_widths)]
            self.bin_error_hi = [err / width for err, width in zip(self.bin_error_hi, self.bin_widths)]
            self.is_normalized_by_width = True
        else:
            print("Already normalized by bin width.")

    def get_error_pairs(self):
        """
        For use with pyplot.errorbar.
        """
        return [self.bin_error_low, self.bin_error_hi]
    
    def get_bin_centers(self):
        """
        Calculates the bin centers for the histogram.
        Returns:
        - List of bin centers.
        """
        return [(self.bin_edges[i] + self.bin_edges[i + 1]) / 2 for i in range(len(self.bin_edges) - 1)]

