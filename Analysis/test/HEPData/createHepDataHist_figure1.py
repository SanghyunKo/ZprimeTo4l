import sys
filenames = sys.argv[1:]

import hepdata_lib

class HepDataTable:
    def __init__(self, name):
        self.filename = name
        self.xaxis = "$\\mathrm{M}_{4\\ell}$"

        if "mergedEle" in self.filename:
            self.bkgs = ["out_bkg_mergedEle3E_OS", "out_bkg_mergedEle3E_SS"]
            self.legends = ["Prompt $\\mathrm{e}_{\\mathrm{ME}}$ misID",
                            "Nonprompt $\\mathrm{e}_{\\mathrm{ME}}$ misID"]
            self.figure = "figure1_2"
        elif "me2e" in self.filename:
            self.bkgs = ["out_bkg_mergedEle2E_OS"]
            self.legends = ["Prompt $\\mathrm{e}_{\\mathrm{ME}}$ misID"]
            self.figure = "figure1_3"
        elif "cme" in self.filename:
            self.bkgs = ["out_bkg_mergedEMu1M_MM"]
            self.legends = ["$\\mu_{\\mathrm{MM}}$ misID"]
            self.figure = "figure1_5"
            self.xaxis = "$\\mathrm{M}_{\\mathrm{T}}$"
        elif "mergedMu" in self.filename:
            self.bkgs = ["out_bkg_mergedMu3M_MM"]
            self.legends = ["$\\mu_{\\mathrm{MM}}$ misID"]
            self.figure = "figure1_4"
            self.xaxis = "$\\mathrm{M}_{\\mathrm{T}}$"
        elif "resolved" in self.filename:
            self.bkgs = ["out_bkg_resolvedEle_ZZ",
                         "out_bkg_resolvedEle_Nonprompt",
                         "out_bkg_resolvedEle_NonpromptDR03"]
            self.legends = ["ZZ",
                            "Resolved $\\ell$ misID ($\\Delta\\mathrm{R} > 0.3$)",
                            "Resolved $\\ell$ misID ($\\Delta\\mathrm{R} < 0.3$)"]
            self.figure = "figure1_1"

    def createTable(self):
        self.table = hepdata_lib.Table(self.figure)
        self.table.description = ""
        self.table.keywords["observables"] = ["$\\mathrm{N}_{\\mathrm{events}}/\\mathrm{GeV}$"]

        reader = hepdata_lib.RootFileReader(self.filename)

        # read common histograms
        totb = reader.read_graph("Total background")
        data = reader.read_hist_1d("Data")
        sig = reader.read_hist_1d("signal")

        # read bkgs
        bkgs = []

        for bkg in self.bkgs:
            bkgs.append(reader.read_hist_1d(bkg))

        # create variables and uncertainties
        # x-axis
        xaxis = hepdata_lib.Variable(self.xaxis, is_independent=True, is_binned=True, units="GeV")
        xaxis.values = data["x_edges"]
        # y-axis: N events/GeV
        var_data = hepdata_lib.Variable("Data", is_independent=False, is_binned=False, units="/GeV")
        var_data.values = data["y"]
        # common components
        unc_data = hepdata_lib.Uncertainty("stat", is_symmetric=True)
        unc_data.values = data["dy"]
        var_data.add_uncertainty(unc_data)

        var_totb = hepdata_lib.Variable("Total background", is_independent=False, is_binned=False, units="/GeV")
        var_totb.values = totb["y"]

        unc_totb = hepdata_lib.Uncertainty("stat+syst", is_symmetric=False)
        unc_totb.values = totb["dy"]
        var_totb.add_uncertainty(unc_totb)

        var_sig = hepdata_lib.Variable("Signal", is_independent=False, is_binned=False, units="/GeV")
        var_sig.values = sig["y"]

        unc_sig = hepdata_lib.Uncertainty("stat", is_symmetric=True)
        unc_sig.values = sig["dy"]
        var_sig.add_uncertainty(unc_sig)

        self.table.add_variable(xaxis)
        self.table.add_variable(var_data)
        self.table.add_variable(var_totb)
        self.table.add_variable(var_sig)

        bkg_variables = []
        for bkg, legend in zip(bkgs, self.legends):
            var_bkg = hepdata_lib.Variable(legend, is_independent=False, is_binned=False, units="/GeV")
            var_bkg.values = bkg["y"]

            unc_bkg = hepdata_lib.Uncertainty("stat+syst", is_symmetric=True)
            unc_bkg.values = bkg["dy"]
            var_bkg.add_uncertainty(unc_bkg)

            bkg_variables.append(var_bkg)
            self.table.add_variable(var_bkg)

from hepdata_lib import Submission
submission = Submission()

for filename in filenames:
    table = HepDataTable(filename)
    table.createTable()
    submission.add_table(table.table)

submission.create_files("figure1",remove_old=True)
