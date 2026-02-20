import sys
filenames = sys.argv[1:]

import hepdata_lib

class HepDataTable:
    def __init__(self, name):
        self.filename = name
        self.xaxis = "$\\mathrm{M}_{\\mathrm{X}}$"
        self.yaxis = "$\\mathrm{M}_{\\mathrm{Y}}$"

    def createTable(self,treename):
        if treename=="limitTreeA_el":
            self.figure = "figureA_1"
        elif treename=="limitTreeA_mu":
            self.figure = "figureA_2"
        elif treename=="limitTreeZ_el":
            self.figure = "figureA_3"
        elif treename=="limitTreeZ_mu":
            self.figure = "figureA_4"

        self.table = hepdata_lib.Table(self.figure)
        self.table.description = ""
        self.table.keywords["observables"] = ["Upper limit on cross section times branching fraction"]

        reader = hepdata_lib.RootFileReader(self.filename)

        # read tree
        # x-axis
        xaxis = hepdata_lib.Variable(self.xaxis, is_independent=True, is_binned=False, units="GeV")
        xaxis.values = reader.read_tree(treename, "mX")
        # y-axis
        yaxis = hepdata_lib.Variable(self.yaxis, is_independent=True, is_binned=False, units="GeV")
        yaxis.values = reader.read_tree(treename, "mY")
        # observed limit
        obs = hepdata_lib.Variable("Observed limit", is_independent=False, is_binned=False, units="fb")
        obs.values = reader.read_tree(treename, "limit")
        # significance
        significance = hepdata_lib.Variable("Local significance", is_independent=False, is_binned=False, units="")
        significance.values = reader.read_tree(treename, "significance")
        # expected limit
        exp = hepdata_lib.Variable("Expected limit", is_independent=False, is_binned=False, units="fb")
        exp.values = reader.read_tree(treename, "nominal")
        # expected limit uncertainties
        # 1 sigma
        exp_unc_1sigma = hepdata_lib.Uncertainty("68% expected", is_symmetric=False)
        q16 = reader.read_tree(treename, "quant0p16")
        q84 = reader.read_tree(treename, "quant0p84")
        exp_unc_1sigma.values = list(zip(q16, q84))
        exp.add_uncertainty(exp_unc_1sigma)
        # 2 sigma
        exp_unc_2sigma = hepdata_lib.Uncertainty("95% expected", is_symmetric=False)
        q025 = reader.read_tree(treename, "quant0p025")
        q975 = reader.read_tree(treename, "quant0p975")
        exp_unc_2sigma.values = list(zip(q025, q975))
        exp.add_uncertainty(exp_unc_2sigma)

        self.table.add_variable(xaxis)
        self.table.add_variable(yaxis)
        self.table.add_variable(obs)
        self.table.add_variable(significance)
        self.table.add_variable(exp)

from hepdata_lib import Submission
submission = Submission()

table = HepDataTable(filenames[0])
table.createTable("limitTreeA_el")
submission.add_table(table.table)
table.createTable("limitTreeA_mu")
submission.add_table(table.table)
table.createTable("limitTreeZ_el")
submission.add_table(table.table)
table.createTable("limitTreeZ_mu")
submission.add_table(table.table)

submission.create_files("figureA",remove_old=True)
