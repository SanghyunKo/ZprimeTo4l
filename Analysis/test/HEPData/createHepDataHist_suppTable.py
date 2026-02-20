import sys
filenames = sys.argv[1:]

import hepdata_lib

class HepDataTable:
    def __init__(self, name):
        self.filename = name
        self.xaxis = "$\\mathrm{M}_{\\mathrm{X}}$"
        self.yaxis = "$\\mathrm{M}_{\\mathrm{Y}}$"

        if "el" in self.filename:
            self.figure = "data15"
        elif "mu" in self.filename:
            self.figure = "data16"
        elif "A" in self.filename:
            if "REFF" in self.filename:
                self.figure = "data1"
            elif "REMuFF" in self.filename:
                self.figure = "data2"
            elif "RMFF" in self.filename:
                self.figure = "data3"
            elif "MEMu2M" in self.filename:
                self.figure = "data4"
            elif "MMFF" in self.filename:
                self.figure = "data5"
            elif "ME2E" in self.filename:
                self.figure = "data6"
            elif "MEMu1M" in self.filename:
                self.figure = "data7"
            else:
                raise ValueError(f"Unknown filename: {self.filename}")
        elif "Z" in self.filename:
            if "REFF" in self.filename:
                self.figure = "data8"
            elif "REMuFF" in self.filename:
                self.figure = "data9"
            elif "RMFF" in self.filename:
                self.figure = "data10"
            elif "ME3E" in self.filename:
                self.figure = "data11"
            elif "MEMu2M" in self.filename:
                self.figure = "data12"
            elif "MM2E" in self.filename:
                self.figure = "data13"
            elif "MMFF" in self.filename:
                self.figure = "data14"
            else:
                raise ValueError(f"Unknown filename: {self.filename}")
        else:
            raise ValueError(f"Unknown filename: {self.filename}")

    def createTable(self,treename):
        self.table = hepdata_lib.Table(self.figure)
        self.table.description = ""
        self.table.keywords["observables"] = ["Efficiency"]

        reader = hepdata_lib.RootFileReader(self.filename)

        # read tree
        # x-axis
        xaxis = hepdata_lib.Variable(self.xaxis, is_independent=True, is_binned=False, units="GeV")
        xaxis.values = reader.read_tree(treename, "mX")
        # log y-axis
        logmY = hepdata_lib.Variable("log10($\\mathrm{M}_{\\mathrm{Y}}$ [GeV])", is_independent=True, is_binned=False, units="")
        logmY.values = reader.read_tree(treename, "logmY")
        # y-axis
        yaxis = hepdata_lib.Variable(self.yaxis, is_independent=True, is_binned=False, units="GeV")
        yaxis.values = reader.read_tree(treename, "mY")
        # efficiency
        efficiency = hepdata_lib.Variable("Efficiency", is_independent=False, is_binned=False, units="")
        efficiency.values = reader.read_tree(treename, "eff")

        self.table.add_variable(xaxis)
        self.table.add_variable(logmY)
        self.table.add_variable(yaxis)
        self.table.add_variable(efficiency)

from hepdata_lib import Submission
submission = Submission()

for filename in filenames:
    table = HepDataTable(filename)
    table.createTable("effTree")
    submission.add_table(table.table)

submission.create_files("suppEffTable",remove_old=True)
