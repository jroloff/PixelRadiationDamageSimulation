
import ROOT
from optparse import OptionParser
import datetime
import math
from array import array

boltzmanConstant = 8.6173303e-5
absoluteZero = 273.15

#Command line 
parser = OptionParser()
parser.add_option("--input_profile", default="input_profile.txt")
parser.add_option("--output_root_file", default="testFile.root")
parser.add_option("--userTrefC", type="float", default=0.)
parser.add_option("--bandGap", type="float", default=1.21)
parser.add_option("--sensorThickness", type="float", default=0.0285)
parser.add_option("--sensorWidth", type="float", default=6.48)
parser.add_option("--sensorLength", type="float", default=1.62)
parser.add_option("--Ndonor_0", type="float", default=1.7e12)
(opt, args) = parser.parse_args()

userTrefK = absoluteZero + opt.userTrefC

# Sheffield model helpers
def theta(T, Tref, Ei):
    return math.exp(-Ei / boltzmanConstant * (1. / T - 1. / Tref))

def tau(T, Ei, tau0):
    return tau0 * math.exp(Ei / (boltzmanConstant * T))

# Data holder
class DataElement:
    def __init__(self):
        self.fill = 0
        self.timestamp = 0
        self.duration = 0
        self.temperature = 0.
        self.doseRate = 0
        self.leakageCurrentData = 0.

# Leakage constants for Sheffield model
class LeakageCurrentConstants:
    def __init__(self):
        self.alpha_ref = 4.0e-17
        self.Ak = 1.0
        self.Tk = 293.15
        self.Ei = 1.11
        self.tau0 = 1e6

# Sensor object
class Sensor:
    def __init__(self, consts, thickness, width, length, Ndonor_0, userTrefK):
        self.consts = consts
        self.depth = thickness
        self.width = width
        self.length = length
        self.volume = self.depth * self.width * self.length
        self.Ndonor = Ndonor_0
        self.userTrefK = userTrefK

        self.time = 0.
        self.totalDose = 0.
        self.temperature_vector = []
        self.duration_vector = []
        self.doseRate_vector = []
        self.fluence_vector = []
        self.alpha_vec = []
        self.leakage_current = []

    def irradiate(self, profile):
        self.time += profile.duration / (24. * 3600.)
        self.totalDose += profile.doseRate * profile.duration

        self.temperature_vector.append(profile.temperature)
        self.duration_vector.append(profile.duration)
        self.doseRate_vector.append(profile.doseRate)
        self.fluence_vector.append(self.totalDose)

        T_i = profile.temperature
        t_i = profile.duration
        θ_i = theta(T_i, self.consts.Tk, self.consts.Ei)
        τ_k = tau(T_i, self.consts.Ei, self.consts.tau0)
        sum_theta_t = sum([theta(T, self.consts.Tk, self.consts.Ei) * t
                           for T, t in zip(self.temperature_vector, self.duration_vector)])

        α_sheffield = (self.consts.alpha_ref * (self.consts.Ak * self.consts.Tk) /
                       (θ_i * t_i) * (1 - math.exp(-θ_i * t_i / τ_k)) *
                       math.exp(-sum_theta_t / τ_k))

        self.alpha_vec.append(α_sheffield)
        current = profile.doseRate * t_i * α_sheffield
        self.leakage_current.append(current * 1000. * self.depth)  # convert A to mA and scale

# Profile loader
def getProfile(filename):
    profile = []
    with open(filename, "r") as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            words = line.split()
            d = DataElement()
            d.fill = int(words[0])
            d.timestamp = int(words[1])
            d.duration = int(words[2])
            d.temperature = float(words[3])
            d.doseRate = float(words[4])
            d.leakageCurrentData = float(words[5])
            profile.append(d)
    return profile

# Plotter
def plot(vecx, vecy, title, yaxis, name, rootfile):
    file = ROOT.TFile(rootfile, "UPDATE")
    t_x = array('d', vecx)
    t_y = array('d', vecy)
    gr = ROOT.TGraph(len(vecx), t_x, t_y)
    gr.SetName(name)
    gr.GetXaxis().SetTitle("Time [days]")
    gr.GetYaxis().SetTitle(yaxis)
    gr.Write()
    file.Close()
    print(f"Saved plot {name} to {rootfile}")

# Main execution
ROOT.gROOT.SetBatch(ROOT.kTRUE)
file = ROOT.TFile(opt.output_root_file, "RECREATE")
file.Close()

consts = LeakageCurrentConstants()
sensor = Sensor(consts, opt.sensorThickness, opt.sensorWidth, opt.sensorLength, opt.Ndonor_0, userTrefK)

profile = getProfile(opt.input_profile)
for p in profile:
    sensor.irradiate(p)

plot(sensor.fluence_vector, sensor.alpha_vec, "Alpha vs Fluence", "#alpha [A/cm]", "alpha_vs_fluence", opt.output_root_file)
plot(sensor.fluence_vector, sensor.leakage_current, "Ileak vs Fluence", "I_{leak} [mA]", "Ileak_vs_fluence", opt.output_root_file)

print(f"{opt.output_root_file} has been created with Sheffield model simulation results.")

