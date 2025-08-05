import ROOT
from optparse import OptionParser
import datetime
import math
from array import array

#Command line 
parser = OptionParser()
parser.add_option("--input_profile", default="/afs/cern.ch/user/s/singhsh/PixelMonitoring/data/radiation_simulation/profiles/per_phase/BPix_BmI_SEC1_LYR1/profile_BPix_BmI_SEC1_LYR1_phase1.txt", help="Input profile file name, should have been made using PixelMonitoring repository")
parser.add_option("--output_root_file", default="SheffieldFile.root")
parser.add_option("--timestep", default=1, help="step size is 1 second, do not change this!")
parser.add_option("--userTrefC", type="float", default=0.)
parser.add_option("--bandGap", type="float", default=1.21)
parser.add_option("--sensorThickness", type="float", default=0.0285)
parser.add_option("--sensorWidth", type="float", default=6.48)
parser.add_option("--sensorLength", type="float", default=1.62)
parser.add_option("--Ndonor_0", type="float", default=1.7e12)


(opt, args) = parser.parse_args()

# Defining constants used in the analysis
boltzmanConstant = 8.6173303e-5
absoluteZero = 273.15
Tref = 293.15  # check if this is correct??????
Ei = 1.21


# set a reference temperature for the volume corrected leakage current plot (only affects a couple plots)
userTrefK = absoluteZero + float(opt.userTrefC);

#Recreating this file so we don't end up with a bunch of graphs from old runs
ROOT.gROOT.SetBatch(ROOT.kTRUE)
file = ROOT.TFile(opt.output_root_file, "RECREATE");
file.Close()

#Information about the conditions for each data point
# leakage current constants will be stored in this struct
class DataElement:
  def __init__(self):
    fill = 0
    timestamp = 0
    duration = 0
    temperature = 0.
    doseRate = 0
    coolPipeTemp = 0.
    leakageCurrentData = 0.

#From the paper
class LeakageCurrentConstants:
    def __init__(self):

        self.alpha_Tref = 4.81e-17 #A/cm  +/- 0.13 (uncertainty quoted from table 10)
        self.Ak = [0.42, 0.10, 0.23, 0.21, 0.04]  # A/cm # we have 5 values in paper doest explain why???
        self.tauk = [1.2e6, 4.1e4, 3.7e3, 124, 8]  # minutes

class Sensor:
  smallNumber = 1e-30;
  electronCharge = 1.6021766208e-13
  epsilon0 = 8.854187817E-12
  epsilon = 11.68

  # Room temperature (useful for scaling of the tau and Theta functions in eq. 22 of https://cds.cern.ch/record/2773267/)
  Tref_leakage = 294.15

  # According to page 5 of https://iopscience.iop.org/article/10.1088/1748-0221/14/06/P06012/pdf, t0 is one minute = 60 seconds
  t0 = 60.

  # This parameter is multiplied to all fluence rates from the input profile
  DoseRateScaling = 0.0859955711871/8.9275E-02;
  def __init__(self, leakageCurrentConstants,  thickness, width, length, Ndonor_0, userTrefK):
        self.time = 0.0;
        self.totalDose = 0.;
        # IBL // initial donor concentration (right now the code only works for originally n-type bulk material!)Taken from hamburg model
    
        
        self.userTrefK = userTrefK
        self.leakageCurrentConstants = leakageCurrentConstants

        self.depth = thickness
        self.width = width
        self.length = length
        
        self.volume=self.depth*self.width*self.length;

        self.fill_vector = [];
        self.tmpTimeVector = [];
        self.fluence_vector = [];
        self.flux_vector = [];
        self.doseRate_vector = [];
        self.duration_vector = [];

        self.leakage_current = [];
        self.leakage_current_per_module = [];
        self.leakage_current_per_volume = [];
        self.alpha_vector = [];
        self.powerconsumption = [];
        self.temperature_vector = [];

        self.tmpTime_vector_data = [];
        self.fluence_vector_data = [];
        self.fill_vector_data = [];
        self.leakageCurrentData = [];
        self.leakage_current_Tref = [];
        self.Theta_vector = [];
        self.alpha1_vector = [];
        self.alpha2_vector = [];
        self.alpha3_vector = [];
        self.alpha4_vector = [];
        self.alpha5_vector = [];
        self.alpha_total_vector = [];


#Getting leakage current per module and and per volume 
def getPerModule(self, leakageCurrent):
    #convert from surface to volume normalisation
     leakageCurrent/=self.depth;

    #scale to user defined temperature -- are we sure we shouldn't use the input temperature???
     volumeNorm = self.userTrefK*self.userTrefK/(self.Tref_leakage*self.Tref_leakage)*math.exp(-(opt.bandGap/(2*self.boltzmanConstant))*(1/self.userTrefK - 1/self.Tref_leakage))
     return volumeNorm*leakageCurrent

def getPerVolume(self, leakageCurrent):
    #convert from surface to volume normalisation
     leakageCurrent*=sensor.width*sensor.length

    #scale to user defined temperature
     unitSurface = self.userTrefK*self.userTrefK/(self.Tref_leakage*self.Tref_leakage)*math.exp(-(opt.bandGap/(2*self.boltzmanConstant))*(1/self.userTrefK - 1/self.Tref_leakage));

     return unitSurface*leakageCurrent

 

# Sheffield model Formulas

def theta(T, Tref, Ei):
    return math.exp(-Ei / boltzmanConstant * (1. / T - 1. / Tref))

#function to read in the temp/rad profile
def leakageCurrentScale(leakageCurrent, T, Tref, bandGap):
    # Taken from equation 25 of https://cds.cern.ch/record/2773267
    return leakageCurrent*(T*T/(Tref*Tref)*math.exp(-(bandGap/(2*boltzmanConstant))*(1/T - 1/Tref)));

def compute_alpha_sheffield(self, temperature_profile):
    alpha_Tref = self.leakageCurrentConstants.alpha_Tref
    Ak = self.leakageCurrentConstants.Ak
    tauk = [tau * 60 for tau in self.leakageCurrentConstants.tauk]
    cumulative_theta_time = 0.0
    for i in range(len(self.duration_vector)):
        T = self.temperature_vector[i]
        t = self.duration_vector[i]
        doseRate = self.doseRate_vector[i]

        theta_i = theta(T, Tref, Ei)
        cumulative_theta_time += theta_i * t

        total_alpha = 0.0
        for k in range(5):
            tau_k = tauk[k]
            A_k = Ak[k]
            factor = alpha_Tref * A_k * tau_k * theta_i * t
            decay = 1 - math.exp(-theta_i * t / tau_k)
            annealing_decay = math.exp(-cumulative_theta_time / tau_k)

            alpha_k = factor * decay * annealing_decay
            total_alpha += alpha_k
            alpha_k_list.append(alpha_k)
             # Store individual alpha_k values
            self.alpha1_vector.append(alpha_k_list[0])
            self.alpha2_vector.append(alpha_k_list[1])
            self.alpha3_vector.append(alpha_k_list[2])
            self.alpha4_vector.append(alpha_k_list[3])
            self.alpha5_vector.append(alpha_k_list[4])
            self.alpha_total_vector.append(total_alpha)

            fluence_i = doseRate * t  # [neq/cm²]
            Ileak = total_alpha * fluence_i * self.volume

            self.leakage_current.append(Ileak)
            self.leakage_current_per_module.append(self.getPerModule(Ileak))
            self.leakage_current_per_volume.append(self.getPerVolume(Ileak))
            self.leakage_current_Tref.append(leakageCurrentScale(Ileak, T, Tref, opt.bandGap))

def irradiate(self, profile):

    # Increment time and dosage
    # Converting to time in days
      self.time += profile.duration / (24.0 * 3600.0);

      self.totalDose += float(profile.doseRate) * float(profile.duration);

    # Append information given by profile
      self.tmpTimeVector.append(self.time);
      self.fluence_vector.append(self.totalDose);
      self.doseRate_vector.append(profile.doseRate);
      self.duration_vector.append(profile.duration);

      self.temperature_vector.append(profile.temperature);
      self.fill_vector.append(profile.fill);

    # Only fill these ones if there is any data, not just for every time step
      if (profile.leakageCurrentData > 0. and profile.doseRate > 0.):
          print(profile.fill, self.totalDose, userTrefK, profile.temperature, profile.leakageCurrentData, leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap))
          self.fill_vector_data.append(profile.fill);
          self.tmpTime_vector_data.append(self.time);
          self.fluence_vector_data.append(self.totalDose);
          self.leakageCurrentData.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap));
          self.flux_vector.append(profile.doseRate);


      G_i_tmp = 0;
    # TODO: not sure if there is a better name for G_i...
      for i in range(len(self.beta_contribution_vector)):
          if i>0 :
              G_i_tmp += int(self.doseRate_vector[i]) * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1]));
          else:
              G_i_tmp += int(self.doseRate_vector[i]) * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1]))
              self.G_i.append(G_i_tmp);
              alpha = G_i_tmp / (self.totalDose + profile.doseRate);
              if (alpha > 5e-16):
                  self.alpha_vec.append(1e-17);
              else:
                  self.alpha_vec.append(alpha);
   # Convert current in Amps to mA
      amps_to_mA = 1000.
      self.leakage_current.append(G_i_tmp * amps_to_mA * self.depth);
      self.powerconsumption.append(self.leakage_current[-1] * self.NtoV_function());

    # TODO: may need to fix these functions, not sure
    # Convert to per volume or per module units
      self.leakage_current_per_volume.append(self.getPerVolume(self.leakage_current[-1]))
      self.leakage_current_per_module.append(self.getPerModule(self.leakage_current[-1]))


# An array of DataElements that contain the conditions at each time step
def getProfile(filename, sensor):
    profile = [];

    myfile = open(filename, "r")

    print( "Reading temperature/radiation file: ", filename)

    for index, line in enumerate( myfile):
      #Ignoring the first line of the file, which just has some extra info
      if index == 0:
        continue
      words = line.split()
      fill = int(words[0])
      timestamp = int(words[1])
      timestep = int(words[2])
      temperature = float(words[3])
      doseRate = int(words[4])
      leakageCurrent = float(words[5])

      profileSnapshot = DataElement();
      profileSnapshot.fill = fill;
      profileSnapshot.timestamp = timestamp;
      profileSnapshot.duration = timestep;
      if(timestep < 0):
        print("The timestep is negative, which should not happen. Check your inputs")
        print(line)
        profileSnapshot.duration = 0.001
      profileSnapshot.doseRate = doseRate * sensor.DoseRateScaling;
      profileSnapshot.temperature = temperature;   #Add 3 degree?????????????????
      # Scaling to milliamps from microamps
      #profileSnapshot.leakageCurrentData = leakageCurrentScale(leakageCurrent/1.e3, userTrefK, profileSnapshot.temperature, opt.bandGap);
      profileSnapshot.leakageCurrentData = leakageCurrent/1.e3

      profile.append(profileSnapshot);

    myfile.close();

    return profile;

def getBeginTime(profile):
    timestamp = profile[0].timestamp;
    beginTime = datetime.datetime.fromtimestamp(timestamp);

    return beginTime;

# Converts the datetime array into a format that can be used for TGraphs for a nice x-axis
def convertDatetime(dateVector, beginTime):
  reformattedVec = []

  for k in range(len(dateVector)):
    # Convert to seconds based on the date
    daysInSeconds = 86400
    d2 = beginTime + datetime.timedelta(seconds=dateVector[k]*daysInSeconds)
    da = ROOT.TDatime(d2.date().year, d2.date().month, d2.date().day, d2.time().hour, d2.time().minute, d2.time().second);
    reformattedVec.append(da.Convert());

  return reformattedVec


#function to plot graphs containing the different inforamtion
def plot_vectors(vecx, vecy, yName, xName, plotname, mode, rootOutputFileName):
    # We have already created this root file at the start of this script,
    # and so we just want to add on to it
    file = ROOT.TFile(rootOutputFileName, "UPDATE");
    file.cd();

    # Hack for converting to arrays for ROOT TGraph (lists do not work)
    t_timevector = array('d')
    t_vdepvector = array('d')
    for i in vecx:
      t_timevector.append(i)
    for i in vecy:
      t_vdepvector.append(i)
    gr = ROOT.TGraph(len(vecx),t_timevector, t_vdepvector);

    #in case of plot as a function of time, convert the days from simulation into a more handy date
    if(mode=="date") :
        gr.GetXaxis().SetTimeDisplay(1);
        gr.GetXaxis().SetNdivisions(6, 2, 0);
        gr.GetXaxis().SetTimeFormat("%d/%m/%Y");
        gr.GetXaxis().SetTimeOffset(0, "gmt");

    gr.GetXaxis().SetTitle(xName);
    gr.GetYaxis().SetTitle(yName);

    gr.SetName(plotname);
    gr.Write();

    print("Writing ", plotname, " to file ", rootOutputFileName)
    file.Close();

# Read constants values from a config file 

leakageCurrentConstants= LeakageCurrentConstants()

sensor = Sensor(leakageCurrentConstants, thickness=opt.sensorThickness, width=opt.sensorWidth, length = opt.sensorLength, Ndonor_0 = opt.Ndonor_0, userTrefK = userTrefK);

# iterate through the profile and irradiate the sensor at each step
profile = getProfile(opt.input_profile, sensor);
print("Profile succesfully read. Length of the profile is: ", len(profile))
#for t in range(len(profile)):
#    sensor.irradiate(profile[t])

#data output, plotting and visualisation
print("Processing finished, writing data...")

beginTime = getBeginTime(profile);
time_vector = convertDatetime(sensor.tmpTimeVector, beginTime)

print(sensor.alpha_vector)
print(time_vector)
print(leakage_current)
# plots as function of time

plot_vectors(time_vector, sensor.alpha_vector, "#alpha [A/cm]", "Date [days]", "alpha", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_module, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_volume, "I_{leak} (@%d C) [mA/cm^{3}]", "Date [days]", "I_leak_volume", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.G_i, "G_{i} [A/cm^{3}]", "Date [days]", "G_i", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.powerconsumption, "P [mW/cm^{2}]", "Date [days]", "powerconsumption", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.temperature_vector, "Temperature [K]", "Date [days]", "temperature", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.fluence_vector, "Fluence [n_{eq}/cm^{2}]", "Date [days]", "fluence", "date", opt.output_root_file)

plot_vectors(time_vector_data, sensor.flux_vector, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_data", "date", opt.output_root_file)



#plots as function of dose
plot_vectors(sensor.fluence_vector, sensor.Neff_vector, "N_{eff} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", "Neff_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.N_donor_vec, "N_{donor} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", "Ndonors_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.Nacceptor_vector, "N_{acceptor} [1/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", "Nacceptors_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.V_dep_vector, "U_{dep} [V]", "Fluence [n_{eq}/cm^{2}]", "U_dep_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.alpha_vec, "#alpha [A/cm]", "Fluence [n_{eq}/cm^{2}]", "alpha_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_module, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Date [days]", "I_leak_per_module_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_volume, "I_{leak} (@%d C) [mA/cm^{3}]"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_volume_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.G_i, "G_{i} [A/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", "G_i_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.powerconsumption, "P [mW/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "powerconsumption_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.temperature_vector, "Temperature [K]", "Fluence [n_{eq}/cm^{2}", "temperature_vs_fluence", "fluence", opt.output_root_file)

plot_vectors(sensor.fluence_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_per_module_data_vs_fluence", "fluence", opt.output_root_file)

plot_vectors(sensor.fill_vector, sensor.leakage_current, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Fill", "I_leak_vs_fill", "fill", opt.output_root_file)

print(opt.output_root_file, " has been created")
