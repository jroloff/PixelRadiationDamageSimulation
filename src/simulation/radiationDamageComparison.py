import ROOT 
from optparse import OptionParser
import datetime
import math
from array import array


parser = OptionParser()
#These are the main options that should be changed
#parser.add_option("--input_profile", default="../../oldPixelMonitoring/data/radiation_simulation/profiles/per_phase/BPix_BmI_SEC1_LYR2/profile_BPix_BmI_SEC1_LYR2_phase1_toy.txt", help="Input profile file name, should have been made using PixelMonitoring repository")
parser.add_option("--input_profile", default="../../oldPixelMonitoring/data/radiation_simulation/profiles/per_year/BPix_BmI_SEC1_LYR2/profile_BPix_BmI_SEC1_LYR2_2017.txt", help="Input profile file name, should have been made using PixelMonitoring repository")

parser.add_option("--output_root_file", default="testFile.root", help="Output ROOT file name")
# These are options that shoudlb e changed if you are using different sensors etc.
#parser.add_option("--timestep", default=1, help="step size is 1 second, do not change this!")
parser.add_option("--userTrefC", type="float",  default=0., help="reference temprature in celsius")
parser.add_option("--bandGap", default=1.21, help="eV used for scaling temperatures")
# Dimensions of sensors for CMS are taken from page 68 of https://cds.cern.ch/record/2773267
parser.add_option("--sensorThickness", default=0.0285, help="Thickness of the sensor in cm. This is correct for CMS")
parser.add_option("--sensorWidth", default=6.48, help="Width of the sensor in cm. 285 is correct for CMS")
parser.add_option("--sensorLength", default=1.62, help="Length of the sensor in cm. 285 is correct for CMS")
parser.add_option("--useBunches", default=0, help="Do the bunches exist in the file?")

(opt, args) = parser.parse_args()


## Defining constants used in the analysis
absoluteZero = 273.15

# set a reference temperature for the volume corrected leakage current plot (only affects a couple plots)
userTrefK = absoluteZero + float(opt.userTrefC); 

#Recreating this file so we don't end up with a bunch of graphs from old runs
ROOT.gROOT.SetBatch(ROOT.kTRUE)
file = ROOT.TFile(opt.output_root_file, "RECREATE");
file.Close()


#function to read in the temp/rad profile
def leakageCurrentScale(leakageCurrent, T, Tref, bandGap, boltzmanConstant):
    # Taken from equation 25 of https://cds.cern.ch/record/2773267
    #return leakageCurrent
    return leakageCurrent*(T*T/(Tref*Tref)*math.exp(-(bandGap/(2*boltzmanConstant))*(1/T - 1/Tref)));

# An array of DataElements that contain the conditions at each time step

def isGoodFill(nBunch, fillLength):
  # Don't want too few bunches?
  if nBunch < 1000:
    return False
  # Arbitrary guess of fill length over 40 minutes
  if fillLength < 10000:
    return False
  return True

#Information about the conditions for each data point
class DataElement:
  def __init__(self):
    fill = 0
    timestamp = 0
    duration = 0
    temperature = 0.
    doseRate = 0
    leakageCurrentData = 0.
    bunches = 0

# leakage current constants will be stored in this struct
# Note: the default values are taken from table 9 of https://cds.cern.ch/record/2773267
# All values can be changed directly when initializing the function, but this is not done by default
class LeakageCurrentConstants:
  def __init__(self, _alpha_1 = 1.23e-17, _alpha_0_star= 7.07e-17, _beta=3.29e-18, _k01=1.2e13, _E1=1.11, _E1_star=1.3):
    #  A/cm  +/- 0.06 (uncertainty quoted from table 9 as well)
    self.alpha_1 = _alpha_1;
    # A/cm
    self.alpha_0_star = _alpha_0_star;
    #A/cm  +/- 0.18
    self.beta = _beta;
    #  /s (Hz)   +5.3/-1.0
    self.k01 = _k01;
    #eV    +/- 0.05
    self.E1 = _E1;
    #  eV    +/- 0.14
    self.E1_star = _E1_star;

    #self.alpha_Tref = 4.81e-17 #A/cm  +/- 0.13 (uncertainty quoted from table 10)
    #This is different than the value in the note, but it is consistent with the value from https://indico.cern.ch/event/162497/contributions/1411637/attachments/190492/267366/111115IEWG_Gibson.pdf
    # when scaled from -7 -> 20 degrees
    self.alpha_Tref = 9.504e-17 #A/cm  +/- 0.13 (uncertainty quoted from table 10)
    self.Ak = [0.42, 0.10, 0.23, 0.21, 0.04]  # A/cm # we have 5 values in paper why???
    self.tauk = [1.2e6*60, 4.1e4*60, 3.7e3*60, 124*60, 8*60]  # minutes * 60 seconds per minute



# all atributes of a sensor (like irradiation, doping concentration, leakage currents, etc) will be stored in this class

class Sensor:
  # Small value -- just something to compare to to make sure numbers are reasonable
  # Epsilon is used for something else...
  smallNumber = 1e-30;

  # Fundamental constants that are relevant for some of these calculations
  boltzmanConstant = 8.6173303e-5


  #T-ref cannot be changed in a trivial way since the other parameters (e.g. alpha*_0) were evaluated at this reference temperature!!!
  # We have 2 different values of Tref based on the difference in the reference temperature for  tables 9 and 10
  # https://cds.cern.ch/record/2773267/
  # These are very close to each other, but this difference is intentional
  # This Tref is useful for scaling of eq. 24 of https://cds.cern.ch/record/2773267/
  # TODO: Actually, isn't this for the Sheffield model???
  Tref = 293.15;
  # Room temperature (useful for scaling of the tau and Theta functions in eq. 22 of https://cds.cern.ch/record/2773267/)
  Tref_leakage = 294.15
  # Different scaling for sheffield...
  Tref_sheff = 293.15

  # According to page 5 of https://iopscience.iop.org/article/10.1088/1748-0221/14/06/P06012/pdf, t0 is one minute = 60 seconds
  t0 = 60.

  # This parameter is multiplied to all fluence rates from the input profile
  # TODO: Not sure where these numbers come from
  # TODO: make this configurable
  # Could this be the scaling factor that is discussed in 5.3.1.3?
  # https://cds.cern.ch/record/2773267/files/10.23731_CYRM-2021-001.59.pdf
  # This number is very close to 1, and it mentioned that this was for bpix
  DoseRateScaling = 0.0859955711871/8.9275E-02;

  def __init__(self, leakageCurrentConstants, thickness, width, length, userTrefK):
        # Things related to conditions
        self.time = 0.0;
        self.totalDose = 0.;
        # IBL       // initial donor concentration (right now the code only works for originally n-type bulk material!)
        # TODO:Figure out how to get the right number
        self.userTrefK = userTrefK
        self.currentTemp = 265

        self.leakageCurrentConstants = leakageCurrentConstants

        self.depth = thickness
        self.width = width
        self.length = length
        # Volume in cm^3
        self.volume=self.depth*self.width*self.length;

        self.fill_vector = [];
        self.tmpTimeVector = [];
        self.fluence_vector = [];
        self.flux_vector = [];
        self.doseRate_vector = [];
        self.duration_vector = [];

 
        self.leakage_current_hamburg = [];
        self.leakage_current_hamburg_first = []
        self.leakage_current_hamburg_second = []
        self.leakage_current_hamburgCutBins_first = []
        self.leakage_current_hamburgCutBins_second = []

        self.leakage_current_hamburg_scale_to_data = []
        self.leakage_current_hamburg_first_scale_to_data = []
        self.leakage_current_hamburg_second_scale_to_data = []
        self.leakage_current_hamburgCutBins_first_scale_to_data = []
        self.leakage_current_hamburgCutBins_second_scale_to_data = []

        self.leakage_current_hamburg_scale_to_dataOldTemp = []
        self.leakage_current_hamburg_first_scale_to_dataOldTemp = []
        self.leakage_current_hamburg_second_scale_to_dataOldTemp = []
        self.leakage_current_hamburgCutBins_first_scale_to_dataOldTemp = []
        self.leakage_current_hamburgCutBins_second_scale_to_dataOldTemp = []

        self.leakage_current_sheffield = [];
        self.leakage_current_sheffield_first = [];
        self.leakage_current_sheffield_second = [];
        self.leakage_current_sheffieldCutBins_first = [];
        self.leakage_current_sheffieldCutBins_second = [];

        self.leakage_current_sheffield_scale_to_data = []
        self.leakage_current_sheffield_first_scale_to_data = []
        self.leakage_current_sheffield_second_scale_to_data = []
        self.leakage_current_sheffieldCutBins_first_scale_to_data = []
        self.leakage_current_sheffieldCutBins_second_scale_to_data = []

        self.leakage_current_sheffield_scale_to_dataOldTemp = []
        self.leakage_current_sheffield_first_scale_to_dataOldTemp = []
        self.leakage_current_sheffield_second_scale_to_dataOldTemp = []
        self.leakage_current_sheffieldCutBins_first_scale_to_dataOldTemp = []
        self.leakage_current_sheffieldCutBins_second_scale_to_dataOldTemp = []

        self.temperature_vector = [];

        # Things for computing simulation behavior
        self.G_i = [];
        self.alpha_vec = [];
        self.Theta_vector = []
        self.alpha1_contribution_vector = []
        self.alpha1_contribution_sum_vector = []
        self.beta_contribution_vector = []
        self.beta_contribution_sum_vector = []

        self.sheffield_alpha0_contribution_vector = []
        self.sheffield_alpha1_contribution_vector = []
        self.sheffield_alpha2_contribution_vector = []
        self.sheffield_alpha3_contribution_vector = []
        self.sheffield_alpha4_contribution_vector = []

        self.sheffield_alpha0_contribution_sum_vector = []
        self.sheffield_alpha1_contribution_sum_vector = []
        self.sheffield_alpha2_contribution_sum_vector = []
        self.sheffield_alpha3_contribution_sum_vector = []
        self.sheffield_alpha4_contribution_sum_vector = []


        # Things for plotting data
        self.tmpTime_vector_data = [];
        self.fluence_vector_data = [];
        self.fill_vector_data = [];
        self.leakageCurrentData = [];

        self.tmpTime_vector_timeStart_data_avgFill = [];
        self.tmpTime_vector_timeEnd_data_avgFill = [];
        self.tmpTime_vector_timeStart_dataCutBins_avgFill = [];
        self.tmpTime_vector_timeEnd_dataCutBins_avgFill = [];

        self.tmpTime_vector_data_avgFill = [];
        self.fluence_vector_data_avgFill = [];
        self.fill_vector_data_avgFill = [];
        self.flux_vector_avgFill = [];
        self.temperature_vector_avgFill = [];
        self.tmpTime_vector_dataCutBins_avgFill = [];
        self.fluence_vector_dataCutBins_avgFill = [];
        self.fill_vector_dataCutBins_avgFill = [];
        self.flux_vectorCutBins_avgFill = [];
        self.temperature_vectorCutBins_avgFill = [];

        self.tmpTime_vector_dataCutBins_avgFill = [];
        self.fluence_vector_dataCutBins_avgFill = [];
        self.fill_vector_dataCutBins_avgFill = [];
        self.flux_vectorCutBins_avgFill = [];
        self.temperature_vectorCutBins_avgFill = [];

        self.leakageCurrentData_avgFill = [];
        self.leakageCurrentDataCutBins_avgFill = [];
        self.leakageCurrentDataNoRescale_avgFill = [];
        self.leakageCurrentDataNoRescaleCutBins_avgFill = [];

        self.fill_vector_data_median = [];
        self.fill_vector_data_first = [];
        self.fill_vector_data_second = [];
        self.fill_vector_dataCutBins_median = [];
        self.fill_vector_dataCutBins_first = [];
        self.fill_vector_dataCutBins_second = [];


        self.tmpTime_vector_data_median = [];
        self.fluence_vector_data_median = [];
        self.flux_vector_median = [];
        self.temperature_vector_median = [];
        self.tmpTime_vector_dataCutBins_median = [];
        self.fluence_vector_dataCutBins_median = [];
        self.flux_vectorCutBins_median = [];
        self.temperature_vectorCutBins_median = [];

        self.leakageCurrentData_median = [];
        self.leakageCurrentDataCutBins_median = [];
        self.leakageCurrentDataNoRescale_median = [];
        self.leakageCurrentDataNoRescaleCutBins_median = [];


        self.tmpTime_vector_data_first = [];
        self.fluence_vector_data_first = [];
        self.flux_vector_first = [];
        self.temperature_vector_first = [];
        self.tmpTime_vector_dataCutBins_first = [];
        self.fluence_vector_dataCutBins_first = [];
        self.flux_vectorCutBins_first = [];
        self.temperature_vectorCutBins_first = [];

        self.leakageCurrentData_first = [];
        self.leakageCurrentDataCutBins_first = [];
        self.leakageCurrentDataNoRescale_first = [];
        self.leakageCurrentDataNoRescaleCutBins_first = [];

        self.tmpTime_vector_data_second = [];
        self.fluence_vector_data_second = [];
        self.flux_vector_second = [];
        self.temperature_vector_second = [];
        self.tmpTime_vector_dataCutBins_second = [];
        self.fluence_vector_dataCutBins_second = [];
        self.flux_vectorCutBins_second = [];
        self.temperature_vectorCutBins_second = [];

        self.leakageCurrentData_second = [];
        self.leakageCurrentDataCutBins_second = [];
        self.leakageCurrentDataNoRescale_second = [];
        self.leakageCurrentDataNoRescaleCutBins_second = [];


        self.tmpTime_vector_data_medianVec = [];
        self.fluence_vector_data_medianVec = [];
        self.flux_vector_medianVec = [];
        self.temperature_vector_medianVec = [];
        self.leakageCurrentData_medianVec = [];
        self.leakageCurrentDataCutBins_medianVec = [];
        self.leakageCurrentDataNoRescale_medianVec = [];
        self.leakageCurrentDataNoRescaleCutBins_medianVec = [];

        self.fillLength = 0
        self.cBunches = 0

  def median(self, vector):
    vector.sort()
    n = len(vector)
    return vector[math.floor(n/2)]


  # TODO: maybe this should be at the current temperature, since we will need to scale the data, not just the predictions????
  def getPerModule(self, leakageCurrent, T, Tref):
    #convert from surface to volume normalisation
    leakageCurrent/=self.depth;

    #scale to user defined temperature -- are we sure we shouldn't use the input temperature???
    volumeNorm = leakageCurrentScale(leakageCurrent, T, Tref, opt.bandGap, self.boltzmanConstant) 

    return volumeNorm

  def getPerVolume(self, leakageCurrent, T, Tref):
    #convert from surface to volume normalisation
    leakageCurrent*=sensor.width*sensor.length
  
    #scale to user defined temperature
    unitSurface = leakageCurrentScale(leakageCurrent, T, Tref, opt.bandGap, self.boltzmanConstant)

    return unitSurface


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

    # Equations taken from page 5 of https://cds.cern.ch/record/2773267
    tauInverse = self.leakageCurrentConstants.k01 * math.exp(-self.leakageCurrentConstants.E1 / (self.boltzmanConstant * int(profile.temperature)))
    # The individual term in the sum over j of t_j / tau(T_j)
    self.alpha1_contribution_vector.append(profile.duration *tauInverse)

    # Equations taken from page 5 of https://cds.cern.ch/record/2773267
    self.Theta_vector.append(math.exp(-self.leakageCurrentConstants.E1_star/self.boltzmanConstant * (1.0 / float(int(profile.temperature)) - 1.0/self.Tref_leakage) ))

    # Equation 3.2 of https://iopscience.iop.org/article/10.1088/1748-0221/14/06/P06012/pdf, t0 is one minute = 60 seconds
    # The individual term in the sum over Theta(T_j)*t_j / t0
    self.beta_contribution_vector.append(self.Theta_vector[-1] * profile.duration / self.t0)

    # The sum of terms from j=0 to n of t_j / tau(T_j)
    if len(self.alpha1_contribution_sum_vector):
      self.alpha1_contribution_sum_vector.append(self.alpha1_contribution_sum_vector[-1] + self.alpha1_contribution_vector[-1] )
    else:
      self.alpha1_contribution_sum_vector.append(self.alpha1_contribution_vector[-1] )
    # The sum of terms from j=0 to n of Theta(T_j)*t_j / t0
    if len(self.beta_contribution_sum_vector):
      self.beta_contribution_sum_vector.append(self.beta_contribution_sum_vector[-1] + self.beta_contribution_vector[-1])
    else:
      self.beta_contribution_sum_vector.append(self.beta_contribution_vector[-1])

    theta2 = math.exp(-opt.bandGap/self.boltzmanConstant * (1.0 / float(int(profile.temperature)) - 1.0/self.Tref_sheff) )
    #theta2 = math.exp(-1.09/self.boltzmanConstant * (1.0 / float(int(profile.temperature)) - 1.0/self.Tref_sheff) )
    self.sheffield_alpha0_contribution_vector.append(theta2 * self.duration_vector[-1] / self.leakageCurrentConstants.tauk[0])
    self.sheffield_alpha1_contribution_vector.append(theta2 * self.duration_vector[-1] / self.leakageCurrentConstants.tauk[1])
    self.sheffield_alpha2_contribution_vector.append(theta2 * self.duration_vector[-1] / self.leakageCurrentConstants.tauk[2])
    self.sheffield_alpha3_contribution_vector.append(theta2 * self.duration_vector[-1] / self.leakageCurrentConstants.tauk[3])
    self.sheffield_alpha4_contribution_vector.append(theta2 * self.duration_vector[-1] / self.leakageCurrentConstants.tauk[4])

    if len(self.sheffield_alpha4_contribution_sum_vector):
      self.sheffield_alpha0_contribution_sum_vector.append(self.sheffield_alpha0_contribution_sum_vector[-1] + self.sheffield_alpha0_contribution_vector[-1])
      self.sheffield_alpha1_contribution_sum_vector.append(self.sheffield_alpha1_contribution_sum_vector[-1] + self.sheffield_alpha1_contribution_vector[-1])
      self.sheffield_alpha2_contribution_sum_vector.append(self.sheffield_alpha2_contribution_sum_vector[-1] + self.sheffield_alpha2_contribution_vector[-1])
      self.sheffield_alpha3_contribution_sum_vector.append(self.sheffield_alpha3_contribution_sum_vector[-1] + self.sheffield_alpha3_contribution_vector[-1])
      self.sheffield_alpha4_contribution_sum_vector.append(self.sheffield_alpha4_contribution_sum_vector[-1] + self.sheffield_alpha4_contribution_vector[-1])
    else:
      self.sheffield_alpha0_contribution_sum_vector.append(self.sheffield_alpha0_contribution_vector[-1])
      self.sheffield_alpha1_contribution_sum_vector.append(self.sheffield_alpha1_contribution_vector[-1])
      self.sheffield_alpha2_contribution_sum_vector.append(self.sheffield_alpha2_contribution_vector[-1])
      self.sheffield_alpha3_contribution_sum_vector.append(self.sheffield_alpha3_contribution_vector[-1])
      self.sheffield_alpha4_contribution_sum_vector.append(self.sheffield_alpha4_contribution_vector[-1])

    G_i_tmp = 0;  
    sheffield_G_i_tmp = 0;  

    sheffield_G_0_tmp = 0;  
    sheffield_G_1_tmp = 0;  
    sheffield_G_2_tmp = 0;  
    sheffield_G_3_tmp = 0;  
    sheffield_G_4_tmp = 0;  

    # TODO: not sure if there is a better name for G_i...
    for i in range(len(self.beta_contribution_vector)): 
        if i>0 :
            # Sum over everything minus sum up to i to give sum from i to n
            G_i_tmp += self.doseRate_vector[i] * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1]));

            sheffield_G_0_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[0] / self.sheffield_alpha0_contribution_vector[-1] * (1-math.exp(-self.sheffield_alpha0_contribution_vector[-1])) * math.exp(-(self.sheffield_alpha0_contribution_sum_vector[-1] - self.sheffield_alpha0_contribution_sum_vector[i-1]))

            sheffield_G_1_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[1] / self.sheffield_alpha1_contribution_vector[-1] * (1-math.exp(-self.sheffield_alpha1_contribution_vector[-1])) * math.exp(-(self.sheffield_alpha1_contribution_sum_vector[-1] - self.sheffield_alpha1_contribution_sum_vector[i-1]))

            sheffield_G_2_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[2] / self.sheffield_alpha2_contribution_vector[-1] * (1-math.exp(-self.sheffield_alpha2_contribution_vector[-1])) * math.exp(-(self.sheffield_alpha2_contribution_sum_vector[-1] - self.sheffield_alpha2_contribution_sum_vector[i-1]))

            sheffield_G_3_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[3] / self.sheffield_alpha3_contribution_vector[-1] * (1-math.exp(-self.sheffield_alpha3_contribution_vector[-1])) * math.exp(-(self.sheffield_alpha3_contribution_sum_vector[-1] - self.sheffield_alpha3_contribution_sum_vector[i-1]))

            sheffield_G_4_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[4] / self.sheffield_alpha4_contribution_vector[-1] * (1-math.exp(-self.sheffield_alpha4_contribution_vector[-1])) * math.exp(-(self.sheffield_alpha4_contribution_sum_vector[-1] - self.sheffield_alpha4_contribution_sum_vector[i-1]))

        else:
            # If we are on the first time increment, then the sum from j=i to i-1 is just the sum from 0 to i-1, so we don't need to do anything special
            G_i_tmp += self.doseRate_vector[i] * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1]))

            sheffield_G_0_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[0] / self.sheffield_alpha0_contribution_vector[i] * (1-math.exp(-self.sheffield_alpha0_contribution_vector[i])) * math.exp(-self.sheffield_alpha0_contribution_sum_vector[-1])
              
            sheffield_G_1_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[1] / self.sheffield_alpha1_contribution_vector[i] * (1-math.exp(-self.sheffield_alpha1_contribution_vector[i])) * math.exp(-self.sheffield_alpha1_contribution_sum_vector[-1])
              
            sheffield_G_2_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[2] / self.sheffield_alpha2_contribution_vector[i] * (1-math.exp(-self.sheffield_alpha2_contribution_vector[i])) * math.exp(-self.sheffield_alpha2_contribution_sum_vector[-1])
    
            sheffield_G_3_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[3] / self.sheffield_alpha3_contribution_vector[i] * (1-math.exp(-self.sheffield_alpha3_contribution_vector[i])) * math.exp(-self.sheffield_alpha3_contribution_sum_vector[-1])
           
            sheffield_G_4_tmp  += self.doseRate_vector[i] * self.duration_vector[i] * self.leakageCurrentConstants.alpha_Tref * self.leakageCurrentConstants.Ak[4] / self.sheffield_alpha4_contribution_vector[i] * (1-math.exp(-self.sheffield_alpha4_contribution_vector[i])) * math.exp(-self.sheffield_alpha4_contribution_sum_vector[-1])


        sheffield_G_i_tmp = sheffield_G_0_tmp + sheffield_G_1_tmp + sheffield_G_2_tmp + sheffield_G_3_tmp + sheffield_G_4_tmp 


    self.G_i.append(G_i_tmp);
    
    # alpha is the current-related damage rate
    # Delta I = alpha * Fluence * volume???
    alpha = G_i_tmp / (self.totalDose + profile.doseRate);

    # Page 102 of https://inspirehep.net/files/9def3c9d10b5a5f1d5c3f249f04e139c mentions a range of validity -- maybe this is related???
    if (alpha > 5e-16):
        self.alpha_vec.append(1e-17);
    else:
        self.alpha_vec.append(alpha);
    

    # Convert current in Amps to mA
    amps_to_mA = 1000.
    self.leakage_current_hamburg.append(G_i_tmp * amps_to_mA * self.depth);
    self.leakage_current_sheffield.append(sheffield_G_i_tmp*amps_to_mA*self.depth);

    # Only fill these ones if there is any data, not just for every time step
    if (profile.leakageCurrentData > 0. and profile.doseRate > 100. and profile.fill > 0):
      tDelta = profile.duration / (24.0 * 3600.0);
      # Undoing a correction to the temperaturs based on fluence
      self.currentTemp = profile.temperature -   3e-8 * profile.doseRate
      self.fill_vector_data.append(profile.fill);
      self.tmpTime_vector_data.append(self.time);
      self.fluence_vector_data.append(self.totalDose);
      self.leakageCurrentData.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
      self.flux_vector.append(profile.doseRate);


      # TODO: may need to fix these functions, not sure
      # Convert to per volume or per module units
      isGood = isGoodFill(self.cBunches, self.fillLength)
      self.leakage_current_hamburg_scale_to_data.append(leakageCurrentScale(self.leakage_current_hamburg[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant));
      self.leakage_current_sheffield_scale_to_data.append(leakageCurrentScale(self.leakage_current_sheffield[-1], self.currentTemp, self.Tref_leakage, opt.bandGap, self.boltzmanConstant));
      self.leakage_current_hamburg_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_hamburg[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant));
      self.leakage_current_sheffield_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_sheffield[-1], profile.temperature, self.Tref_leakage, opt.bandGap, self.boltzmanConstant));


      if(len(self.fill_vector_data_avgFill)==0 or self.fill_vector_data_avgFill[-1] != profile.fill):
        if(len(self.fill_vector_data_avgFill)!=0):
          timeDelta = (self.tmpTime_vector_timeEnd_data_avgFill[-1])
          self.fillLength = (self.time - self.tmpTime_vector_timeStart_data_avgFill[-1])*24*3600.
          isGood = isGoodFill(self.cBunches, self.fillLength)

          if(timeDelta):
            self.tmpTime_vector_data_avgFill[-1] = self.tmpTime_vector_data_avgFill[-1] / timeDelta
            self.fluence_vector_data_avgFill[-1] = self.fluence_vector_data_avgFill[-1]
            self.leakageCurrentData_avgFill[-1] = self.leakageCurrentData_avgFill[-1] / timeDelta
            self.flux_vector_avgFill[-1] = self.flux_vector_avgFill[-1] / timeDelta
            self.temperature_vector_avgFill[-1] = self.temperature_vector_avgFill[-1]/ timeDelta

            self.fill_vector_data_median.append(profile.fill);
            self.tmpTime_vector_data_median.append(self.median(self.tmpTime_vector_data_medianVec));
            self.fluence_vector_data_median.append(self.fluence_vector_data_avgFill[-1]);
            self.leakageCurrentData_median.append(self.median(self.leakageCurrentData_medianVec));
            self.leakageCurrentDataNoRescale_median.append(self.median(self.leakageCurrentDataNoRescale_medianVec));
            self.flux_vector_median.append(self.median(self.flux_vector_medianVec));
            self.temperature_vector_median.append(self.median(self.temperature_vector_medianVec));


            self.fill_vector_dataCutBins_median.append(profile.fill);
            self.tmpTime_vector_dataCutBins_median.append(self.median(self.tmpTime_vector_data_medianVec));
            self.fluence_vector_dataCutBins_median.append(self.fluence_vector_data_avgFill[-1]);
            self.leakageCurrentDataCutBins_median.append(self.median(self.leakageCurrentData_medianVec));
            self.leakageCurrentDataNoRescaleCutBins_median.append(self.median(self.leakageCurrentDataNoRescale_medianVec));
            self.flux_vectorCutBins_median.append(self.median(self.flux_vector_medianVec));
            self.temperature_vectorCutBins_median.append(self.median(self.temperature_vector_medianVec));

            if not isGood or len(self.leakageCurrentDataCutBins_second) != len(self.leakageCurrentDataCutBins_first):
              if len(self.leakageCurrentDataCutBins_second) >= len(self.leakageCurrentDataCutBins_first):
                self.fill_vector_dataCutBins_second.pop()
                self.tmpTime_vector_dataCutBins_second.pop()
                self.fluence_vector_dataCutBins_second.pop()
                self.leakageCurrentDataCutBins_second.pop()
                self.leakageCurrentDataNoRescaleCutBins_second.pop()
                self.flux_vectorCutBins_second.pop()
                self.temperature_vectorCutBins_second.pop()
                self.leakage_current_hamburgCutBins_second_scale_to_data.pop()
                self.leakage_current_sheffieldCutBins_second_scale_to_data.pop()
                self.leakage_current_hamburgCutBins_second_scale_to_dataOldTemp.pop()
                self.leakage_current_sheffieldCutBins_second_scale_to_dataOldTemp.pop()

              if len(self.leakageCurrentDataCutBins_second) < len(self.leakageCurrentDataCutBins_first):
                self.fill_vector_dataCutBins_first.pop()
                self.tmpTime_vector_dataCutBins_first.pop()
                self.fluence_vector_dataCutBins_first.pop()
                self.leakageCurrentDataCutBins_first.pop()
                self.leakageCurrentDataNoRescaleCutBins_first.pop()
                self.flux_vectorCutBins_first.pop()
                self.temperature_vectorCutBins_first.pop()
                self.leakage_current_hamburgCutBins_first_scale_to_data.pop()
                self.leakage_current_sheffieldCutBins_first_scale_to_data.pop()
                self.leakage_current_hamburgCutBins_first_scale_to_dataOldTemp.pop()
                self.leakage_current_sheffieldCutBins_first_scale_to_dataOldTemp.pop()

              if not isGood:
                self.tmpTime_vector_timeStart_dataCutBins_avgFill.pop()
                self.tmpTime_vector_timeEnd_dataCutBins_avgFill.pop()
                self.fill_vector_dataCutBins_avgFill.pop()
                self.tmpTime_vector_dataCutBins_avgFill.pop()
                self.fluence_vector_dataCutBins_avgFill.pop()
                self.leakageCurrentDataCutBins_avgFill.pop()
                self.leakageCurrentDataNoRescaleCutBins_avgFill.pop()
                self.flux_vectorCutBins_avgFill.pop()
                self.temperature_vectorCutBins_avgFill.pop()
  
  
                self.fill_vector_dataCutBins_median.pop()
                self.tmpTime_vector_dataCutBins_median.pop()
                self.fluence_vector_dataCutBins_median.pop()
                self.leakageCurrentDataCutBins_median.pop()
                self.leakageCurrentDataNoRescaleCutBins_median.pop()
                self.flux_vectorCutBins_median.pop()
                self.temperature_vectorCutBins_median.pop()




        self.tmpTime_vector_timeStart_data_avgFill.append(self.time);
        self.tmpTime_vector_timeEnd_data_avgFill.append(tDelta);
        self.fill_vector_data_avgFill.append(profile.fill);
        self.tmpTime_vector_data_avgFill.append(self.time*tDelta);
        self.fluence_vector_data_avgFill.append(self.totalDose);
        self.leakageCurrentData_avgFill.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant)*tDelta);
        self.leakageCurrentDataNoRescale_avgFill.append(profile.leakageCurrentData);
        self.flux_vector_avgFill.append(profile.doseRate*tDelta);
        self.temperature_vector_avgFill.append(profile.temperature*tDelta);

        self.tmpTime_vector_timeStart_dataCutBins_avgFill.append(self.time);
        self.tmpTime_vector_timeEnd_dataCutBins_avgFill.append(tDelta);
        self.fill_vector_dataCutBins_avgFill.append(profile.fill);
        self.tmpTime_vector_dataCutBins_avgFill.append(self.time*tDelta);
        self.fluence_vector_dataCutBins_avgFill.append(self.totalDose);
        self.leakageCurrentDataCutBins_avgFill.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant)*tDelta);
        self.leakageCurrentDataNoRescaleCutBins_avgFill.append(profile.leakageCurrentData);
        self.flux_vectorCutBins_avgFill.append(profile.doseRate*tDelta);
        self.temperature_vectorCutBins_avgFill.append(profile.temperature*tDelta);

        self.fill_vector_data_first.append(profile.fill);
        self.tmpTime_vector_data_first.append(self.time);
        self.fluence_vector_data_first.append(self.totalDose);
        self.leakageCurrentData_first.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.leakageCurrentDataNoRescale_first.append(profile.leakageCurrentData);
        self.flux_vector_first.append(profile.doseRate);
        self.temperature_vector_first.append(profile.temperature);

        self.fill_vector_dataCutBins_first.append(profile.fill);
        self.tmpTime_vector_dataCutBins_first.append(self.time);
        self.fluence_vector_dataCutBins_first.append(self.totalDose);
        self.leakageCurrentDataCutBins_first.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.leakageCurrentDataNoRescaleCutBins_first.append(profile.leakageCurrentData);
        self.flux_vectorCutBins_first.append(profile.doseRate);
        self.temperature_vectorCutBins_first.append(profile.temperature);

        self.tmpTime_vector_data_medianVec.clear()
        self.fluence_vector_data_medianVec.clear()
        self.leakageCurrentData_medianVec.clear()
        self.leakageCurrentDataNoRescale_medianVec.clear()
        self.flux_vector_medianVec.clear()
        self.temperature_vector_medianVec.clear()

        self.tmpTime_vector_data_medianVec.append(self.time);
        self.fluence_vector_data_medianVec.append(self.totalDose);
        self.leakageCurrentData_medianVec.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.leakageCurrentDataNoRescale_medianVec.append(profile.leakageCurrentData);
        self.flux_vector_medianVec.append(profile.doseRate);
        self.temperature_vector_medianVec.append(profile.temperature);


        self.tmpTime_vector_data_medianVec.append(self.time);
        self.fluence_vector_data_medianVec.append(self.totalDose);
        self.leakageCurrentData_medianVec.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.leakageCurrentDataNoRescale_medianVec.append(profile.leakageCurrentData);
        self.flux_vector_medianVec.append(profile.doseRate);
        self.temperature_vector_medianVec.append(profile.temperature);

        self.leakage_current_hamburg_first.append(self.leakage_current_hamburg[-1])
        self.leakage_current_sheffield_first.append(self.leakage_current_sheffield[-1])
        self.leakage_current_hamburgCutBins_first.append(self.leakage_current_hamburg[-1])
        self.leakage_current_sheffieldCutBins_first.append(self.leakage_current_sheffield[-1])

        self.leakage_current_hamburg_first_scale_to_data.append(leakageCurrentScale(self.leakage_current_hamburg[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_sheffield_first_scale_to_data.append(leakageCurrentScale(self.leakage_current_sheffield[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_hamburgCutBins_first_scale_to_data.append(leakageCurrentScale(self.leakage_current_hamburg[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_sheffieldCutBins_first_scale_to_data.append(leakageCurrentScale(self.leakage_current_sheffield[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))

        self.leakage_current_hamburg_first_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_hamburg[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_sheffield_first_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_sheffield[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_hamburgCutBins_first_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_hamburg[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.leakage_current_sheffieldCutBins_first_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_sheffield[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
        self.cBunches = profile.bunches

      else:

        if(len(self.fill_vector_data_first) > len(self.fill_vector_data_second)):

          self.fill_vector_data_second.append(profile.fill);
          self.tmpTime_vector_data_second.append(self.time);
          self.fluence_vector_data_second.append(self.totalDose);
          self.leakageCurrentData_second.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
          self.leakageCurrentDataNoRescale_second.append(profile.leakageCurrentData);
          self.flux_vector_second.append(profile.doseRate);
          self.temperature_vector_second.append(profile.temperature);
          self.leakage_current_hamburg_second_scale_to_data.append(leakageCurrentScale(self.leakage_current_hamburg[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_sheffield_second_scale_to_data.append(leakageCurrentScale(self.leakage_current_sheffield[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))

          self.leakage_current_hamburg_second_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_hamburg[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_sheffield_second_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_sheffield[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))

          self.leakage_current_hamburg_second.append(self.leakage_current_hamburg[-1])
          self.leakage_current_sheffield_second.append(self.leakage_current_sheffield[-1])

          self.fill_vector_dataCutBins_second.append(profile.fill);
          self.tmpTime_vector_dataCutBins_second.append(self.time);
          self.fluence_vector_dataCutBins_second.append(self.totalDose);
          self.leakageCurrentDataCutBins_second.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
          self.leakageCurrentDataNoRescaleCutBins_second.append(profile.leakageCurrentData);
          self.flux_vectorCutBins_second.append(profile.doseRate);
          self.temperature_vectorCutBins_second.append(profile.temperature);
          self.leakage_current_hamburgCutBins_second_scale_to_data.append(leakageCurrentScale(self.leakage_current_hamburg[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_sheffieldCutBins_second_scale_to_data.append(leakageCurrentScale(self.leakage_current_sheffield[-1], self.currentTemp, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_hamburgCutBins_second_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_hamburg[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_sheffieldCutBins_second_scale_to_dataOldTemp.append(leakageCurrentScale(self.leakage_current_sheffield[-1], profile.temperature, self.userTrefK, opt.bandGap, self.boltzmanConstant))
          self.leakage_current_hamburgCutBins_second.append(self.leakage_current_hamburg[-1])
          self.leakage_current_sheffieldCutBins_second.append(self.leakage_current_sheffield[-1])



        self.tmpTime_vector_timeEnd_data_avgFill[-1] += tDelta;
        self.tmpTime_vector_data_avgFill[-1] += self.time*tDelta;
        self.fluence_vector_data_avgFill[-1] += self.totalDose;
        self.leakageCurrentData_avgFill[-1] += leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant)*tDelta;
        self.leakageCurrentDataNoRescale_avgFill[-1] += profile.leakageCurrentData*tDelta;
        self.flux_vector_avgFill[-1] += profile.doseRate*tDelta;
        self.temperature_vector_avgFill[-1] += profile.temperature*tDelta;


        self.tmpTime_vector_timeEnd_dataCutBins_avgFill[-1] += tDelta;
        self.tmpTime_vector_dataCutBins_avgFill[-1] += self.time*tDelta;
        self.fluence_vector_dataCutBins_avgFill[-1] += self.totalDose;
        self.leakageCurrentDataCutBins_avgFill[-1] += leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant)*tDelta;
        self.leakageCurrentDataNoRescaleCutBins_avgFill[-1] += profile.leakageCurrentData*tDelta;
        self.flux_vectorCutBins_avgFill[-1] += profile.doseRate*tDelta;
        self.temperature_vectorCutBins_avgFill[-1] += profile.temperature*tDelta;

        self.tmpTime_vector_data_medianVec.append(self.time);
        self.fluence_vector_data_medianVec.append(self.totalDose);
        self.leakageCurrentData_medianVec.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.leakageCurrentDataNoRescale_medianVec.append(profile.leakageCurrentData);
        self.flux_vector_medianVec.append(profile.doseRate);
        self.temperature_vector_medianVec.append(profile.temperature);



  def finalize(self):
          timeDelta = (self.tmpTime_vector_timeEnd_data_avgFill[-1])
          if(timeDelta):
            self.tmpTime_vector_data_avgFill[-1] = self.tmpTime_vector_data_avgFill[-1] / timeDelta
            self.fluence_vector_data_avgFill[-1] = self.fluence_vector_data_avgFill[-1] / timeDelta
            self.leakageCurrentData_avgFill[-1] = self.leakageCurrentData_avgFill[-1] / timeDelta
            self.flux_vector_avgFill[-1] = self.flux_vector_avgFill[-1] / timeDelta



def getProfile(filename, sensor):
    profile = [];

    myfile = open(filename, "r")

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
      bunches = 2000
      if len(words) > 5:
        bunches = float(words[6])

      profileSnapshot = DataElement();
      profileSnapshot.fill = fill;
      profileSnapshot.timestamp = timestamp;
      profileSnapshot.duration = timestep;
      profileSnapshot.bunches = bunches;
      if(timestep < 0):
        print("The timestep is negative, which should not happen. Check your inputs")
        print(line)
        profileSnapshot.duration = 0.001
      profileSnapshot.doseRate = doseRate * sensor.DoseRateScaling;
      profileSnapshot.temperature = temperature;   #Add 3 degree?????????????????
      # Scaling to milliamps from microamps
      profileSnapshot.leakageCurrentData = leakageCurrent/1.e3

      profile.append(profileSnapshot);

    myfile.close();

    return profile;
#....................................................................................................................................................    
    
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
        start_time = min(t_timevector)
        end_time   = max(t_timevector)
        gr.GetXaxis().SetRangeUser(start_time, end_time)
    gr.GetXaxis().SetTitle(xName);
    gr.GetYaxis().SetTitle(yName);

    gr.SetName(plotname);
    gr.Write();
    print("Writing ", plotname, " to file ", rootOutputFileName)
    file.Close();
    

# Read constants values from a config file 
leakageCurrentConstants= LeakageCurrentConstants()

sensor = Sensor(leakageCurrentConstants, thickness=opt.sensorThickness, width=opt.sensorWidth, length = opt.sensorLength, userTrefK = userTrefK);

# iterate through the profile and irradiate the sensor at each step
profile = getProfile(opt.input_profile, sensor);
print("Profile succesfully read. Length of the profile is: ", len(profile))
for t in range(len(profile)):
    sensor.irradiate(profile[t])
sensor.finalize()

#data output, plotting and visualisation
print("Processing finished, writing data...")

beginTime = getBeginTime(profile);
time_vector = convertDatetime(sensor.tmpTimeVector, beginTime)
time_vector_data = convertDatetime(sensor.tmpTime_vector_data, beginTime)
time_vector_data_fill_average = convertDatetime(sensor.tmpTime_vector_data_avgFill, beginTime)
time_vector_data_fill_median = convertDatetime(sensor.tmpTime_vector_data_median, beginTime)
time_vector_data_fill_first = convertDatetime(sensor.tmpTime_vector_data_first, beginTime)
time_vector_data_fill_second = convertDatetime(sensor.tmpTime_vector_data_second, beginTime)

time_vector_data_fillCutBins_average = convertDatetime(sensor.tmpTime_vector_dataCutBins_avgFill, beginTime)
time_vector_data_fillCutBins_median = convertDatetime(sensor.tmpTime_vector_dataCutBins_median, beginTime)
time_vector_data_fillCutBins_first = convertDatetime(sensor.tmpTime_vector_dataCutBins_first, beginTime)
time_vector_data_fillCutBins_second = convertDatetime(sensor.tmpTime_vector_dataCutBins_second, beginTime)

###################################################
# Hamburg model
###################################################

# All simulation

plot_vectors(time_vector, sensor.leakage_current_hamburg, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_hamburg_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_hamburg_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_scale_to_dataOldTemp", "date", opt.output_root_file)

#plots as function of dose
plot_vectors(sensor.fill_vector, sensor.leakage_current_hamburg, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Fill", "I_leak_vs_fill", "fill", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_hamburg, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_hamburg", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_sheffield, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_sheffield", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_hamburg_scale_to_data, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_hamburg_scale_to_data", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_sheffield_scale_to_data, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_sheffield_scale_to_data", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_hamburg_scale_to_dataOldTemp, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_hamburg_scale_to_dataOldTemp", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_sheffield_scale_to_dataOldTemp, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_sheffield_scale_to_dataOldTemp", "fluence", opt.output_root_file)



# First measurement in run
plot_vectors(time_vector_data_fill_first, sensor.leakage_current_hamburg_first, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_first, sensor.leakage_current_hamburg_first_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_first_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_first, sensor.leakage_current_hamburg_first_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_first_scale_to_dataOldTemp", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_hamburgCutBins_first, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_hamburgCutBins_first_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_first_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_hamburgCutBins_first_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_first_scale_to_dataOldTemp", "date", opt.output_root_file)

# Measurement at 20 minutes of data taking
plot_vectors(time_vector_data_fill_second, sensor.leakage_current_hamburg_second, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_second, sensor.leakage_current_hamburg_second_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_second_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_second, sensor.leakage_current_hamburg_second_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg_second_scale_to_dataOldTemp", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_hamburgCutBins_second, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_hamburgCutBins_second_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_second_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_hamburgCutBins_second_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburgCutBins_second_scale_to_dataOldTemp", "date", opt.output_root_file)

###################################################
# Sheffield model
###################################################
plot_vectors(time_vector, sensor.leakage_current_sheffield, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_sheffield_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_sheffield_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_scale_to_dataOldTemp", "date", opt.output_root_file)

plot_vectors(sensor.fill_vector, sensor.leakage_current_sheffield, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Fill", "I_leak_vs_fill_sheffield", "fill", opt.output_root_file)

plot_vectors(time_vector_data, sensor.leakage_current_sheffield_first, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_first", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_sheffield_first_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_first_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakage_current_sheffield_first_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_first_scale_to_dataOldTemp", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_sheffieldCutBins_first, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_sheffieldCutBins_first_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_first_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakage_current_sheffieldCutBins_first_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_first_scale_to_dataOldTemp", "date", opt.output_root_file)

plot_vectors(time_vector_data, sensor.leakage_current_sheffield_second, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_average, sensor.leakage_current_sheffield_second_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_second_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_average, sensor.leakage_current_sheffield_second_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield_second_scale_to_dataOldTemp", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_sheffieldCutBins_second, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_sheffieldCutBins_second_scale_to_data, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_second_scale_to_data", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakage_current_sheffieldCutBins_second_scale_to_dataOldTemp, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffieldCutBins_second_scale_to_dataOldTemp", "date", opt.output_root_file)

###################################################
# Data
###################################################

# Average
plot_vectors(time_vector_data_fill_average, sensor.flux_vector_avgFill, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux_fill_average", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_average, sensor.leakageCurrentData_avgFill, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_data_fill_average", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_average, sensor.leakageCurrentDataNoRescale_avgFill, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), 
                                            "Date [days]", "I_leak_per_module_data_fill_avgFill_no_rescale", "date", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_avgFill, sensor.leakageCurrentData_avgFill, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data_fill_average", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_avgFill, sensor.leakageCurrentDataNoRescale_avgFill, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data_fill_average_no_rescale", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_avgFill, sensor.temperature_vector_avgFill, "Temperature [K]", "Fill", "Temperature_vs_fill_data_fill_average", "fill", opt.output_root_file)

plot_vectors(sensor.fluence_vector_data_avgFill, sensor.leakageCurrentData_avgFill, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_per_module_data_vs_fluence_fill_average", "fluence", opt.output_root_file)

# Median
plot_vectors(time_vector_data_fill_median, sensor.flux_vector_median, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux_median", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_median, sensor.temperature_vector_median, "Temperature [K]", "Date [days]", "Temperature_data_median", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_median, sensor.leakageCurrentData_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_data_median", "date", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_median, sensor.flux_vector_median, "Flux [n_{eq}/cm^{2}/s]", "Fill", "flux_vs_fill_median", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_median, sensor.temperature_vector_median, "Temperature [K]", "Fill", "Temperature_vs_fill_data_median", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_median, sensor.leakageCurrentData_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data_median", "fill", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_median, sensor.leakageCurrentData_median, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_data_median", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_median, sensor.leakageCurrentDataCutBins_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_data_fillCutBins_median", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_median, sensor.leakageCurrentDataNoRescale_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_data_median_no_rescale", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_median, sensor.leakageCurrentDataNoRescaleCutBins_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_dataCutBins_median_no_rescale", "fill", opt.output_root_file)

plot_vectors(time_vector_data_fill_median, sensor.leakageCurrentData_median, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC),        
                                                    "Date [days]", "I_leak_vs_date_data_median", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_median, sensor.leakageCurrentDataCutBins_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Date [days]", "I_leak_vs_date_data_fillCutBins_median", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_median, sensor.leakageCurrentDataNoRescale_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Date [days]", "I_leak_vs_date_data_median_no_rescale", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_median, sensor.leakageCurrentDataNoRescaleCutBins_median, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Date [days]", "I_leak_vs_date_dataCutBins_median_no_rescale", "date", opt.output_root_file)



#First measurement in run
plot_vectors(time_vector_data_fill_first, sensor.flux_vector_first, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_first, sensor.temperature_vector_first, "Temperature [K]", "Date [days]", "Temperature_data_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_first, sensor.leakageCurrentData_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_data_first", "date", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_first, sensor.flux_vector_first, "Flux [n_{eq}/cm^{2}/s]", "Fill", "flux_vs_fill_first", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_first, sensor.temperature_vector_first, "Temperature [K]", "Fill", "Temperature_vs_fill_data_first", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_first, sensor.leakageCurrentData_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data_first", "fill", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_first, sensor.leakageCurrentData_first, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC),        
                                                    "Fill", "I_leak_vs_fill_data_first", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_first, sensor.leakageCurrentDataCutBins_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Fill", "I_leak_vs_fill_data_fillCutBins_first", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_first, sensor.leakageCurrentDataNoRescale_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_data_first_no_rescale", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_first, sensor.leakageCurrentDataNoRescaleCutBins_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_dataCutBins_first_no_rescale", "fill", opt.output_root_file)

plot_vectors(time_vector_data_fill_first, sensor.leakageCurrentData_first, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC),        
                                                    "Date [days]", "I_leak_vs_date_data_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakageCurrentDataCutBins_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Date [days]", "I_leak_vs_date_data_fillCutBins_first", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_first, sensor.leakageCurrentDataNoRescale_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Date [days]", "I_leak_vs_date_data_first_no_rescale", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_first, sensor.leakageCurrentDataNoRescaleCutBins_first, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Date [days]", "I_leak_vs_date_dataCutBins_first_no_rescale", "date", opt.output_root_file)



# Measurement at 20 minutes into run
plot_vectors(time_vector_data_fill_second, sensor.flux_vector_second, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_second, sensor.temperature_vector_second, "Temperature [K]", "Date [days]", "Temperature_data_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_second, sensor.leakageCurrentData_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_data_second", "date", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_second, sensor.flux_vector_second, "Flux [n_{eq}/cm^{2}/s]", "Fill", "flux_vs_fill_second", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_second, sensor.temperature_vector_second, "Temperature [K]", "Fill", "Temperature_vs_fill_data_second", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_second, sensor.leakageCurrentData_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data_second", "fill", opt.output_root_file)

plot_vectors(sensor.fill_vector_data_second, sensor.leakageCurrentData_second, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC),        
                                                    "Fill", "I_leak_vs_fill_data_second", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_second, sensor.leakageCurrentDataCutBins_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Fill", "I_leak_vs_fill_data_fillCutBins_second", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data_second, sensor.leakageCurrentDataNoRescale_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_data_second_no_rescale", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_dataCutBins_second, sensor.leakageCurrentDataNoRescaleCutBins_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Fill", "I_leak_vs_fill_dataCutBins_second_no_rescale", "fill", opt.output_root_file)

plot_vectors(time_vector_data_fill_second, sensor.leakageCurrentData_second, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC),        
                                                    "Date [days]", "I_leak_vs_date_data_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakageCurrentDataCutBins_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Date [days]", "I_leak_vs_date_data_fillCutBins_second", "date", opt.output_root_file)
plot_vectors(time_vector_data_fill_second, sensor.leakageCurrentDataNoRescale_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), 
                                                    "Date [days]", "I_leak_vs_date_data_second_no_rescale", "date", opt.output_root_file)
plot_vectors(time_vector_data_fillCutBins_second, sensor.leakageCurrentDataNoRescaleCutBins_second, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC),
                                                    "Date [days]", "I_leak_vs_date_dataCutBins_second_no_rescale", "date", opt.output_root_file)



###################################################
# Control plots
###################################################

plot_vectors(time_vector, sensor.temperature_vector, "Temperature [K]", "Date [days]", "temperature", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.fluence_vector, "Fluence [n_{eq}/cm^{2}]", "Date [days]", "fluence", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.flux_vector, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_data", "date", opt.output_root_file)

plot_vectors(sensor.fluence_vector, sensor.temperature_vector, "Temperature [K]", "Fluence [n_{eq}/cm^{2}", "temperature_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fill_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data, sensor.temperature_vector, "Temperature [K]", "Fill", "Temperature_vs_fill_data", "fill", opt.output_root_file)
plot_vectors(sensor.fluence_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_per_module_data_vs_fluence", "fluence", opt.output_root_file)



print(opt.output_root_file, " has been created")

