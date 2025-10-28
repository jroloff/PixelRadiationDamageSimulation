import ROOT 
from optparse import OptionParser
import datetime
import math
from array import array


parser = OptionParser()
#These are the main options that should be changed
parser.add_option("--input_profile", default="/Users/jroloff/Work/trackerMonitoring/profile_BPix_BmI_SEC1_LYR1_phase1_newL1.txt", help="Input profile file name, should have been made using PixelMonitoring repository")

parser.add_option("--output_root_file", default="testFile.root", help="Output ROOT file name")
# These are options that shoudlb e changed if you are using different sensors etc.
parser.add_option("--timestep", default=1, help="step size is 1 second, do not change this!")
parser.add_option("--userTrefC", type="float",  default=0., help="reference temprature in celsius")
parser.add_option("--bandGap", default=1.21, help="eV used for scaling temperatures")
# Dimensions of sensors for CMS are taken from page 68 of https://cds.cern.ch/record/2773267
parser.add_option("--sensorThickness", default=0.0285, help="Thickness of the sensor in cm. This is correct for CMS")
parser.add_option("--sensorWidth", default=6.48, help="Width of the sensor in cm. 285 is correct for CMS")
parser.add_option("--sensorLength", default=1.62, help="Length of the sensor in cm. 285 is correct for CMS")

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
    return leakageCurrent*(T*T/(Tref*Tref)*math.exp(-(bandGap/(2*boltzmanConstant))*(1/T - 1/Tref)));

# An array of DataElements that contain the conditions at each time step


#Information about the conditions for each data point
class DataElement:
  def __init__(self):
    fill = 0
    timestamp = 0
    duration = 0
    temperature = 0.
    doseRate = 0
    leakageCurrentData = 0.

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

 
        self.G_i = [];
        self.leakage_current_hamburg = [];
        self.leakage_current_sheffield = [];
        self.leakage_current_per_module = [];
        self.leakage_current_per_volume = [];
        self.leakage_current_per_module_sheffield = [];
        self.leakage_current_per_volume_sheffield = [];
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

        self.temperature_vector = [];

        self.tmpTime_vector_data = [];
        self.fluence_vector_data = [];
        self.fill_vector_data = [];
        self.leakageCurrentData = [];




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

    # Only fill these ones if there is any data, not just for every time step
    if (profile.leakageCurrentData > 0. and profile.doseRate > 0.):
        self.fill_vector_data.append(profile.fill);
        self.tmpTime_vector_data.append(self.time);
        self.fluence_vector_data.append(self.totalDose);
        self.leakageCurrentData.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap, self.boltzmanConstant));
        self.flux_vector.append(profile.doseRate);


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

            if(self.doseRate_vector[i] > 0):
              test1 = math.exp(-self.sheffield_alpha0_contribution_vector[-1])
              #print(self.leakageCurrentConstants.alpha_Tref, self.leakageCurrentConstants.Ak[0], 1/self.sheffield_alpha0_contribution_vector[-1] , (1-test1), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[-1] - self.sheffield_alpha0_contribution_sum_vector[i-1])), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[-1])), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[i-1] )), sheffield_G_0_tmp, sheffield_G_0_tmp + sheffield_G_1_tmp + sheffield_G_2_tmp + sheffield_G_3_tmp + sheffield_G_4_tmp)
              #print(self.leakageCurrentConstants.alpha_1,  math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])), self.leakageCurrentConstants.alpha_0_star, self.leakageCurrentConstants.beta , math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1]), G_i_tmp)

              #print(self.leakageCurrentConstants.alpha_Tref, self.leakageCurrentConstants.Ak[0], 1/self.sheffield_alpha0_contribution_vector[-1] , (1-test1), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[-1] - self.sheffield_alpha0_contribution_sum_vector[i-1])), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[-1])), math.exp(-(self.sheffield_alpha0_contribution_sum_vector[i-1] )), sheffield_G_0_tmp, sheffield_G_0_tmp + sheffield_G_1_tmp + sheffield_G_2_tmp + sheffield_G_3_tmp + sheffield_G_4_tmp)
              #print(math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])), (self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1])), math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1]), self.leakageCurrentConstants.alpha_0_star,  self.leakageCurrentConstants.beta,  self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])))
    

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


    # TODO: may need to fix these functions, not sure
    # Convert to per volume or per module units
    self.leakage_current_per_volume.append(self.getPerVolume(self.leakage_current_hamburg[-1], self.userTrefK, self.Tref_leakage))
    self.leakage_current_per_module.append(self.getPerModule(self.leakage_current_hamburg[-1], self.userTrefK, self.Tref_leakage))
    self.leakage_current_per_volume_sheffield.append(self.getPerVolume(self.leakage_current_sheffield[-1], self.userTrefK, self.Tref_sheff))
    self.leakage_current_per_module_sheffield.append(self.getPerModule(self.leakage_current_sheffield[-1], self.userTrefK, self.Tref_sheff))





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
        #print("X-axis range:", start_time, end_time)
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

#data output, plotting and visualisation
print("Processing finished, writing data...")

beginTime = getBeginTime(profile);
time_vector = convertDatetime(sensor.tmpTimeVector, beginTime)
time_vector_data = convertDatetime(sensor.tmpTime_vector_data, beginTime)

# plots as function of time
plot_vectors(time_vector, sensor.alpha_vec, "#alpha [A/cm]", "Date [days]", "alpha", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_hamburg, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_hamburg", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_module, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_hamburg", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_volume, "I_{leak} (@%d C) [mA/cm^{3}]", "Date [days]", "I_leak_volume_hamburg", "date", opt.output_root_file)

plot_vectors(time_vector, sensor.leakage_current_sheffield, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak_sheffield", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_module_sheffield, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_sheffield", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_volume_sheffield, "I_{leak} (@%d C) [mA/cm^{3}]", "Date [days]", "I_leak_volume_sheffield", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.G_i, "G_{i} [A/cm^{3}]", "Date [days]", "G_i", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.temperature_vector, "Temperature [K]", "Date [days]", "temperature", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.fluence_vector, "Fluence [n_{eq}/cm^{2}]", "Date [days]", "fluence", "date", opt.output_root_file)

plot_vectors(time_vector_data, sensor.flux_vector, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_data", "date", opt.output_root_file)
plot_vectors(sensor.fill_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA], 1 ROG"%(opt.userTrefC), "Fill", "I_leak_vs_fill_data", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector_data, sensor.temperature_vector, "Temperature [K]", "Fill", "Temperature_vs_fill_data", "fill", opt.output_root_file)

#plots as function of dose
plot_vectors(sensor.fluence_vector, sensor.alpha_vec, "#alpha [A/cm]", "Fluence [n_{eq}/cm^{2}]", "alpha_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_hamburg, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_hamburg", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_module, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Date [days]", "I_leak_per_module_vs_fluence_hamburg", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_volume, "I_{leak} (@%d C) [mA/cm^{3}]"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_volume_vs_fluence_hamburg", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_sheffield, "I_{leak} [mA/cm^{2}]", "Fluence [n_{eq}/cm^{2}]", "I_leak_vs_fluence_sheffield", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_module_sheffield, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Date [days]", "I_leak_per_module_vs_fluence_sheffield", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.leakage_current_per_volume_sheffield, "I_{leak} (@%d C) [mA/cm^{3}]"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_volume_vs_fluence_sheffield", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.G_i, "G_{i} [A/cm^{3}]", "Fluence [n_{eq}/cm^{2}]", "G_i_vs_fluence", "fluence", opt.output_root_file)
plot_vectors(sensor.fluence_vector, sensor.temperature_vector, "Temperature [K]", "Fluence [n_{eq}/cm^{2}", "temperature_vs_fluence", "fluence", opt.output_root_file)

plot_vectors(sensor.fluence_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA] per module"%(opt.userTrefC), "Fluence [n_{eq}/cm^{2}]", "I_leak_per_module_data_vs_fluence", "fluence", opt.output_root_file)

plot_vectors(sensor.fill_vector, sensor.leakage_current_hamburg, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Fill", "I_leak_vs_fill", "fill", opt.output_root_file)
plot_vectors(sensor.fill_vector, sensor.leakage_current_sheffield, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Fill", "I_leak_vs_fill_sheffield", "fill", opt.output_root_file)


print(opt.output_root_file, " has been created")

