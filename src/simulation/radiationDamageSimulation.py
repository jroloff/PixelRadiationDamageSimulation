import ROOT 
from optparse import OptionParser
import datetime
import math
from array import array


parser = OptionParser()
#These are the main options that should be changed
parser.add_option("--input_profile", default="/afs/cern.ch/user/s/singhsh/PixelMonitoring/data/radiation_simulation/profiles/per_year/BPix_BmI_SEC1_LYR1/profile_BPix_BmI_SEC1_LYR1_2022.txt", help="Input profile file name, should have been made using PixelMonitoring repository")
parser.add_option("--inputAnnealingConstants", default="config/annealing_constants.py", help="Input annealing constants file name")
parser.add_option("--output_root_file", default="testFile.root", help="Output ROOT file name")
# These are options that shoudlb e changed if you are using different sensors etc.
parser.add_option("--timestep", default=1, help="step size is 1 second, do not change this!")
parser.add_option("--donorremovalfraction", default=0.99, help="fraction of donors which can be removed by donor removal")
parser.add_option("--userTrefC", type="float",  default=0., help="reference temprature in celsius")
parser.add_option("--bandGap", default=1.21, help="eV used for scaling temperatures")
# Dimensions of sensors for CMS are taken from page 68 of https://cds.cern.ch/record/2773267
parser.add_option("--sensorThickness", default=0.0285, help="Thickness of the sensor in cm. This is correct for CMS")
parser.add_option("--sensorWidth", default=6.48, help="Width of the sensor in cm. 285 is correct for CMS")
parser.add_option("--sensorLength", default=1.62, help="Length of the sensor in cm. 285 is correct for CMS")
parser.add_option("--Ndonor_0", default=1.7e12, help="Initial number of donors for the sensor")


(opt, args) = parser.parse_args()


# Defining constants used in the analysis
boltzmanConstant = 8.6173303e-5
absoluteZero = 273.15

# set a reference temperature for the volume corrected leakage current plot (only affects a couple plots)
userTrefK = absoluteZero + opt.userTrefC; 

#Recreating this file so we don't end up with a bunch of graphs from old runs
ROOT.gROOT.SetBatch(ROOT.kTRUE)
file = ROOT.TFile(opt.output_root_file, "RECREATE");
file.Close()


#Information about the conditions for each data point
class DataElement:
  def __init__(self):
    fill = 0
    timestamp = 0
    duration = 0
    temperature = 0.
    doseRate = 0
    coolPipeTemp = 0.
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



#Class to store a set of annealing specific constants and compute temperature dependent values
class AnnealingConstants:
  def __init__(self, filename):
    self.gA = 0.0
    self.gY = 0.0
    self.gC = 0.0
    self.ka = 0.0
    self.ky = 0.0
    self.Ea = 0.0
    self.Ey = 0.0
    self.cc = 0.0

    print("Reading annealing constant file: ", filename)
    myfile = open(filename, "r")

    quantity = 0;
    valueStr = "";
    value = 0.;

    for line in myfile:
        found =  line.find("#");
        if (found != -1) :
          line = line[:found];
        found =  line.find(":");
        if (found ==-1): 
          continue;
        quantity = line[: found].strip();
        valueStr = line[found+1: len(line)];
        found =  valueStr.find(",");
        if (found != -1) :
            valueStr = valueStr[:found].strip();

        value = float(valueStr);

        if (quantity == "\"gA\"" or quantity == "'gA'"):
          self.gA = value;
        if (quantity == "\"gY\"" or quantity == "'gY'"):
          self.gY = value
        if (quantity == "\"gC\"" or quantity == "'gC'"):
          self.gC = value
        if (quantity == "\"ka\"" or quantity == "'ka'"):
          self.ka = value
        if (quantity == "\"ky\"" or quantity == "'ky'"):
          self.ky = value
        if (quantity == "\"Ea\"" or quantity == "'Ea'"):
          self.Ea = value
        if (quantity == "\"Ey\"" or quantity == "'Ey'"):
          self.Ey = value
        if (quantity == "\"cc\"" or quantity == "'cc'"):
          self.cc = value

    if not self.gA or not self.gY or not self.gC or not self.ka or not self.ky or not self.Ea or not self.Ey or not self.cc:
      print("WARNING: Problem with initializing annealing constants")

    myfile.close();




# all atributes of a sensor (like irradiation, doping concentration, leakage currents, etc) will be stored in this class

class Sensor:
  # Small value -- just something to compare to to make sure numbers are reasonable
  # Epsilon is used for something else...
  smallNumber = 1e-30;

  # Fundamental constants that are relevant for some of these calculations
  boltzmanConstant = 8.6173303e-5
  # Electron charge = 1.6 x 10^-19 C???
  # Not sure the reason for the discrepancy of 10e6
  electronCharge = 1.6021766208e-13

  #Permittivity of vacuum
  #epsilon0 = 8.85E-12 C^2 kg^-1 m^-3 s^2
  epsilon0 = 8.854187817E-12
  #Permittivity of silicon
  epsilon = 11.68


  #T-ref cannot be changed in a trivial way since the other parameters (e.g. alpha*_0) were evaluated at this reference temperature!!!
  # We have 2 different values of Tref based on the difference in the reference temperature for  tables 9 and 10
  # https://cds.cern.ch/record/2773267/
  # These are very close to each other, but this difference is intentional
  # This Tref is useful for scaling of eq. 24 of https://cds.cern.ch/record/2773267/
  # TODO: Actually, isn't this for the Sheffield model???
  Tref = 293.15;
  # Room temperature (useful for scaling of the tau and Theta functions in eq. 22 of https://cds.cern.ch/record/2773267/)
  Tref_leakage = 294.15

  # According to page 5 of https://iopscience.iop.org/article/10.1088/1748-0221/14/06/P06012/pdf, t0 is one minute = 60 seconds
  t0 = 60.

  # This parameter is multiplied to all fluence rates from the input profile
  # TODO: Not sure where these numbers come from
  # TODO: make this configurable
  # Could this be the scaling factor that is discussed in 5.3.1.3?
  # https://cds.cern.ch/record/2773267/files/10.23731_CYRM-2021-001.59.pdf
  # This number is very close to 1, and it mentioned that this was for bpix
  DoseRateScaling = 0.0859955711871/8.9275E-02;

  def __init__(self, leakageCurrentConstants, annealingConstants, thickness, width, length, Ndonor_0, userTrefK):
        # Things related to conditions
        self.time = 0.0;
        self.totalDose = 0.;
        # IBL       // initial donor concentration (right now the code only works for originally n-type bulk material!)
        # TODO:Figure out how to get the right number
        self.Ndonor_0 = Ndonor_0
        self.userTrefK = userTrefK

        self.annealingConstants = annealingConstants
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

        # Initialize without any defects
        self.Nacceptor = 0.0;
        self.Nacceptor_stable_constdamage = 0.0
        self.Nacceptor_stable_reverseannealing = 0.0;
        self.Nacceptor_reversible=0;
        self.Ndonor = Ndonor_0;
        self.Neffective = 0.0

        self.Nbenef_anneal_g1 = 0.0;
        self.Nrevers_anneal_g1 = 0.0;
        self.Nnadefects_g1 = 0.0;
        self.Nconstdamage_g1 = 0.0;
        self.Nneutral_reversible=0;
        self.Nnadefects_g1=0;
 
        self.G_i = [];
        self.leakage_current = [];
        self.leakage_current_per_module = [];
        self.leakage_current_per_volume = [];
        self.alpha_vec = [];
        self.powerconsumption = [];

        self.N_benef_anneal_g1_vec = [];
        self.N_revers_anneal_g1_vec = [];
        self.N_nadefects_g1_vec = [];
        self.N_constdamage_g1_vec = [];
        self.N_donor_vec = [];


        self.Theta_vector = []
        self.alpha1_contribution_vector = []
        self.alpha1_contribution_sum_vector = []
        self.beta_contribution_vector = []
        self.beta_contribution_sum_vector = []

        self.V_dep_vector = [];
        self.Neff_vector = [];
        self.Nacceptor_vector = [];
        self.temperature_vector = [];

        self.tmpTime_vector_data = [];
        self.fluence_vector_data = [];
        self.fill_vector_data = [];
        self.leakageCurrentData = [];

        self.tmpTime_vector_timeStart_data_avgFill = [];
        self.tmpTime_vector_timeEnd_data_avgFill = [];
        self.tmpTime_vector_data_avgFill = [];
        self.fluence_vector_data_avgFill = [];
        self.fill_vector_data_avgFill = [];
        self.leakageCurrentData_avgFill = [];
        self.flux_vector_avgFill = [];


        self.Ndonor_stable_donorremoval= opt.donorremovalfraction*abs(self.Nacceptor-self.Ndonor);
        self.Ndonor_const=self.Ndonor-self.Ndonor_stable_donorremoval;

        #Definition of Hamburg model depletion voltage functions
        #Allowing large range of values to determine lots of scaling 
        # Most of the parameters need to be evaluated with considerations for the irradiation, so they get set later (and at each time step)
        self.Naccept_rev_TF1 = ROOT.TF1("Naccept_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
        self.Naccept_rev_TF1.SetParameter(0, self.annealingConstants.gA);

        self.Naccept_const_TF1 = ROOT.TF1("Naccept_const_TF1","[0] * [1]*x",0,10000000);
        self.Naccept_const_TF1.SetParameter(0, self.annealingConstants.gC);

        self.Nneutrals_rev_TF1 = ROOT.TF1("Ndonor_rev_TF1","[0]*[1]*(1-TMath::Exp(-[2]*x))/[2] + [3]*TMath::Exp(-[2]*x)",0,10000000);
        self.Nneutrals_rev_TF1.SetParameter(0, self.annealingConstants.gY);

        self.Ndonor_const_TF1 = ROOT.TF1("Ndonor_const_TF1","-[0]*(1-TMath::Exp(-[1]*[2]*x))",0,10000000);
        self.Ndonor_const_TF1.SetParameter(1, self.annealingConstants.cc);

        self.Naccept_rev_TF1_approx = ROOT.TF1("Naccept_rev1_TF1","[0]*[1]*x+[3]*TMath::Exp(-[2]*x)",0,10000000);
        self.Naccept_rev_TF1_approx.SetParameter(0, self.annealingConstants.gA);

        self.Ndonor_neutrals_rev_TF1_approx = ROOT.TF1("Ndonor_rev1_TF1","[0]*[1]*x+[3]*TMath::Exp(-[2]*x)",0,10000000);
        self.Ndonor_neutrals_rev_TF1_approx.SetParameter(0, self.annealingConstants.gY);

        self.Ndonor_TF1 = ROOT.TF1("Ndonor_TF1","[0]*[1]/[2] * ([2]*x+TMath::Exp(-[2]*x)-1) +[3]*(1-TMath::Exp(-[2]*x)) ",0,10000000);
        self.Ndonor_TF1.SetParameter(0, self.annealingConstants.gY);

        self.Ndonor_g1_TF1 = ROOT.TF1("Ndonor_TF1","[0]/[1] * ([1]*x+TMath::Exp(-[1]*x)-1) +[2]*(1-TMath::Exp(-[1]*x)) ",0,10000000);


  def get_ka(self, T):
    return self.annealingConstants.ka*math.exp(-self.annealingConstants.Ea/(self.boltzmanConstant*T));

  def get_ky(self, T):
    return self.annealingConstants.ky*math.exp(-self.annealingConstants.Ey/(self.boltzmanConstant*T));

  # TODO: maybe this should be at the current temperature, since we will need to scale the data, not just the predictions????
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


  #calculate depletion voltage from given effective doping concentration
  def NtoV_function(self):
    # q/2*epsilon*epsilon_0 * Neff*D^2
    #print(self.Neffective, self.electronCharge/(2*self.epsilon*self.epsilon0), self.electronCharge, self.epsilon, self.epsilon0, self.depth*self.depth)
    # Note: converting depth from cm to m (1e-2 x 1e-2 = 1e-4) to get units correct
    value = (self.electronCharge/(2*self.epsilon*self.epsilon0))*abs(self.Neffective)*(self.depth*1e-2 * self.depth*1e-2) 
    return value;

    #return 1.6021766208e-13/(2*11.68*8.854187817)*fabs(doping_conc)*thickness*thickness;                // q/2*epsilon*epsilon_0 * Neff*D^2
    # 



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
        #if( profile.temperature - userTrefK < 0.5 and userTrefK- profile.temperature< 0.5):
        print(profile.fill, self.totalDose, userTrefK, profile.temperature, profile.leakageCurrentData, leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap))
        self.fill_vector_data.append(profile.fill);
        self.tmpTime_vector_data.append(self.time);
        self.fluence_vector_data.append(self.totalDose);
        self.leakageCurrentData.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap));
        self.flux_vector.append(profile.doseRate);

    # Only fill these ones if there is any data, not just for every time step
    if (profile.leakageCurrentData > 0. and profile.doseRate > 0.):
      if(len(self.fill_vector_data_avgFill)==0 or self.fill_vector_data[-1] != profile.fill):
        if(len(self.fill_vector_data_avgFill)!=0):
          timeDelta = (self.tmpTime_vector_timeEnd_data_avgFill[-1] - self.tmpTime_vector_timeStart_data_avgFill[-1])
          if(timeDelta):
            self.tmpTime_vector_data_avgFill[-1] = self.tmpTime_vector_data_avgFill[-1] / timeDelta
            self.fluence_vector_data_avgFill[-1] = self.fluence_vector_data_avgFill[-1] / timeDelta
            self.leakageCurrentData_avgFill[-1] = self.leakageCurrentData_avgFill[-1] / timeDelta
            self.flux_vector_avgFill[-1] = self.flux_vector_avgFill[-1] / timeDelta

        self.tmpTime_vector_timeStart_data_avgFill.append(self.time);
        self.tmpTime_vector_timeEnd_data_avgFill.append(self.time);
        self.fill_vector_data_avgFill.append(profile.fill);
        self.tmpTime_vector_data_avgFill.append(self.time);
        self.fluence_vector_data_avgFill.append(self.totalDose);
        self.leakageCurrentData_avgFill.append(leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap));
        self.flux_vector_avgFill.append(profile.doseRate);


      else:
        self.tmpTime_vector_timeEnd_data_avgFill[-1] = self.time;
        self.tmpTime_vector_data_avgFill[-1] += self.time;
        self.fluence_vector_data_avgFill[-1] += self.totalDose;
        self.leakageCurrentData_avgFill[-1] += leakageCurrentScale(profile.leakageCurrentData, userTrefK, profile.temperature, opt.bandGap);
        self.flux_vector_avgFill[-1] += profile.doseRate;

    
    ##################################################################
    # Set all the parameters based on current conditions
    ##################################################################
    self.Naccept_rev_TF1.SetParameter(1, profile.doseRate);
    self.Naccept_rev_TF1.SetParameter(2, self.get_ka(profile.temperature));
    self.Naccept_rev_TF1.SetParameter(3, self.Nacceptor_reversible);
  
    self.Naccept_rev_TF1_approx.SetParameter(1, profile.doseRate);
    self.Naccept_rev_TF1_approx.SetParameter(2, self.get_ka(profile.temperature));
    self.Naccept_rev_TF1_approx.SetParameter(3, self.Nacceptor_reversible);
  
    self.Naccept_const_TF1.SetParameter(1, profile.doseRate);
  
    self.Nneutrals_rev_TF1.SetParameter(1, profile.doseRate);
    self.Nneutrals_rev_TF1.SetParameter(2, self.get_ky(profile.temperature));
    self.Nneutrals_rev_TF1.SetParameter(3, self.Nneutral_reversible);
  
    self.Ndonor_neutrals_rev_TF1_approx.SetParameter(1, profile.doseRate);
    self.Ndonor_neutrals_rev_TF1_approx.SetParameter(2, self.get_ky(profile.temperature));
    self.Ndonor_neutrals_rev_TF1_approx.SetParameter(3, self.Nneutral_reversible);
  
    self.Ndonor_const_TF1.SetParameter(0, self.Ndonor_stable_donorremoval);
    self.Ndonor_const_TF1.SetParameter(2, profile.doseRate);
  
    self.Ndonor_TF1.SetParameter(1, profile.doseRate);
    self.Ndonor_TF1.SetParameter(2, self.get_ky(profile.temperature));
    self.Ndonor_TF1.SetParameter(3, self.Nneutral_reversible);
  
    self.Ndonor_g1_TF1.SetParameter(0, profile.doseRate);
    self.Ndonor_g1_TF1.SetParameter(1, self.get_ky(profile.temperature));
    self.Ndonor_g1_TF1.SetParameter(2, self.Nnadefects_g1);
  
    # Determine the impact of the current parameters on the acceptors and donors
    if(self.get_ka(profile.temperature)>self.smallNumber and self.get_ky(profile.temperature)>self.smallNumber):
      self.Nacceptor_reversible             =  self.Naccept_rev_TF1.Eval(profile.duration);
      self.Nneutral_reversible              =  self.Nneutrals_rev_TF1.Eval(profile.duration);
    else:
      print("Potential numerical problem due to ultra low temp and therby caused very small ky and ka values. Using an approximation in order to perform calculation. In general, no problem!")
      self.Nacceptor_reversible             =  self.Naccept_rev_TF1_approx.Eval(profile.duration);
      self.Nneutral_reversible              =  self.Ndonor_neutrals_rev_TF1_approx.Eval(profile.duration);

    ##################################################################
    # Evaluate the functions to determine the annealing at this time step
    # Append relevant values to vectors
    ##################################################################
    self.Ndonor_stable_donorremoval        += self.Ndonor_const_TF1.Eval(profile.duration);
    self.Nacceptor_stable_constdamage      +=  self.Naccept_const_TF1.Eval(profile.duration);
    self.Nacceptor_stable_reverseannealing +=  self.Ndonor_TF1.Eval(profile.duration);

    self.Nrevers_anneal_g1 += self.Ndonor_g1_TF1.Eval(profile.duration);
    self.N_revers_anneal_g1_vec.append(self.Nrevers_anneal_g1);

    self.Nconstdamage_g1  = self.Nacceptor_stable_constdamage/self.annealingConstants.gC;
    self.N_constdamage_g1_vec.append(self.Nconstdamage_g1);

    self.Nbenef_anneal_g1 = self.Nacceptor_reversible/self.annealingConstants.gA;
    self.N_benef_anneal_g1_vec.append(self.Nbenef_anneal_g1);

    self.Nnadefects_g1     = self.Nneutral_reversible/self.annealingConstants.gY;
    self.N_nadefects_g1_vec.append(self.Nnadefects_g1);
  
    self.Nacceptor =  self.Nacceptor_reversible + self.Nacceptor_stable_constdamage + self.Nacceptor_stable_reverseannealing;
    self.Nacceptor_vector.append(self.Nacceptor);

    self.Ndonor    =  self.Ndonor_const         + self.Ndonor_stable_donorremoval;
    self.N_donor_vec.append(self.Ndonor);

    self.Neffective = self.Nacceptor-self.Ndonor
    self.Neff_vector.append(self.Neffective);

    # Append all the acceptor and donor information for this current time step to the vectors
    self.V_dep_vector.append(self.NtoV_function());

    # Equations taken from page 5 of https://cds.cern.ch/record/2773267
    tauInverse = self.leakageCurrentConstants.k01 * math.exp(-self.leakageCurrentConstants.E1 / (self.boltzmanConstant * int(profile.temperature)))
    # The individual term in the sum over j of t_j / tau(T_j)
    self.alpha1_contribution_vector.append(profile.duration *tauInverse)

    # Equations taken from page 5 of https://cds.cern.ch/record/2773267
    self.Theta_vector.append(math.exp(-self.leakageCurrentConstants.E1_star/self.boltzmanConstant * (1.0 / float(int(profile.temperature)) - 1.0/self.Tref) ))

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

    G_i_tmp = 0;  
    # TODO: not sure if there is a better name for G_i...
    for i in range(len(self.beta_contribution_vector)): 
        if i>0 :
            # Sum over everything minus sum up to i to give sum from i to n
            G_i_tmp += int(self.doseRate_vector[i]) * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1] - self.alpha1_contribution_sum_vector[i-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1] - self.beta_contribution_sum_vector[i-1]));
        else:
            # If we are on the first time increment, then the sum from j=i to i-1 is just the sum from 0 to i-1, so we don't need to do anything special
            G_i_tmp += int(self.doseRate_vector[i]) * self.duration_vector[i]*( self.leakageCurrentConstants.alpha_1 * math.exp(-(self.alpha1_contribution_sum_vector[-1])) + self.leakageCurrentConstants.alpha_0_star  - self.leakageCurrentConstants.beta * math.log(self.beta_contribution_sum_vector[-1]))
          

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
    self.leakage_current.append(G_i_tmp * amps_to_mA * self.depth);
    self.powerconsumption.append(self.leakage_current[-1] * self.NtoV_function());

    # TODO: may need to fix these functions, not sure
    # Convert to per volume or per module units
    self.leakage_current_per_volume.append(self.getPerVolume(self.leakage_current[-1]))
    self.leakage_current_per_module.append(self.getPerModule(self.leakage_current[-1]))



#function to read in the temp/rad profile
def leakageCurrentScale(leakageCurrent, T, Tref, bandGap):
    # Taken from equation 25 of https://cds.cern.ch/record/2773267
    return leakageCurrent*(T*T/(Tref*Tref)*math.exp(-(bandGap/(2*boltzmanConstant))*(1/T - 1/Tref)));

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
annealingConstants = AnnealingConstants(opt.inputAnnealingConstants);
leakageCurrentConstants= LeakageCurrentConstants()

sensor = Sensor(leakageCurrentConstants, annealingConstants, thickness=opt.sensorThickness, width=opt.sensorWidth, length = opt.sensorLength, Ndonor_0 = opt.Ndonor_0, userTrefK = userTrefK);

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
plot_vectors(time_vector, sensor.Neff_vector, "N_{eff} [1/cm^{3}]", "Date [days]", "Neff", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.N_donor_vec, "N_{donor} [1/cm^{3}]", "Date [days]", "Ndonors", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.Nacceptor_vector, "N_{acceptor} [1/cm^{3}]", "Date [days]", "Nacceptors", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.V_dep_vector, "U_{dep} [V]", "Date [days]", "U_dep", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.alpha_vec, "#alpha [A/cm]", "Date [days]", "alpha", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current, "I_{leak} (@%d C) [mA/cm^{2}]"%(opt.userTrefC), "Date [days]", "I_leak", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_module, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.leakage_current_per_volume, "I_{leak} (@%d C) [mA/cm^{3}]", "Date [days]", "I_leak_volume", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.G_i, "G_{i} [A/cm^{3}]", "Date [days]", "G_i", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.powerconsumption, "P [mW/cm^{2}]", "Date [days]", "powerconsumption", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.temperature_vector, "Temperature [K]", "Date [days]", "temperature", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.fluence_vector, "Fluence [n_{eq}/cm^{2}]", "Date [days]", "fluence", "date", opt.output_root_file)

plot_vectors(time_vector_data, sensor.flux_vector, "Flux [n_{eq}/cm^{2}/s]", "Date [days]", "flux", "date", opt.output_root_file)
plot_vectors(time_vector_data, sensor.leakageCurrentData, "I_{leak} (@%d C) [mA],  1 ROG"%(opt.userTrefC), "Date [days]", "I_leak_per_module_data", "date", opt.output_root_file)

plot_vectors(time_vector, sensor.N_benef_anneal_g1_vec, "N_{dep, benef_anneal_g1} [V]", "Date [days]", "N_benef_anneal_g1", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.N_revers_anneal_g1_vec, "N_{dep, revers_anneal_g1} [V]", "Date [days]", "N_revers_anneal_g1", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.N_nadefects_g1_vec, "N_{dep, nadefects_g1} [V]", "Date [days]", "N_nadefects_g1", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.N_constdamage_g1_vec, "N_{dep, constdamage_g1} [V]", "Date [days]", "N_constdamage_g1", "date", opt.output_root_file)
plot_vectors(time_vector, sensor.N_donor_vec, "N_{dep, donor} [V]", "Date [days]", "N_donor", "date", opt.output_root_file)

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

