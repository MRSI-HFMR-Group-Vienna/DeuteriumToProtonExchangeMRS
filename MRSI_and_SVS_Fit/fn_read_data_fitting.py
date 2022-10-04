#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


############# Author: Fabian Niess
############# Fitting 1H MRSI data after given Deuterium (Bednarik et al. Data)


import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy import stats

def exp_decrease(t,m0,tau,const):
    m = m0 * np.exp(-t/tau) + const
    return m
    
def exp_increase(t,m0,tau,d):
    m = m0 * (1-d*np.exp(-t/tau))
    return m    

def fit(function,time,signals,param):
    	
    	try:
    		popt, pcov = curve_fit(function, time, signals, param)
    		perr = np.sqrt(np.diag(pcov))
    		return [popt, perr]
    	except:
    		return [[0,0,0],[0,0,0]]

##### Data Points SVS for 5 Subjects
glu4_deu=[]

#### old fit
#glu4_deu.append(np.array([7.591597578,7.061428165,6.651216276,6.837578857,6.552679738,6.529116653,6.204588709]))
#glu4_deu.append(np.array([7.28942269,6.971726386,6.8020124,6.312532426,6.532953641,6.459479903]))
#glu4_deu.append(np.array([7.499749229,7.135764737,6.880245538,6.993925508,6.725891083,6.537119756,6.571536628,6.489144723]))
#glu4_deu.append(np.array([7.40554402,6.834823994,6.599946665,6.650960746,6.254538829]))
#glu4_deu.append(np.array([7.361203632,7.151363164,6.722239406,6.328788527,6.34767417,6.398035882,5.997240587]))

#### new fit

glu4_deu.append(np.array([7.4705690053061,7.06785446107741,6.60301905630281,6.82151311753327,6.59016646446572,6.50234042024564,6.12318896105161]))
glu4_deu.append(np.array([7.20456569716295,6.9417159862926,6.80615176533191,6.31770663332086,6.44913148875604,6.39842426106845]))
glu4_deu.append(np.array([7.42465787271561,7.0836179620009,6.83957105341606,6.92300589395789,6.67687311435951,6.47245775503204,6.5006170137149,6.46724307749817]))
glu4_deu.append(np.array([7.32264613933337,6.79018667404948,6.53724185820603,6.55424655170811,6.17908050131847]))
glu4_deu.append(np.array([7.24369296966067,7.06742697619268,6.61207315973371,6.27108239856052,6.27842681495502,6.36655981168901,5.92169801865076,6.04130708564689]))

glu4_dex=[]
#### old fit
#glu4_dex.append(np.array([8.087082454,7.904463642,7.858245301,7.800754193,8.004791261,7.95180926,7.899954536]))
#glu4_dex.append(np.array([7.470822804,7.720892186,7.559389043,7.711514585,7.678172]))
#glu4_dex.append(np.array([7.997229466,8.042075693,8.155755663,8.071277887,8.10778063,8.238147568]))
#glu4_dex.append(np.array([7.669753575,7.559502167,7.718638633,7.450290867,7.792486273,7.680154651]))
#glu4_dex.append(np.array([7.704162352,7.778735703,7.672652485,7.9131778,7.864862671,7.838604449,7.831252147]))

#### new fit

glu4_dex.append(np.array([7.92137279157849,7.71620844718297,7.67337193571578,7.68915380836158,7.90671819555024,7.85260891790746,7.79399053379446]))
glu4_dex.append(np.array([7.43227044067572,7.65837484073553,7.46040324621312,7.63961963704393,7.58960576053301]))
glu4_dex.append(np.array([7.86894839860082,7.89815059279046,8.01495936954902,7.92943865799365,7.95238323914265,8.16409914701753]))
glu4_dex.append(np.array([7.57614388895612,7.43676946815435,7.6146278708193,7.34523999777706,7.69783648025319,7.57718399657404]))
glu4_dex.append(np.array([7.53400907299776,7.58862617487925,7.49934821988067,7.80919524016985,7.70941399634791,7.75037682275902,7.63273998793737]))


glu4_extra=[]

glu4_extra.append(np.array([7.14173290249229,6.84003024517518,6.36592606939114,6.38855376868992,6.07068846901653,6.05668084564109]))
glu4_extra.append(np.array([6.88333706136686,6.66199459237625,6.52236266100367,6.51512248678435,6.17173136666808,6.21000085897019]))
glu4_extra.append(np.array([7.17268679452285,6.71059237731585,6.80771399155183,6.37424489222491,6.2566766223603,6.17386766706435]))
glu4_extra.append(np.array([7.46569512471825,7.04779529765991,6.76635255698797,6.64695260639987,6.56912942432013,6.38896342745059]))
glu4_extra.append(np.array([7.21759808381999,6.72492217775279,6.72180397581565,6.38399876595945,6.30812251882251,6.27278289686833]))


gln4_deu=[]
#### old fit
#gln4_deu.append(np.array([3.3855869,3.223858452,3.149956049,3.166021789,2.962522418,3.170305986,2.764378294]))
#gln4_deu.append(np.array([2.392553274,2.246640639,2.110041577,2.279755563,2.145226184,2.147295866]))
#gln4_deu.append(np.array([2.921262354,2.868072644,2.87120145,2.828441094,2.672000768,2.687644801,2.588565928,2.590651799]))
#gln4_deu.append(np.array([2.577273859,2.569834306,2.338145357,2.545390059,2.117084341]))
#gln4_deu.append(np.array([2.560053715,2.703794436,2.563201322,2.040698555,1.920040286,2.128831552,1.821415266]))

#### new fit

gln4_deu.append(np.array([3.2142190085912,3.15209814804528,3.0792667943018,3.12210876709208,2.93574618545434,3.08569309022034,2.64549181980015]))
gln4_deu.append(np.array([2.24974516311868,2.11004157663247,2.05726466618212,2.20110761819385,1.9982787074435,2.01587101092695]))
gln4_deu.append(np.array([2.71371818862284,2.74187744730571,2.67512957487225,2.68764480095352,2.55519199159338,2.5541490560866,2.4123098271655,2.43316853730096]))
gln4_deu.append(np.array([2.4157292706393,2.42423161739034,2.17128680154689,2.34877328997486,2.01080500662101]))
gln4_deu.append(np.array([2.35231165092399,2.53802046547062,2.37959091181784,1.88121979933395,1.7437742925226,2.01761610380323,1.5790495248174,1.95361476093687]))


gln4_dex=[]
#### old fit
#gln4_dex.append(np.array([3.295029553,3.088737932,3.009828569,3.083101549,3.110156188,3.118047124,3.075210613]))
#gln4_dex.append(np.array([2.362113709,2.40379194,2.482980578,2.101624769,2.51423925]))
#gln4_dex.append(np.array([3.009911873,2.993224904,3.058929841,2.854514482,3.0870891,2.886845483]))
#gln4_dex.append(np.array([2.727162174,2.573226247,2.863416272,2.64603378,2.815571322,2.529541727]))
#gln4_dex.append(np.array([2.818032391,2.816982062,2.700395556,2.98818567,2.85689456,2.743459041,2.823284036]))

#### new fit

gln4_dex.append(np.array([3.02899227137774,2.80353694786619,2.7257548612547,2.77986413889747,2.85088256580361,2.89710090712348,2.77535503242724]))
gln4_dex.append(np.array([2.21415599136912,2.2443727084278,2.32877112503999,1.96825443185706,2.30168027526324]))
gln4_dex.append(np.array([2.7919383516308,2.79611009365789,2.82426935234076,2.65218499372325,2.83156990088817,2.73353396325153]))
gln4_dex.append(np.array([2.49521817539892,2.32880095653113,2.63563270381862,2.42241064214426,2.60546958289883,2.30071805084719]))
gln4_dex.append(np.array([2.48087682007817,2.43781333590238,2.36218965637418,2.77811989377932,2.52078931760695,2.48507813560752,2.39580018060894]))


gln4_extra=[]

gln4_extra.append(np.array([2.34250563216932,2.35759076503517,2.17333664217365,2.07636078803601,1.96429980103251,2.24660728752209]))
gln4_extra.append(np.array([1.96725876644925,1.98691066790169,1.91864616811954,1.72626439600621,2.02414584960104,1.96312152403821]))
gln4_extra.append(np.array([2.80937048337355,2.46688900159402,2.63659624331163,2.41372769695958,2.49346965391124,2.407593700271]))
gln4_extra.append(np.array([2.11934912293872,2.08416878035473,2.15026518157314,2.06071521863207,2.03086523098504,1.97542953964057]))
gln4_extra.append(np.array([2.42180350450757,2.16299274372543,2.12141671789698,2.25653880183946,2.214962776011,2.17962315405681]))


glc6_deu=[]


glc6_deu.append(np.array([3.05784580790665,2.76866249157222,2.93039093885555,2.85006223987377,2.60800509360865,2.18708271094409,2.8543464371528]))
glc6_deu.append(np.array([2.7495735503249,2.84788348155594,2.46292248768282,2.67092560534007,2.30769628047592,2.2538845286442]))
glc6_deu.append(np.array([2.78880954511048,2.06918404543725,2.06188349688984,1.794892007156,1.5341581304628,1.88979913827233,2.12028788526912,1.840781169454]))
glc6_deu.append(np.array([2.58258782562846,2.15003093466929,2.29244524274921,2.17978914829793,1.98317237968013]))
glc6_deu.append(np.array([2.39113213758063,2.131979159089,2.24739141671685,2.55061089357548,2.52752844204991,2.09840468414272,2.56844733339069,1.74272509018053]))

glc6_dex=[]


glc6_dex.append(np.array([3.61404883589023,4.53277927919981,4.50008825729064,4.50459736376087,4.175432591434,3.58473964383372,3.44495734325656]))
glc6_dex.append(np.array([3.51868460369577,3.49055179815837,3.51660069217448,3.51555873641383,3.59370541846216]))
glc6_dex.append(np.array([3.26855987822597,3.09856139062201,3.26960281373275,3.04432874426982,3.2331000709957,2.57917950824915]))
glc6_dex.append(np.array([3.08495919476165,3.6538980617659,3.55612794568108,3.46147815245002,3.10784156235597,2.99967037009191]))
glc6_dex.append(np.array([3.37365637006396,4.03326290817112,4.17400697840417,3.77383167423407,3.48394090270927,3.09741987400954,3.76858002982239]))

glc6_extra=[]

glc6_extra.append(np.array([2.43840397681654,2.46749673305784,2.1765691706449,2.44486903375905,2.31987793287053,1.9933925572738]))
glc6_extra.append(np.array([3.41115636790202,3.00053505860635,2.94468228605732,2.47200234059606,2.36650265911456,2.5278551131451]))
glc6_extra.append(np.array([2.47813466218976,2.154055170476,2.00786158273131,2.39328104133096,2.13667551319167,1.95265561253401]))
glc6_extra.append(np.array([2.70995244995485,2.81869169066901,2.38693294077455,2.42424542533333,2.05858307665728,2.03726165690941]))
glc6_extra.append(np.array([1.99149163718305,1.65264702668113,1.52272194596721,1.79192671320646,1.51648554209294,1.494658128533]))


#### Time Points
subject_time_deu=[]

subject_time_deu.append(np.array([21,39,49,60,71,81,93]))
subject_time_deu.append(np.array([24,42,55,66,78,90]))
subject_time_deu.append(np.array([24,42,54,66,77,87,98,109]))
subject_time_deu.append(np.array([20,43,55,67,78]))
subject_time_deu.append(np.array([21,38,49,60,71,82,97,112]))


subject_time_dex=[]

subject_time_dex.append(np.array([23,42,53,64,75,86,97]))
subject_time_dex.append(np.array([19,37,48,60,70]))
subject_time_dex.append(np.array([24,44,54,66,77,88]))
subject_time_dex.append(np.array([22,39,50,61,72,83]))
subject_time_dex.append(np.array([19,37,48,59,72,83,95]))

subject_time_extra=[]

subject_time_extra.append(np.array([26,45,57,69,82,94]))
subject_time_extra.append(np.array([22,41,54,67,78,90]))
subject_time_extra.append(np.array([25,46,57,69,80,92]))
subject_time_extra.append(np.array([27,48,61,75,86,97]))
subject_time_extra.append(np.array([24,45,57,69,80,91]))


################################################ Direct 2H MRSI


### For absolute quantification
TR=290e-3;
TE=1.51e-3;
T1_D2O=350e-3;
T1_glx=150e-3;
T1_glc=70e-3;
T2_D2O=30e-3;
T2_glx=40e-3;
T2_glc=40e-3;

R_D2O=np.exp(-TE/T2_D2O)*(1-np.exp(-TR/T1_D2O));
R_glx=np.exp(-TE/T2_glx)*(1-np.exp(-TR/T1_glx));
R_glc=np.exp(-TE/T2_glc)*(1-np.exp(-TR/T1_glc));

#molar D2O Concentration
mol_D2O=55e3*0.000156;

dGM=0.84; ### bednarik numberss
dWM=0.7;
dCSF=0.97;

##fraction

fGM=[0.78,0.79,0.85,0.79,0.8]
fWM=[0.15,0.16,0.11,0.16,0.15]
fCSF=[0.06,0.03,0.03,0.05,0.03]


abs_factor_glx=[(fGM[i]*dGM*R_D2O+fWM[i]*dWM*R_D2O+fCSF[i]*dCSF*R_D2O)/((1-fCSF[i])*R_glx) * mol_D2O for i in np.arange(0,5)] 
abs_factor_glc=[(fGM[i]*dGM*R_D2O+fWM[i]*dWM*R_D2O+fCSF[i]*dCSF*R_D2O)/((1-fCSF[i])*R_glc) * mol_D2O for i in np.arange(0,5)] 

#### PCC Voxel

glx_2H=[]

glx_2H.append(np.array([0.0000,0.1016,0.1065,0.1669,0.1669,0.2395,0.2419,0.2903,0.2903,0.3871]))
glx_2H.append(np.array([0.0600,0.0407,0.0386,0.1071,0.1436,0.1693,0.1843,0.2571,0.2786,0.3000]))
glx_2H.append(np.array([0.0205,0.0474,0.0789,0.1827,0.1827,0.2932,0.2256,0.3609,0.4286,0.4060]))
glx_2H.append(np.array([0.0000,0.0202,0.0273,0.0841,0.1045,0.1864,0.2273,0.2159,0.2500,0.2955]))
glx_2H.append(np.array([0.0649,0.0000,0.0912,0.1236,0.1338,0.2635,0.3243,0.3243,0.4257]))


glx_2H=[glx_2H[i]*abs_factor_glx[i] for i in np.arange(0,np.size(glx_2H))]

glc_2H=[]

glc_2H.append(np.array([0.0387,0.0871,0.1427,0.1669,0.1621,0.1839,0.2177,0.2226,0.1815,0.2129]))
glc_2H.append(np.array([0.0009,0.0300,0.0621,0.0900,0.0943,0.1093,0.1007,0.0964,0.1071,0.1307]))
glc_2H.append(np.array([0.0609,0.1015,0.1105,0.1489,0.1353,0.1556,0.1128,0.1376,0.1376,0.1421]))
glc_2H.append(np.array([0.0159,0.0318,0.0455,0.0659,0.1273,0.1705,0.1727,0.2091,0.1864,0.2068]))
glc_2H.append(np.array([0.0405,0.0628,0.1297,0.1703,0.1784,0.1804,0.1926,0.1642,0.1703]))

glc_2H=[glc_2H[i]*abs_factor_glc[i] for i in np.arange(0,np.size(glc_2H))]

subject_time_2H=[]
subject_time_2H.append(np.array([10,16.6,23.2,29.8,36.4,43,49.6,56.2,62.8,69.4]))
subject_time_2H.append(np.array([10,16.6,23.2,29.8,36.4,43,49.6,56.2,62.8,69.4]))
subject_time_2H.append(np.array([10,16.6,23.2,29.8,36.4,43,49.6,56.2,62.8,69.4]))
subject_time_2H.append(np.array([10,16.6,23.2,29.8,36.4,43,49.6,56.2,62.8,69.4]))
subject_time_2H.append(np.array([10,16.6,23.2,29.8,36.4,43,49.6,56.2,62.8]))

################# Exponential Fitting glu4 DEU

##Starting Values ### glu4/tcr ratios
m0=[0.5,0.5,0.5,0.5,0.5,0.5]  
tau=[20,20,20,20,20,20]
const=[0,0,0,0,0,0]

### SVS
fitting_glu4=[fit(exp_decrease,subject_time_deu[i],glu4_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_deu))]
m0_deu=[fitting_glu4[i][0][0] for i in np.arange(0,np.size(glu4_deu))]
tau_deu=[fitting_glu4[i][0][1] for i in np.arange(0,np.size(glu4_deu))]
const_deu=[fitting_glu4[i][0][2] for i in np.arange(0,np.size(glu4_deu))]

fit_subject_time_deu=np.arange(0,125,1)
fit_subject_time_dex=np.arange(0,125,1)
fit_subject_time_2H=np.arange(0,125,1)

plot_fit_glu4=[exp_decrease(fit_subject_time_deu,m0_deu[i],tau_deu[i],const_deu[i]) for i in np.arange(0,np.size(m0_deu))]


### SVS
fitting_glu4_extra=[fit(exp_decrease,subject_time_extra[i],glu4_extra[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.shape(glu4_extra)[0])]
m0_extra=[fitting_glu4_extra[i][0][0] for i in np.arange(0,np.shape(glu4_extra)[0])]
tau_extra=[fitting_glu4_extra[i][0][1] for i in np.arange(0,np.shape(glu4_extra)[0])]
const_extra=[fitting_glu4_extra[i][0][2] for i in np.arange(0,np.shape(glu4_extra)[0])]

fit_subject_time_deu=np.arange(0,125,1)
fit_subject_time_dex=np.arange(0,125,1)
fit_subject_time_extra=np.arange(0,125,1)
plot_fit_glu4_extra=[exp_decrease(fit_subject_time_extra,m0_extra[i],tau_extra[i],const_extra[i]) for i in np.arange(0,np.size(m0_extra))]




########################################## Linear Fitting SVS 

################# Linear Fitting glu4 DEU

glu4_deu_slopes=[stats.linregress(subject_time_deu[i],glu4_deu[i])[0] for i in np.arange(0,np.size(glu4_deu))]
glu4_deu_intercept=[stats.linregress(subject_time_deu[i],glu4_deu[i])[1] for i in np.arange(0,np.size(glu4_deu))]
glu4_deu_rvalue=[stats.linregress(subject_time_deu[i],glu4_deu[i])[2] for i in np.arange(0,np.size(glu4_deu))]
glu4_deu_pvalue=[stats.linregress(subject_time_deu[i],glu4_deu[i])[3] for i in np.arange(0,np.size(glu4_deu))]
glu4_deu_r2value=[glu4_deu_rvalue[i]*glu4_deu_rvalue[i] for i in np.arange(0,np.size(glu4_deu_rvalue))]

linreg_plt_glu4_deu=[fit_subject_time_deu*glu4_deu_slopes[i]+glu4_deu_intercept[i] for i in np.arange(0,np.size(glu4_deu_slopes))]

################# Linear Fitting glu4 extra

glu4_extra_slopes=[stats.linregress(subject_time_extra[i],glu4_extra[i])[0] for i in np.arange(0,np.shape(glu4_extra)[0])]
glu4_extra_intercept=[stats.linregress(subject_time_extra[i],glu4_extra[i])[1] for i in np.arange(0,np.shape(glu4_extra)[0])]
glu4_extra_rvalue=[stats.linregress(subject_time_extra[i],glu4_extra[i])[2] for i in np.arange(0,np.shape(glu4_extra)[0])]
glu4_extra_pvalue=[stats.linregress(subject_time_extra[i],glu4_extra[i])[3] for i in np.arange(0,np.shape(glu4_extra)[0])]
glu4_extra_r2value=[glu4_extra_rvalue[i]*glu4_extra_rvalue[i] for i in np.arange(0,np.size(glu4_extra_rvalue))]

linreg_plt_glu4_extra=[fit_subject_time_extra*glu4_extra_slopes[i]+glu4_extra_intercept[i] for i in np.arange(0,np.size(glu4_extra_slopes))]




################# Linear Fitting glu4 DEX

glu4_dex_slopes=[stats.linregress(subject_time_dex[i],glu4_dex[i])[0] for i in np.arange(0,np.size(glu4_dex))]
glu4_dex_intercept=[stats.linregress(subject_time_dex[i],glu4_dex[i])[1] for i in np.arange(0,np.size(glu4_dex))]
glu4_dex_rvalue=[stats.linregress(subject_time_dex[i],glu4_dex[i])[2] for i in np.arange(0,np.size(glu4_dex))]
glu4_dex_pvalue=[stats.linregress(subject_time_dex[i],glu4_dex[i])[3] for i in np.arange(0,np.size(glu4_dex))]
glu4_dex_r2value=[glu4_dex_rvalue[i]*glu4_dex_rvalue[i] for i in np.arange(0,np.size(glu4_dex_rvalue))]

linreg_plt_glu4_dex=[fit_subject_time_dex*glu4_dex_slopes[i]+glu4_dex_intercept[i] for i in np.arange(0,np.size(glu4_dex_slopes))]

################# Linear Fitting gln4 DEU

gln4_deu_slopes=[stats.linregress(subject_time_deu[i],gln4_deu[i])[0] for i in np.arange(0,np.size(gln4_deu))]
gln4_deu_intercept=[stats.linregress(subject_time_deu[i],gln4_deu[i])[1] for i in np.arange(0,np.size(gln4_deu))]
gln4_deu_rvalue=[stats.linregress(subject_time_deu[i],gln4_deu[i])[2] for i in np.arange(0,np.size(gln4_deu))]
gln4_deu_pvalue=[stats.linregress(subject_time_deu[i],gln4_deu[i])[3] for i in np.arange(0,np.size(gln4_deu))]
gln4_deu_r2value=[gln4_deu_rvalue[i]*gln4_deu_rvalue[i] for i in np.arange(0,np.size(gln4_deu_rvalue))]

linreg_plt_gln4_deu=[fit_subject_time_deu*gln4_deu_slopes[i]+gln4_deu_intercept[i] for i in np.arange(0,np.size(gln4_deu_slopes))]


################# Linear Fitting gln4 DEX

gln4_dex_slopes=[stats.linregress(subject_time_dex[i],gln4_dex[i])[0] for i in np.arange(0,np.size(gln4_dex))]
gln4_dex_intercept=[stats.linregress(subject_time_dex[i],gln4_dex[i])[1] for i in np.arange(0,np.size(gln4_dex))]
gln4_dex_rvalue=[stats.linregress(subject_time_dex[i],gln4_dex[i])[2] for i in np.arange(0,np.size(gln4_dex))]
gln4_dex_pvalue=[stats.linregress(subject_time_dex[i],gln4_dex[i])[3] for i in np.arange(0,np.size(gln4_dex))]
gln4_dex_r2value=[gln4_dex_rvalue[i]*gln4_dex_rvalue[i] for i in np.arange(0,np.size(gln4_dex_rvalue))]

linreg_plt_gln4_dex=[fit_subject_time_dex*gln4_dex_slopes[i]+gln4_dex_intercept[i] for i in np.arange(0,np.size(gln4_dex_slopes))]

################# Linear Fitting gln4 extra

gln4_extra_slopes=[stats.linregress(subject_time_extra[i],gln4_extra[i])[0] for i in np.arange(0,np.shape(gln4_extra)[0])]
gln4_extra_intercept=[stats.linregress(subject_time_extra[i],gln4_extra[i])[1] for i in np.arange(0,np.shape(gln4_extra)[0])]
gln4_extra_rvalue=[stats.linregress(subject_time_extra[i],gln4_extra[i])[2] for i in np.arange(0,np.shape(gln4_extra)[0])]
gln4_extra_pvalue=[stats.linregress(subject_time_extra[i],gln4_extra[i])[3] for i in np.arange(0,np.shape(gln4_extra)[0])]
gln4_extra_r2value=[gln4_extra_rvalue[i]*gln4_extra_rvalue[i] for i in np.arange(0,np.size(gln4_extra_rvalue))]

linreg_plt_gln4_extra=[fit_subject_time_extra*gln4_extra_slopes[i]+gln4_extra_intercept[i] for i in np.arange(0,np.size(gln4_extra_slopes))]


################# Linear Fitting glc6 DEU

glc6_deu_slopes=[stats.linregress(subject_time_deu[i],glc6_deu[i])[0] for i in np.arange(0,np.size(glc6_deu))]
glc6_deu_intercept=[stats.linregress(subject_time_deu[i],glc6_deu[i])[1] for i in np.arange(0,np.size(glc6_deu))]
glc6_deu_rvalue=[stats.linregress(subject_time_deu[i],glc6_deu[i])[2] for i in np.arange(0,np.size(glc6_deu))]
glc6_deu_pvalue=[stats.linregress(subject_time_deu[i],glc6_deu[i])[3] for i in np.arange(0,np.size(glc6_deu))]
glc6_deu_r2value=[glc6_deu_rvalue[i]*glc6_deu_rvalue[i] for i in np.arange(0,np.size(glc6_deu_rvalue))]

linreg_plt_glc6_deu=[fit_subject_time_deu*glc6_deu_slopes[i]+glc6_deu_intercept[i] for i in np.arange(0,np.size(glc6_deu_slopes))]


################# Linear Fitting glc6 DEX

glc6_dex_slopes=[stats.linregress(subject_time_dex[i],glc6_dex[i])[0] for i in np.arange(0,np.size(glc6_dex))]
glc6_dex_intercept=[stats.linregress(subject_time_dex[i],glc6_dex[i])[1] for i in np.arange(0,np.size(glc6_dex))]
glc6_dex_rvalue=[stats.linregress(subject_time_dex[i],glc6_dex[i])[2] for i in np.arange(0,np.size(glc6_dex))]
glc6_dex_pvalue=[stats.linregress(subject_time_dex[i],glc6_dex[i])[3] for i in np.arange(0,np.size(glc6_dex))]
glc6_dex_r2value=[glc6_dex_rvalue[i]*glc6_dex_rvalue[i] for i in np.arange(0,np.size(glc6_dex_rvalue))]

linreg_plt_glc6_dex=[fit_subject_time_dex*glc6_dex_slopes[i]+glc6_dex_intercept[i] for i in np.arange(0,np.size(glc6_dex_slopes))]

################# Linear Fitting glc6 extra

glc6_extra_slopes=[stats.linregress(subject_time_extra[i],glc6_extra[i])[0] for i in np.arange(0,np.shape(glc6_extra)[0])]
glc6_extra_intercept=[stats.linregress(subject_time_extra[i],glc6_extra[i])[1] for i in np.arange(0,np.shape(glc6_extra)[0])]
glc6_extra_rvalue=[stats.linregress(subject_time_extra[i],glc6_extra[i])[2] for i in np.arange(0,np.shape(glc6_extra)[0])]
glc6_extra_pvalue=[stats.linregress(subject_time_extra[i],glc6_extra[i])[3] for i in np.arange(0,np.shape(glc6_extra)[0])]
glc6_extra_r2value=[glc6_extra_rvalue[i]*glc6_extra_rvalue[i] for i in np.arange(0,np.size(glc6_extra_rvalue))]

linreg_plt_glc6_extra=[fit_subject_time_extra*glc6_extra_slopes[i]+glc6_extra_intercept[i] for i in np.arange(0,np.size(glc6_extra_slopes))]


################# Linear Fitting glx_2H

glx_2H_slopes=[stats.linregress(subject_time_2H[i],glx_2H[i])[0] for i in np.arange(0,np.size(glx_2H))]
glx_2H_intercept=[stats.linregress(subject_time_2H[i],glx_2H[i])[1] for i in np.arange(0,np.size(glx_2H))]
glx_2H_rvalue=[stats.linregress(subject_time_2H[i],glx_2H[i])[2] for i in np.arange(0,np.size(glx_2H))]
glx_2H_pvalue=[stats.linregress(subject_time_2H[i],glx_2H[i])[3] for i in np.arange(0,np.size(glx_2H))]
glx_2H_r2value=[glx_2H_rvalue[i]*glx_2H_rvalue[i] for i in np.arange(0,np.size(glx_2H_rvalue))]

linreg_plt_glx_2H=[fit_subject_time_2H*glx_2H_slopes[i]+glx_2H_intercept[i] for i in np.arange(0,np.size(glx_2H_slopes))]

################# Linear Fitting glx_2H

glc_2H_slopes=[stats.linregress(subject_time_2H[i],glc_2H[i])[0] for i in np.arange(0,np.size(glc_2H))]
glc_2H_intercept=[stats.linregress(subject_time_2H[i],glc_2H[i])[1] for i in np.arange(0,np.size(glc_2H))]
glc_2H_rvalue=[stats.linregress(subject_time_2H[i],glc_2H[i])[2] for i in np.arange(0,np.size(glc_2H))]
glc_2H_pvalue=[stats.linregress(subject_time_2H[i],glc_2H[i])[3] for i in np.arange(0,np.size(glc_2H))]
glc_2H_r2value=[glc_2H_rvalue[i]*glc_2H_rvalue[i] for i in np.arange(0,np.size(glc_2H_rvalue))]

linreg_plt_glc_2H=[fit_subject_time_2H*glc_2H_slopes[i]+glc_2H_intercept[i] for i in np.arange(0,np.size(glc_2H_slopes))]


########################################################### ALL PLOTS SVS
#### glu4_deu exponential
plt.figure('glu4_deu')
plt.title('glu4_deu')

[plt.plot(subject_time_deu[i],glu4_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' tau='+str(np.round(tau_deu[i],2))+' min') for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,plot_fit_glu4[i]) for i in np.arange(0,np.size(m0_deu))]
plt.ylim(5,8.5)
plt.text(50,7.5,'Intersubject COV [tau] = '+str(np.round(np.std(tau_deu)/np.mean(tau_deu)*100,2))+'%')
plt.legend()


plt.savefig('glu4_svs_deu_exp.eps',format='eps')

#### glu4_deu linear

[plt.plot(subject_time_deu[i],glu4_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glu4_deu_r2value[i],2))+' p='+str(np.round(glu4_deu_pvalue[i],2))+ ' slope='+str(np.round(glu4_deu_slopes[i],4))) for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,linreg_plt_glu4_deu[i]) for i in np.arange(0,np.size(glu4_deu_r2value))]
plt.legend()
plt.ylim(5,8.5)



plt.savefig('glu4_svs_deu.pdf',format='pdf')



########################################################### ALL PLOTS SVS
#### glu4_extra exponential
plt.figure('glu4_extra')
plt.title('glu4_extra')

[plt.plot(subject_time_extra[i],glu4_extra[i],ls='none',marker='o',label='sub '+str(i+1)+' tau='+str(np.round(tau_extra[i],2))+' min') for i in np.arange(0,np.shape(subject_time_extra)[0])]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_extra,plot_fit_glu4_extra[i]) for i in np.arange(0,np.size(m0_extra))]
plt.ylim(5,8.5)
plt.text(50,7.5,'Intersubject COV [tau] = '+str(np.round(np.std(tau_extra)/np.mean(tau_extra)*100,2))+'%')
plt.legend()

plt.savefig('glu4_svs_extra_exp.pdf',format='pdf')

#### glu4_extra linear

[plt.plot(subject_time_extra[i],glu4_extra[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glu4_extra_r2value[i],2))+' p='+str(np.round(glu4_extra_pvalue[i],2))+ ' slope='+str(np.round(glu4_extra_slopes[i],4))) for i in np.arange(0,np.shape(subject_time_extra)[0])]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_extra,linreg_plt_glu4_extra[i]) for i in np.arange(0,np.size(glu4_extra_r2value))]
plt.legend()
plt.ylim(5,8.5)

plt.savefig('glu4_svs_extra.pdf',format='pdf')



#### glu4_dex_linear
plt.figure('glu4_dex')
plt.title('glu4_dex')

[plt.plot(subject_time_dex[i],glu4_dex[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glu4_dex_r2value[i],2))+' p='+str(np.round(glu4_dex_pvalue[i],2))+ 'slope='+str(np.round(glu4_dex_slopes[i],4))) for i in np.arange(0,np.size(subject_time_dex))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_dex,linreg_plt_glu4_dex[i]) for i in np.arange(0,np.size(glu4_dex_r2value))]
plt.legend()
plt.ylim(5,8.5)

plt.savefig('glu4_svs_dex.eps',format='eps')

#### gln4_deu_linear
plt.figure('gln4_deu')
plt.title('gln4_deu')

[plt.plot(subject_time_deu[i],gln4_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(gln4_deu_r2value[i],2))+' p='+str(np.round(gln4_deu_pvalue[i],2))+ ' slope='+str(np.round(gln4_deu_slopes[i],4))) for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,linreg_plt_gln4_deu[i]) for i in np.arange(0,np.size(gln4_deu_r2value))]
plt.legend()
plt.ylim(1.2,4)

plt.savefig('gln4_svs_deu.pdf',format='pdf')


#### gln4_extra_linear
plt.figure('gln4_extra')
plt.title('gln4_extra')

[plt.plot(subject_time_extra[i],gln4_extra[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(gln4_extra_r2value[i],2))+' p='+str(np.round(gln4_extra_pvalue[i],2))+ ' slope='+str(np.round(gln4_extra_slopes[i],4))) for i in np.arange(0,np.shape(subject_time_extra)[0])]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_extra,linreg_plt_gln4_extra[i]) for i in np.arange(0,np.size(gln4_extra_r2value))]
plt.legend()
plt.ylim(1.2,4)

plt.savefig('gln4_svs_extra.pdf',format='pdf')


#### gln4_dex_linear

plt.figure('gln4_dex')
plt.title('gln4_dex')

[plt.plot(subject_time_dex[i],gln4_dex[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(gln4_dex_r2value[i],2))+' p='+str(np.round(gln4_dex_pvalue[i],2))+ ' slope='+str(np.round(gln4_dex_slopes[i],4))) for i in np.arange(0,np.size(subject_time_dex))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_dex,linreg_plt_gln4_dex[i]) for i in np.arange(0,np.size(gln4_dex_r2value))]
plt.legend()
plt.ylim(1.2,4)

plt.savefig('gln4_svs_dex.eps',format='eps')


#### glc6_deu_linear
plt.figure('glc6_deu')
plt.title('glc6_deu')

[plt.plot(subject_time_deu[i],glc6_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glc6_deu_r2value[i],2))+' p='+str(np.round(glc6_deu_pvalue[i],2))+ ' slope='+str(np.round(glc6_deu_slopes[i],4))) for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,linreg_plt_glc6_deu[i]) for i in np.arange(0,np.size(glc6_deu_r2value))]
plt.legend()
plt.ylim(1.2,5)

plt.savefig('glc6_svs_deu.eps',format='eps')


#### glc6_extra_linear
plt.figure('glc6_extra')
plt.title('glc6_extra')

[plt.plot(subject_time_extra[i],glc6_extra[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glc6_extra_r2value[i],2))+' p='+str(np.round(glc6_extra_pvalue[i],2))+ ' slope='+str(np.round(glc6_extra_slopes[i],4))) for i in np.arange(0,np.shape(subject_time_extra)[0])]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_extra,linreg_plt_glc6_extra[i]) for i in np.arange(0,np.size(glc6_extra_r2value))]
plt.legend()
plt.ylim(1.2,5)

plt.savefig('glc6_svs_extra.eps',format='eps')



#### glc6_dex_linear

plt.figure('glc6_dex')
plt.title('glc6_dex')

[plt.plot(subject_time_dex[i],glc6_dex[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glc6_dex_r2value[i],2))+' p='+str(np.round(glc6_dex_pvalue[i],2))+ ' slope='+str(np.round(glc6_dex_slopes[i],4))) for i in np.arange(0,np.size(subject_time_dex))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_dex,linreg_plt_glc6_dex[i]) for i in np.arange(0,np.size(glc6_dex_r2value))]
plt.legend()
plt.ylim(1.2,5)

plt.savefig('glc6_svs_dex.eps',format='eps')


#### glx_2H_linear

plt.figure('glx_2H')
plt.title('glx_2H')

[plt.plot(subject_time_2H[i],glx_2H[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glx_2H_r2value[i],2))+' p='+str(glx_2H_pvalue[i])+ ' slope='+str(np.round(glx_2H_slopes[i],4))) for i in np.arange(0,np.size(subject_time_2H))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_2H,linreg_plt_glx_2H[i]) for i in np.arange(0,np.size(glx_2H_r2value))]
plt.legend()
#plt.ylim(-2,11)
#plt.ylim(-0.3,3.5)

plt.savefig('glx_PCC_2H.pdf',format='pdf')

#### glc_2H_linear

plt.figure('glc_2H')
plt.title('glc_2H')

[plt.plot(subject_time_2H[i],glc_2H[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glc_2H_r2value[i],2))+' p='+str(glc_2H_pvalue[i])+ ' slope='+str(np.round(glc_2H_slopes[i],4))) for i in np.arange(0,np.size(subject_time_2H))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_2H,linreg_plt_glc_2H[i]) for i in np.arange(0,np.size(glc_2H_r2value))]
plt.legend()
plt.ylim(-0.3,2)

plt.savefig('glc_PCC_2H.pdf',format='pdf')


###################################### MRSI


glu4_GM_deu=[]


glu4_GM_deu.append(np.array([0.552410157253131,0.512011808473785,0.491329914427083,0.484089421862397,0.478732139924454,0.479995729168255,0.460140380235207]))
glu4_GM_deu.append(np.array([0.513499789886539,0.499711087035031,0.470204227236244,0.481128879855948,0.474579497473524,0.472186150984452]))
glu4_GM_deu.append(np.array([0.51779996506957,0.485622577877703,0.479469399172592,0.462508937082936,0.466406635421806,0.44481780167264,0.44618496427095,0.442449856531019]))
glu4_GM_deu.append(np.array([0.510889450300549,0.473980534501069,0.465039907477476,0.455054811205847,0.437326586063629,0.437633592080099,0.437884631164355,0.434889700906753]))
glu4_GM_deu.append(np.array([0.539985772297478,0.52708562264263,0.522318909447299,0.514248782540732,0.502581523772467,0.498119378577269,0.495225531787315,0.490649397809115,0.489530441510711,0.482381751394495,0.471052552597617,0.468255401095264,0.460649084183253,0.455208176762363]))

glu4_WM_deu=[]


glu4_WM_deu.append(np.array([0.51310270832656,0.483450514673496,0.463658781920244,0.452808603134653,0.450231819213436,0.450044252138853,0.441160683788487]))
glu4_WM_deu.append(np.array([0.464895635673624,0.44534833575097,0.432814984294412,0.442445488866631,0.426016106110848,0.421797958949662]))
glu4_WM_deu.append(np.array([0.51779996506957,0.485622577877703,0.479469399172592,0.462508937082936,0.466406635421806,0.44481780167264,0.44618496427095,0.442449856531019]))
glu4_WM_deu.append(np.array([0.470976898210203,0.445782943054398,0.431234241341439,0.419211136890951,0.416184971098266,0.409884653578795,0.410188351028733,0.40027131039974]))
glu4_WM_deu.append(np.array([0.505997695067119,0.488940666841691,0.485848512321767,0.477209625173531,0.467413357966448,0.459606056633456,0.463672477756985,0.454533317757789,0.46193819736522,0.446770682750943,0.432512007368906,0.429864629183623,0.426903185426755,0.418728644294416]))


glu4_GM_extra=[]


glu4_GM_extra.append(np.array([0.532938342534823,0.499257681867942,0.481753902387872,0.46933537463977,0.453021535524667,0.439610427480026]))
glu4_GM_extra.append(np.array([0.532636560836127,0.50201242453408,0.490095451670404,0.480086997853861,0.466342697760014,0.454143724410561]))
glu4_GM_extra.append(np.array([0.539173133412249,0.510751489753763,0.495262666899889,0.492578403810578,0.474273390161241,0.464176004212867]))
glu4_GM_extra.append(np.array([0.509124161209203,0.490773474744228,0.47335654036685,0.467985936173939,0.461708720672864,0.443388037977219]))
glu4_GM_extra.append(np.array([0.511955397172089,0.481172440262122,0.470773819767108,0.454699684929359,0.458159357272439,0.447925765755445]))

glu4_WM_extra=[]


glu4_WM_extra.append(np.array([0.49811669967854,0.465976667858187,0.452621507224754,0.44052008222128,0.425100889110407,0.412945157821494]))
glu4_WM_extra.append(np.array([0.445586086628647,0.432611855919884,0.418205361853149,0.415118596130359,0.402553505062156,0.391819583897265]))
glu4_WM_extra.append(np.array([0.499496137050722,0.480993206082174,0.464800637958533,0.461437678751178,0.446536576048895,0.446548574152038]))
glu4_WM_extra.append(np.array([0.477268304366701,0.459076677693762,0.447011551983928,0.438942279171857,0.425426568035934,0.421136554901006]))
glu4_WM_extra.append(np.array([0.474471409892528,0.445391766686068,0.435507765830346,0.423778171689035,0.419863767607319,0.41437125748503]))

glu4_GM_dex=[]


glu4_GM_dex.append(np.array([0.530463090032873,0.521588323000863,0.51603144842364,0.51603144842364,0.528055450773293,0.530169842714972,0.526752345659544]))
glu4_GM_dex.append(np.array([0.57980395783244,0.57270126802055,0.577864616517295,0.568044397763079,0.582190396442699]))
glu4_GM_dex.append(np.array([0.582839414028772,0.643462793790043,0.578633701633008,0.637334705588822,0.609523669681663,0.580487358844237]))
glu4_GM_dex.append(np.array([0.503259178266538,0.5077753448348,0.504539595756434,0.510097347354328,0.510204081632653,0.508957643440068,0.5072,0.510359609430817]))
glu4_GM_dex.append(np.array([0.522980723881392,0.526959968531574,0.533286325244876,0.530381891908249,0.532018636905407,0.525660894031265,0.532814456159404,0.536276433061987]))

glu4_WM_dex=[]


glu4_WM_dex.append(np.array([0.531618922514694,0.52457831721057,0.521564510215645,0.521564510215645,0.527966241932815,0.529112297350732,0.525121769999174]))
glu4_WM_dex.append(np.array([0.500055098153425,0.497994884587451,0.501205687085056,0.495978510163687,0.495902029924712]))
glu4_WM_dex.append(np.array([0.525081753817824,0.580730515875953,0.518578609149329,0.555951754962534,0.54910541098335,0.563272727272727]))
glu4_WM_dex.append(np.array([0.472649693318901,0.469542157398019,0.479521806508745,0.467002721425037,0.476009434840632,0.477035959195678,0.479811097992916,0.479506800086744]))
glu4_WM_dex.append(np.array([0.481573394108213,0.485949049519287,0.483391342547224,0.49136857690884,0.486995866573679,0.485919602514297,0.489636567751881,0.507786892206314]))


#### Time Points
mrsi_subject_time_deu=[]

mrsi_subject_time_deu.append(np.array([17,34,45,56,67,77,89]))
mrsi_subject_time_deu.append(np.array([24,42,55,66,78,90]))
mrsi_subject_time_deu.append(np.array([24,42,54,66,77,87,98,109]))
mrsi_subject_time_deu.append(np.array([21,38,49,60,71,82,97,112]))
mrsi_subject_time_deu.append(np.array([5,10,15,20,25,30,35,40,45,50,55,60,73,78]))


mrsi_subject_time_dex=[]

mrsi_subject_time_dex.append(np.array([19,38,49,60,71,82,93]))
mrsi_subject_time_dex.append(np.array([15,33,44,55,66]))
mrsi_subject_time_dex.append(np.array([21,39,50,61,72,84]))
mrsi_subject_time_dex.append(np.array([15,32,44,55,68,79,91,103]))
mrsi_subject_time_dex.append(np.array([20,40,47,63,74,85,96,108]))

mrsi_subject_time_extra=[]

mrsi_subject_time_extra.append(np.array([21,40,53,65,78,88]))
mrsi_subject_time_extra.append(np.array([17,37,49,62,74,85]))
mrsi_subject_time_extra.append(np.array([21,41,53,65,75,87]))
mrsi_subject_time_extra.append(np.array([22,44,57,71,82,93]))
mrsi_subject_time_extra.append(np.array([20,41,54,65,76,87]))

#################################### direct 2H MRSI
glx_GM_mrsi=[]

glx_GM_mrsi.append(np.array([0.041739074810495,0.0508545049415227,0.0777850037847692,0.127679989765131,0.156369606703839,0.198268601341194,0.238152606160112,0.26700214292568,0.320959092508289,0.353262897534036]));
glx_GM_mrsi.append(np.array([0.0417035696407255,0.0358696907851669,0.0553352073626247,0.077862066309831,0.118525934768378,0.14746428433902,0.176806962147175,0.193182255766491,0.22772343948554,0.251867611382803]));
glx_GM_mrsi.append(np.array([0.0388691103849326,0.0507083720493198,0.0740736873975549,0.115542424020964,0.166626645646931,0.198730463652214,0.231867867993276,0.265036593132393,0.309950617541735,0.327208377269453]));
glx_GM_mrsi.append(np.array([0.042965327325072,0.044529978720741,0.0487545374890474,0.070409312805107,0.10958818375266,0.147421454499937,0.182813869069971,0.20428088621855,0.243960445612717,0.277256227312555]));
glx_GM_mrsi.append(np.array([0.0340065044062107,0.0437604909777591,0.0747665757448594,0.102808959295006,0.142291754091481,0.194940725975661,0.24326741502308,0.270063994964331,0.315725975660932]));

glx_WM_mrsi=[]

glx_WM_mrsi.append(np.array([0.0325885762363279,0.0465270636627092,0.0681219033373843,0.116135365055623,0.141039543797326,0.179713938487426,0.229942974665794,0.266654202112742,0.30524446106385,0.33155090212209]));
glx_WM_mrsi.append(np.array([0.0318796683693113,0.0365255781719243,0.0574450062886624,0.0740008727123386,0.111039811083447,0.138247901640186,0.16466028388819,0.190944326086399,0.224928771272363,0.24179265381555]));
glx_WM_mrsi.append(np.array([0.0333202884776653,0.0436909487337172,0.0758372002012635,0.0973612120534466,0.150360597081679,0.180382400626153,0.212109353161514,0.251747078884106,0.289707608877956,0.317297478615754]));
glx_WM_mrsi.append(np.array([0.0284130347422643,0.0322215321667336,0.0345230701786432,0.0612099514119753,0.0930204946480108,0.135379753771965,0.172752347203449,0.200617396704782,0.231797756913747,0.267416797574252]));
glx_WM_mrsi.append(np.array([0.0285342785814876,0.0398668544360517,0.0677410062732045,0.0970914415567789,0.135952502880553,0.177334048137242,0.225136826270644,0.255601075406478,0.29338432979132]));


glc_GM_mrsi=[]

glc_GM_mrsi.append(np.array([0.0402678124033818,0.0695651246841583,0.109161273814728,0.158416580487649,0.183140185720226,0.205688881307504,0.218546435387059,0.229900743094129,0.227278057933623,0.221776815889634]));
glc_GM_mrsi.append(np.array([0.0176749200970388,0.03844006315222,0.0612268473949709,0.083609303400208,0.105962878817051,0.107926758827833,0.113240787092302,0.118554815356772,0.12551503715969,0.124533097154299]));
glc_GM_mrsi.append(np.array([0.0405291126817911,0.0803691678063957,0.112003173840869,0.12991867032772,0.137968115427581,0.14329265109675,0.143543217481181,0.138218681812013,0.142760197529833,0.138250002610067]));
glc_GM_mrsi.append(np.array([0.0270058830892477,0.0363312054074352,0.0640255351107773,0.0997621729878583,0.139911127800726,0.166854424834147,0.191669795969458,0.202121667292527,0.210007510326699,0.214983101764927]));
glc_GM_mrsi.append(np.array([0.0316906210658833,0.0712048887956358,0.121401594628619,0.157055182543013,0.178674464960134,0.192110784725136,0.193427402433907,0.176930339907679,0.171566827528326]));

glc_WM_mrsi=[]

glc_WM_mrsi.append(np.array([0.0321959427876975,0.0614751799569973,0.0945685706272787,0.146115733383192,0.171833224268487,0.190118724876133,0.207787230064504,0.216621482658689,0.216986070860989,0.211180704870524]));
glc_WM_mrsi.append(np.array([0.0187376472702071,0.0360892220026181,0.0593701070356015,0.0786211145049924,0.0957160091378116,0.103108396006058,0.111527503272671,0.117200133473652,0.120023614569162,0.120254626658795]));
glc_WM_mrsi.append(np.array([0.034075026555599,0.069519762956337,0.102560518812545,0.117543467322637,0.12919997763739,0.138145021524012,0.138033208475429,0.134203611561469,0.134762676804383,0.131464191871191]));
glc_WM_mrsi.append(np.array([0.0252621196069119,0.0336462937931538,0.0559492930990392,0.09247250940708,0.12907792350126,0.159162313228364,0.180396741314434,0.193548387096774,0.208919373104884,0.205604062397253]));
glc_WM_mrsi.append(np.array([0.0301353859941109,0.0649780757905518,0.115817436947894,0.146610549225451,0.169232492638587,0.188534598642939,0.184897900396876,0.173176449878377,0.176577902957368]));


################# Linear Fitting glu4 mrsi DEU

glu4_GM_deu_slopes=[stats.linregress(mrsi_subject_time_deu[i],glu4_GM_deu[i])[0] for i in np.arange(0,np.size(glu4_GM_deu))]
glu4_GM_deu_intercept=[stats.linregress(mrsi_subject_time_deu[i],glu4_GM_deu[i])[1] for i in np.arange(0,np.size(glu4_GM_deu))]
glu4_GM_deu_rvalue=[stats.linregress(mrsi_subject_time_deu[i],glu4_GM_deu[i])[2] for i in np.arange(0,np.size(glu4_GM_deu))]
glu4_GM_deu_pvalue=[stats.linregress(mrsi_subject_time_deu[i],glu4_GM_deu[i])[3] for i in np.arange(0,np.size(glu4_GM_deu))]
glu4_GM_deu_r2value=[glu4_GM_deu_rvalue[i]*glu4_GM_deu_rvalue[i] for i in np.arange(0,np.size(glu4_GM_deu_rvalue))]

linreg_plt_glu4_GM_deu=[fit_subject_time_deu*glu4_GM_deu_slopes[i]+glu4_GM_deu_intercept[i] for i in np.arange(0,np.size(glu4_GM_deu_slopes))]


glu4_WM_deu_slopes=[stats.linregress(mrsi_subject_time_deu[i],glu4_WM_deu[i])[0] for i in np.arange(0,np.size(glu4_WM_deu))]
glu4_WM_deu_intercept=[stats.linregress(mrsi_subject_time_deu[i],glu4_WM_deu[i])[1] for i in np.arange(0,np.size(glu4_WM_deu))]
glu4_WM_deu_rvalue=[stats.linregress(mrsi_subject_time_deu[i],glu4_WM_deu[i])[2] for i in np.arange(0,np.size(glu4_WM_deu))]
glu4_WM_deu_pvalue=[stats.linregress(mrsi_subject_time_deu[i],glu4_WM_deu[i])[3] for i in np.arange(0,np.size(glu4_WM_deu))]
glu4_WM_deu_r2value=[glu4_WM_deu_rvalue[i]*glu4_WM_deu_rvalue[i] for i in np.arange(0,np.size(glu4_WM_deu_rvalue))]

linreg_plt_glu4_WM_deu=[fit_subject_time_deu*glu4_WM_deu_slopes[i]+glu4_WM_deu_intercept[i] for i in np.arange(0,np.size(glu4_WM_deu_slopes))]



################# Linear Fitting glu4 mrsi extra

glu4_GM_extra_slopes=[stats.linregress(mrsi_subject_time_extra[i],glu4_GM_extra[i])[0] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
glu4_GM_extra_intercept=[stats.linregress(mrsi_subject_time_extra[i],glu4_GM_extra[i])[1] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
glu4_GM_extra_rvalue=[stats.linregress(mrsi_subject_time_extra[i],glu4_GM_extra[i])[2] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
glu4_GM_extra_pvalue=[stats.linregress(mrsi_subject_time_extra[i],glu4_GM_extra[i])[3] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
glu4_GM_extra_r2value=[glu4_GM_extra_rvalue[i]*glu4_GM_extra_rvalue[i] for i in np.arange(0,np.shape(glu4_GM_extra_rvalue)[0])]

linreg_plt_glu4_GM_extra=[fit_subject_time_extra*glu4_GM_extra_slopes[i]+glu4_GM_extra_intercept[i] for i in np.arange(0,np.shape(glu4_GM_extra_slopes)[0])]


glu4_WM_extra_slopes=[stats.linregress(mrsi_subject_time_extra[i],glu4_WM_extra[i])[0] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
glu4_WM_extra_intercept=[stats.linregress(mrsi_subject_time_extra[i],glu4_WM_extra[i])[1] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
glu4_WM_extra_rvalue=[stats.linregress(mrsi_subject_time_extra[i],glu4_WM_extra[i])[2] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
glu4_WM_extra_pvalue=[stats.linregress(mrsi_subject_time_extra[i],glu4_WM_extra[i])[3] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
glu4_WM_extra_r2value=[glu4_WM_extra_rvalue[i]*glu4_WM_extra_rvalue[i] for i in np.arange(0,np.shape(glu4_WM_extra_rvalue)[0])]

linreg_plt_glu4_WM_extra=[fit_subject_time_extra*glu4_WM_extra_slopes[i]+glu4_WM_extra_intercept[i] for i in np.arange(0,np.shape(glu4_WM_extra_slopes)[0])]


################# Linear Fitting glu4 mrsi DEX

glu4_GM_dex_slopes=[stats.linregress(mrsi_subject_time_dex[i],glu4_GM_dex[i])[0] for i in np.arange(0,np.size(glu4_GM_dex))]
glu4_GM_dex_intercept=[stats.linregress(mrsi_subject_time_dex[i],glu4_GM_dex[i])[1] for i in np.arange(0,np.size(glu4_GM_dex))]
glu4_GM_dex_rvalue=[stats.linregress(mrsi_subject_time_dex[i],glu4_GM_dex[i])[2] for i in np.arange(0,np.size(glu4_GM_dex))]
glu4_GM_dex_pvalue=[stats.linregress(mrsi_subject_time_dex[i],glu4_GM_dex[i])[3] for i in np.arange(0,np.size(glu4_GM_dex))]
glu4_GM_dex_r2value=[glu4_GM_dex_rvalue[i]*glu4_GM_dex_rvalue[i] for i in np.arange(0,np.size(glu4_GM_dex_rvalue))]

linreg_plt_glu4_GM_dex=[fit_subject_time_dex*glu4_GM_dex_slopes[i]+glu4_GM_dex_intercept[i] for i in np.arange(0,np.size(glu4_GM_dex_slopes))]


glu4_WM_dex_slopes=[stats.linregress(mrsi_subject_time_dex[i],glu4_WM_dex[i])[0] for i in np.arange(0,np.size(glu4_WM_dex))]
glu4_WM_dex_intercept=[stats.linregress(mrsi_subject_time_dex[i],glu4_WM_dex[i])[1] for i in np.arange(0,np.size(glu4_WM_dex))]
glu4_WM_dex_rvalue=[stats.linregress(mrsi_subject_time_dex[i],glu4_WM_dex[i])[2] for i in np.arange(0,np.size(glu4_WM_dex))]
glu4_WM_dex_pvalue=[stats.linregress(mrsi_subject_time_dex[i],glu4_WM_dex[i])[3] for i in np.arange(0,np.size(glu4_WM_dex))]
glu4_WM_dex_r2value=[glu4_WM_dex_rvalue[i]*glu4_WM_dex_rvalue[i] for i in np.arange(0,np.size(glu4_WM_dex_rvalue))]

linreg_plt_glu4_WM_dex=[fit_subject_time_dex*glu4_WM_dex_slopes[i]+glu4_WM_dex_intercept[i] for i in np.arange(0,np.size(glu4_WM_dex_slopes))]


#### pooled mrsi

pooled_glu4_GM_deu=[]
[[pooled_glu4_GM_deu.append(glu4_GM_deu[y][x]) for x in np.arange(0,np.size(glu4_GM_deu[y]))] for y in np.arange(0,np.size(glu4_GM_deu))]

pooled_glu4_WM_deu=[]
[[pooled_glu4_WM_deu.append(glu4_WM_deu[y][x]) for x in np.arange(0,np.size(glu4_WM_deu[y]))] for y in np.arange(0,np.size(glu4_WM_deu))]

pooled_mrsi_subject_time_deu=[]
[[pooled_mrsi_subject_time_deu.append(mrsi_subject_time_deu[y][x]) for x in np.arange(0,np.size(mrsi_subject_time_deu[y]))] for y in np.arange(0,np.size(mrsi_subject_time_deu))]

pooled_glu4_GM_extra=[]
[[pooled_glu4_GM_extra.append(glu4_GM_extra[y][x]) for x in np.arange(0,np.size(glu4_GM_extra[y]))] for y in np.arange(0,np.shape(glu4_GM_extra)[0])]

pooled_glu4_WM_extra=[]
[[pooled_glu4_WM_extra.append(glu4_WM_extra[y][x]) for x in np.arange(0,np.size(glu4_WM_extra[y]))] for y in np.arange(0,np.shape(glu4_WM_extra)[0])]

pooled_mrsi_subject_time_extra=[]
[[pooled_mrsi_subject_time_extra.append(mrsi_subject_time_extra[y][x]) for x in np.arange(0,np.size(mrsi_subject_time_extra[y]))] for y in np.arange(0,np.shape(mrsi_subject_time_extra)[0])]

pooled_glu4_GM_dex=[]
[[pooled_glu4_GM_dex.append(glu4_GM_dex[y][x]) for x in np.arange(0,np.size(glu4_GM_dex[y]))] for y in np.arange(0,np.shape(glu4_GM_dex)[0])]

pooled_glu4_WM_dex=[]
[[pooled_glu4_WM_dex.append(glu4_WM_dex[y][x]) for x in np.arange(0,np.size(glu4_WM_dex[y]))] for y in np.arange(0,np.shape(glu4_WM_dex)[0])]

pooled_mrsi_subject_time_dex=[]
[[pooled_mrsi_subject_time_dex.append(mrsi_subject_time_dex[y][x]) for x in np.arange(0,np.size(mrsi_subject_time_dex[y]))] for y in np.arange(0,np.shape(mrsi_subject_time_dex)[0])]

#### 2H MRSI 

pooled_glx_GM_mrsi=[]
[[pooled_glx_GM_mrsi.append(glx_GM_mrsi[y][x]) for x in np.arange(0,np.size(glx_GM_mrsi[y]))] for y in np.arange(0,np.shape(glx_GM_mrsi)[0])]

pooled_glx_WM_mrsi=[]
[[pooled_glx_WM_mrsi.append(glx_WM_mrsi[y][x]) for x in np.arange(0,np.size(glx_WM_mrsi[y]))] for y in np.arange(0,np.shape(glx_WM_mrsi)[0])]

pooled_glc_GM_mrsi=[]
[[pooled_glc_GM_mrsi.append(glc_GM_mrsi[y][x]) for x in np.arange(0,np.size(glc_GM_mrsi[y]))] for y in np.arange(0,np.shape(glc_GM_mrsi)[0])]

pooled_glc_WM_mrsi=[]
[[pooled_glc_WM_mrsi.append(glc_WM_mrsi[y][x]) for x in np.arange(0,np.size(glc_WM_mrsi[y]))] for y in np.arange(0,np.shape(glc_WM_mrsi)[0])]


pooled_subject_time_2H=[]
[[pooled_subject_time_2H.append(subject_time_2H[y][x]) for x in np.arange(0,np.size(subject_time_2H[y]))] for y in np.arange(0,np.shape(subject_time_2H)[0])]


#### pooled fitting

pooled_glu4_GM_deu_slopes=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_GM_deu)[0]
pooled_glu4_GM_deu_intercept=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_GM_deu)[1]
pooled_glu4_GM_deu_rvalue=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_GM_deu)[2]
pooled_glu4_GM_deu_pvalue=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_GM_deu)[3]
pooled_glu4_GM_deu_r2value=pooled_glu4_GM_deu_rvalue*pooled_glu4_GM_deu_rvalue

pooled_linreg_plt_glu4_GM_deu=fit_subject_time_deu*pooled_glu4_GM_deu_slopes+pooled_glu4_GM_deu_intercept

pooled_glu4_WM_deu_slopes=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_WM_deu)[0]
pooled_glu4_WM_deu_intercept=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_WM_deu)[1]
pooled_glu4_WM_deu_rvalue=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_WM_deu)[2]
pooled_glu4_WM_deu_pvalue=stats.linregress(pooled_mrsi_subject_time_deu,pooled_glu4_WM_deu)[3]
pooled_glu4_WM_deu_r2value=pooled_glu4_WM_deu_rvalue*pooled_glu4_WM_deu_rvalue

pooled_linreg_plt_glu4_WM_deu=fit_subject_time_deu*pooled_glu4_WM_deu_slopes+pooled_glu4_WM_deu_intercept

pooled_glu4_GM_extra_slopes=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_GM_extra)[0]
pooled_glu4_GM_extra_intercept=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_GM_extra)[1]
pooled_glu4_GM_extra_rvalue=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_GM_extra)[2]
pooled_glu4_GM_extra_pvalue=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_GM_extra)[3]
pooled_glu4_GM_extra_r2value=pooled_glu4_GM_extra_rvalue*pooled_glu4_GM_extra_rvalue

pooled_linreg_plt_glu4_GM_extra=fit_subject_time_extra*pooled_glu4_GM_extra_slopes+pooled_glu4_GM_extra_intercept

pooled_glu4_WM_extra_slopes=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_WM_extra)[0]
pooled_glu4_WM_extra_intercept=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_WM_extra)[1]
pooled_glu4_WM_extra_rvalue=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_WM_extra)[2]
pooled_glu4_WM_extra_pvalue=stats.linregress(pooled_mrsi_subject_time_extra,pooled_glu4_WM_extra)[3]
pooled_glu4_WM_extra_r2value=pooled_glu4_WM_extra_rvalue*pooled_glu4_WM_extra_rvalue

pooled_linreg_plt_glu4_WM_extra=fit_subject_time_extra*pooled_glu4_WM_extra_slopes+pooled_glu4_WM_extra_intercept


pooled_glu4_GM_dex_slopes=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_GM_dex)[0]
pooled_glu4_GM_dex_intercept=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_GM_dex)[1]
pooled_glu4_GM_dex_rvalue=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_GM_dex)[2]
pooled_glu4_GM_dex_pvalue=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_GM_dex)[3]
pooled_glu4_GM_dex_r2value=pooled_glu4_GM_dex_rvalue*pooled_glu4_GM_dex_rvalue

pooled_linreg_plt_glu4_GM_dex=fit_subject_time_dex*pooled_glu4_GM_dex_slopes+pooled_glu4_GM_dex_intercept

pooled_glu4_WM_dex_slopes=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_WM_dex)[0]
pooled_glu4_WM_dex_intercept=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_WM_dex)[1]
pooled_glu4_WM_dex_rvalue=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_WM_dex)[2]
pooled_glu4_WM_dex_pvalue=stats.linregress(pooled_mrsi_subject_time_dex,pooled_glu4_WM_dex)[3]
pooled_glu4_WM_dex_r2value=pooled_glu4_WM_dex_rvalue*pooled_glu4_WM_dex_rvalue

pooled_linreg_plt_glu4_WM_dex=fit_subject_time_dex*pooled_glu4_WM_dex_slopes+pooled_glu4_WM_dex_intercept

##### 2H MRSI

pooled_glx_GM_mrsi_slopes=stats.linregress(pooled_subject_time_2H,pooled_glx_GM_mrsi)[0]
pooled_glx_GM_mrsi_intercept=stats.linregress(pooled_subject_time_2H,pooled_glx_GM_mrsi)[1]
pooled_glx_GM_mrsi_rvalue=stats.linregress(pooled_subject_time_2H,pooled_glx_GM_mrsi)[2]
pooled_glx_GM_mrsi_pvalue=stats.linregress(pooled_subject_time_2H,pooled_glx_GM_mrsi)[3]
pooled_glx_GM_mrsi_r2value=pooled_glx_GM_mrsi_rvalue*pooled_glx_GM_mrsi_rvalue

pooled_linreg_plt_glx_GM_mrsi=fit_subject_time_deu*pooled_glx_GM_mrsi_slopes+pooled_glx_GM_mrsi_intercept

pooled_glx_WM_mrsi_slopes=stats.linregress(pooled_subject_time_2H,pooled_glx_WM_mrsi)[0]
pooled_glx_WM_mrsi_intercept=stats.linregress(pooled_subject_time_2H,pooled_glx_WM_mrsi)[1]
pooled_glx_WM_mrsi_rvalue=stats.linregress(pooled_subject_time_2H,pooled_glx_WM_mrsi)[2]
pooled_glx_WM_mrsi_pvalue=stats.linregress(pooled_subject_time_2H,pooled_glx_WM_mrsi)[3]
pooled_glx_WM_mrsi_r2value=pooled_glx_WM_mrsi_rvalue*pooled_glx_WM_mrsi_rvalue

pooled_linreg_plt_glx_WM_mrsi=fit_subject_time_deu*pooled_glx_WM_mrsi_slopes+pooled_glx_WM_mrsi_intercept

pooled_glc_GM_mrsi_slopes=stats.linregress(pooled_subject_time_2H,pooled_glc_GM_mrsi)[0]
pooled_glc_GM_mrsi_intercept=stats.linregress(pooled_subject_time_2H,pooled_glc_GM_mrsi)[1]
pooled_glc_GM_mrsi_rvalue=stats.linregress(pooled_subject_time_2H,pooled_glc_GM_mrsi)[2]
pooled_glc_GM_mrsi_pvalue=stats.linregress(pooled_subject_time_2H,pooled_glc_GM_mrsi)[3]
pooled_glc_GM_mrsi_r2value=pooled_glc_GM_mrsi_rvalue*pooled_glc_GM_mrsi_rvalue

pooled_linreg_plt_glc_GM_mrsi=fit_subject_time_deu*pooled_glc_GM_mrsi_slopes+pooled_glc_GM_mrsi_intercept

pooled_glc_WM_mrsi_slopes=stats.linregress(pooled_subject_time_2H,pooled_glc_WM_mrsi)[0]
pooled_glc_WM_mrsi_intercept=stats.linregress(pooled_subject_time_2H,pooled_glc_WM_mrsi)[1]
pooled_glc_WM_mrsi_rvalue=stats.linregress(pooled_subject_time_2H,pooled_glc_WM_mrsi)[2]
pooled_glc_WM_mrsi_pvalue=stats.linregress(pooled_subject_time_2H,pooled_glc_WM_mrsi)[3]
pooled_glc_WM_mrsi_r2value=pooled_glc_WM_mrsi_rvalue*pooled_glc_WM_mrsi_rvalue

pooled_linreg_plt_glc_WM_mrsi=fit_subject_time_deu*pooled_glc_WM_mrsi_slopes+pooled_glc_WM_mrsi_intercept


########################## direct 2H MRSI


################# Linear Fitting glu4 mrsi DEU

glx_GM_deu_slopes=[stats.linregress(subject_time_2H[i],glx_GM_mrsi[i])[0] for i in np.arange(0,np.size(glx_GM_mrsi))]
glx_GM_deu_intercept=[stats.linregress(subject_time_2H[i],glx_GM_mrsi[i])[1] for i in np.arange(0,np.size(glx_GM_mrsi))]
glx_GM_deu_rvalue=[stats.linregress(subject_time_2H[i],glx_GM_mrsi[i])[2] for i in np.arange(0,np.size(glx_GM_mrsi))]
glx_GM_deu_pvalue=[stats.linregress(subject_time_2H[i],glx_GM_mrsi[i])[3] for i in np.arange(0,np.size(glx_GM_mrsi))]
glx_GM_deu_r2value=[glx_GM_deu_rvalue[i]*glx_GM_deu_rvalue[i] for i in np.arange(0,np.size(glx_GM_deu_rvalue))]

linreg_plt_glx_GM_deu=[fit_subject_time_deu*glx_GM_deu_slopes[i]+glx_GM_deu_intercept[i] for i in np.arange(0,np.size(glx_GM_deu_slopes))]


glx_WM_deu_slopes=[stats.linregress(subject_time_2H[i],glx_WM_mrsi[i])[0] for i in np.arange(0,np.size(glx_WM_mrsi))]
glx_WM_deu_intercept=[stats.linregress(subject_time_2H[i],glx_WM_mrsi[i])[1] for i in np.arange(0,np.size(glx_WM_mrsi))]
glx_WM_deu_rvalue=[stats.linregress(subject_time_2H[i],glx_WM_mrsi[i])[2] for i in np.arange(0,np.size(glx_WM_mrsi))]
glx_WM_deu_pvalue=[stats.linregress(subject_time_2H[i],glx_WM_mrsi[i])[3] for i in np.arange(0,np.size(glx_WM_mrsi))]
glx_WM_deu_r2value=[glx_WM_deu_rvalue[i]*glx_WM_deu_rvalue[i] for i in np.arange(0,np.size(glx_WM_deu_rvalue))]

linreg_plt_glx_WM_deu=[fit_subject_time_deu*glx_WM_deu_slopes[i]+glx_WM_deu_intercept[i] for i in np.arange(0,np.size(glx_WM_deu_slopes))]



################# Linear Fitting glu4 mrsi DEU

glc_GM_deu_slopes=[stats.linregress(subject_time_2H[i],glc_GM_mrsi[i])[0] for i in np.arange(0,np.size(glc_GM_mrsi))]
glc_GM_deu_intercept=[stats.linregress(subject_time_2H[i],glc_GM_mrsi[i])[1] for i in np.arange(0,np.size(glc_GM_mrsi))]
glc_GM_deu_rvalue=[stats.linregress(subject_time_2H[i],glc_GM_mrsi[i])[2] for i in np.arange(0,np.size(glc_GM_mrsi))]
glc_GM_deu_pvalue=[stats.linregress(subject_time_2H[i],glc_GM_mrsi[i])[3] for i in np.arange(0,np.size(glc_GM_mrsi))]
glc_GM_deu_r2value=[glc_GM_deu_rvalue[i]*glc_GM_deu_rvalue[i] for i in np.arange(0,np.size(glc_GM_deu_rvalue))]

linreg_plt_glc_GM_deu=[fit_subject_time_deu*glc_GM_deu_slopes[i]+glc_GM_deu_intercept[i] for i in np.arange(0,np.size(glc_GM_deu_slopes))]


glc_WM_deu_slopes=[stats.linregress(subject_time_2H[i],glc_WM_mrsi[i])[0] for i in np.arange(0,np.size(glc_WM_mrsi))]
glc_WM_deu_intercept=[stats.linregress(subject_time_2H[i],glc_WM_mrsi[i])[1] for i in np.arange(0,np.size(glc_WM_mrsi))]
glc_WM_deu_rvalue=[stats.linregress(subject_time_2H[i],glc_WM_mrsi[i])[2] for i in np.arange(0,np.size(glc_WM_mrsi))]
glc_WM_deu_pvalue=[stats.linregress(subject_time_2H[i],glc_WM_mrsi[i])[3] for i in np.arange(0,np.size(glc_WM_mrsi))]
glc_WM_deu_r2value=[glc_WM_deu_rvalue[i]*glc_WM_deu_rvalue[i] for i in np.arange(0,np.size(glc_WM_deu_rvalue))]

linreg_plt_glc_WM_deu=[fit_subject_time_deu*glc_WM_deu_slopes[i]+glc_WM_deu_intercept[i] for i in np.arange(0,np.size(glc_WM_deu_slopes))]






################# Exponential Fitting glu4 mrsi DEU

##Starting Values ### glu4/tcr ratios
m0=[0.5,0.5,0.5,0.5,0.5,0.5]  
tau=[20,20,20,20,20,20]
const=[0.4,0.4,0.4,0.4,0.4,0.4]

fitting_glu4_GM_deu=[fit(exp_decrease,mrsi_subject_time_deu[i],glu4_GM_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_GM_deu))]
m0_GM_deu=[fitting_glu4_GM_deu[i][0][0] for i in np.arange(0,np.size(glu4_GM_deu))]
tau_GM_deu=[fitting_glu4_GM_deu[i][0][1] for i in np.arange(0,np.size(glu4_GM_deu))]
const_GM_deu=[fitting_glu4_GM_deu[i][0][2] for i in np.arange(0,np.size(glu4_GM_deu))]

fitting_glu4_WM_deu=[fit(exp_decrease,mrsi_subject_time_deu[i],glu4_WM_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_WM_deu))]
m0_WM_deu=[fitting_glu4_WM_deu[i][0][0] for i in np.arange(0,np.size(glu4_WM_deu))]
tau_WM_deu=[fitting_glu4_WM_deu[i][0][1] for i in np.arange(0,np.size(glu4_WM_deu))]
const_WM_deu=[fitting_glu4_WM_deu[i][0][2] for i in np.arange(0,np.size(glu4_WM_deu))]

#plot_fit_glu4_GM=[exp_decrease(fit_subject_time_deu,m0_GM_deu[i],tau_GM_deu[i],const_GM_deu[i]) for i in np.arange(0,np.size(m0_GM_deu))]
#plot_fit_glu4_WM=[exp_decrease(fit_subject_time_deu,m0_WM_deu[i],tau_WM_deu[i],const_WM_deu[i]) for i in np.arange(0,np.size(m0_WM_deu))]

plot_fit_glu4_GM_deu=exp_decrease(fit_subject_time_deu,np.mean(m0_GM_deu),np.mean(tau_GM_deu),np.mean(const_GM_deu))
plot_fit_glu4_WM_deu=exp_decrease(fit_subject_time_deu,np.mean(m0_WM_deu),np.mean(tau_WM_deu),np.mean(const_WM_deu))


plt.figure('glu4_exp_deu')
plt.title('glu4_exp_deu')

[plt.plot(mrsi_subject_time_deu[i],glu4_GM_deu[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.size(mrsi_subject_time_deu))]

[plt.plot(mrsi_subject_time_deu[i],glu4_WM_deu[i],ls='none',marker='o',color='red') for i in np.arange(0,np.size(mrsi_subject_time_deu))]

#[plt.plot(fit_subject_time_deu,plot_fit_glu4_GM[x],color='blue') for x in np.arange(0,np.shape(plot_fit_glu4_GM)[0])]
#[plt.plot(fit_subject_time_deu,plot_fit_glu4_WM[x],color='red') for x in np.arange(0,np.shape(plot_fit_glu4_WM)[0])]

plt.plot(fit_subject_time_deu,plot_fit_glu4_GM_deu,color='blue',label='tau GM ='+str(np.round(np.mean(tau_GM_deu),2))+' min')
plt.plot(fit_subject_time_deu,plot_fit_glu4_WM_deu,color='red',label='tau WM ='+str(np.round(np.mean(tau_WM_deu),2))+' min')
plt.legend()

################# Exponential Fitting glu4 mrsi EXTRA

##Starting Values ### glu4/tcr ratios
m0=[0.5,0.5,0.5,0.5,0.5,0.5]  
tau=[20,20,20,20,20,20]
const=[0.4,0.4,0.4,0.4,0.4,0.4]

fitting_glu4_GM_extra=[fit(exp_decrease,mrsi_subject_time_extra[i],glu4_GM_extra[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
m0_GM_extra=[fitting_glu4_GM_extra[i][0][0] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
tau_GM_extra=[fitting_glu4_GM_extra[i][0][1] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]
const_GM_extra=[fitting_glu4_GM_extra[i][0][2] for i in np.arange(0,np.shape(glu4_GM_extra)[0])]

fitting_glu4_WM_extra=[fit(exp_decrease,mrsi_subject_time_extra[i],glu4_WM_extra[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
m0_WM_extra=[fitting_glu4_WM_extra[i][0][0] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
tau_WM_extra=[fitting_glu4_WM_extra[i][0][1] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]
const_WM_extra=[fitting_glu4_WM_extra[i][0][2] for i in np.arange(0,np.shape(glu4_WM_extra)[0])]

#plot_fit_glu4_GM=[exp_decrease(fit_subject_time_deu,m0_GM_deu[i],tau_GM_deu[i],const_GM_deu[i]) for i in np.arange(0,np.size(m0_GM_deu))]
#plot_fit_glu4_WM=[exp_decrease(fit_subject_time_deu,m0_WM_deu[i],tau_WM_deu[i],const_WM_deu[i]) for i in np.arange(0,np.size(m0_WM_deu))]

plot_fit_glu4_GM_extra=exp_decrease(fit_subject_time_extra,np.mean(m0_GM_extra),np.mean(tau_GM_extra),np.mean(const_GM_extra))
plot_fit_glu4_WM_extra=exp_decrease(fit_subject_time_extra,np.mean(m0_WM_extra),np.mean(tau_WM_extra),np.mean(const_WM_extra))


plt.figure('glu4_exp_extra')
plt.title('glu4_exp_extra')

[plt.plot(mrsi_subject_time_extra[i],glu4_GM_extra[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.shape(mrsi_subject_time_extra)[0])]

[plt.plot(mrsi_subject_time_extra[i],glu4_WM_extra[i],ls='none',marker='o',color='red') for i in np.arange(0,np.shape(mrsi_subject_time_extra)[0])]

#[plt.plot(fit_subject_time_deu,plot_fit_glu4_GM[x],color='blue') for x in np.arange(0,np.shape(plot_fit_glu4_GM)[0])]
#[plt.plot(fit_subject_time_deu,plot_fit_glu4_WM[x],color='red') for x in np.arange(0,np.shape(plot_fit_glu4_WM)[0])]

plt.plot(fit_subject_time_extra,plot_fit_glu4_GM_extra,color='blue',label='tau GM ='+str(np.round(np.mean(tau_GM_extra),2))+' min')
plt.plot(fit_subject_time_extra,plot_fit_glu4_WM_extra,color='red',label='tau WM ='+str(np.round(np.mean(tau_WM_extra),2))+' min')
plt.legend()



######## glu4_mrsi_linear deu pooled

plt.figure('glu4_mrsi_deu')
plt.title('glu4_mrsi_deu')

[plt.plot(mrsi_subject_time_deu[i],glu4_GM_deu[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.size(mrsi_subject_time_deu))]

[plt.plot(mrsi_subject_time_deu[i],glu4_WM_deu[i],ls='none',marker='o',color='red') for i in np.arange(0,np.size(mrsi_subject_time_deu))]

plt.plot(fit_subject_time_deu,pooled_linreg_plt_glu4_GM_deu,label='GM: R²='+str(np.round(np.mean(pooled_glu4_GM_deu_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_GM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_GM_deu_slopes),5)),color='blue')
plt.plot(fit_subject_time_deu,pooled_linreg_plt_glu4_WM_deu,label='WM: R²='+str(np.round(np.mean(pooled_glu4_WM_deu_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_WM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_WM_deu_slopes),5)),color='red')
plt.legend()

plt.ylim(0.3,0.7)
#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_deu[x],color='blue') for x in np.arange(0,np.size(linreg_plt_glu4_GM_deu))]
#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_deu[x],color='red') for x in np.arange(0,np.size(linreg_plt_glu4_WM_deu))]

#linreg_plt_glu4_GM_avg_deu=fit_subject_time_deu*np.mean(glu4_GM_deu_slopes)+np.mean(glu4_GM_deu_intercept)
#linreg_plt_glu4_WM_avg_deu=fit_subject_time_deu*np.mean(glu4_WM_deu_slopes)+np.mean(glu4_WM_deu_intercept)
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_avg_deu,color='blue',label='R²='+str(np.round(np.mean(glu4_GM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_GM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_GM_deu_slopes),5)))
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_avg_deu,color='red',label='R²='+str(np.round(np.mean(glu4_WM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_WM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_WM_deu_slopes),5)))


plt.savefig('../MRSI/glu4_mrsi_deu.pdf',format='pdf')



######## glu4_mrsi_linear EXTRA pooled

plt.figure('glu4_mrsi_extra')
plt.title('glu4_mrsi_extra')

[plt.plot(mrsi_subject_time_extra[i],glu4_GM_extra[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.shape(mrsi_subject_time_extra)[0])]

[plt.plot(mrsi_subject_time_extra[i],glu4_WM_extra[i],ls='none',marker='o',color='red') for i in np.arange(0,np.shape(mrsi_subject_time_extra)[0])]

plt.plot(fit_subject_time_extra,pooled_linreg_plt_glu4_GM_extra,label='GM: R²='+str(np.round(np.mean(pooled_glu4_GM_extra_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_GM_extra_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_GM_extra_slopes),5)),color='blue')
plt.plot(fit_subject_time_extra,pooled_linreg_plt_glu4_WM_extra,label='WM: R²='+str(np.round(np.mean(pooled_glu4_WM_extra_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_WM_extra_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_WM_extra_slopes),5)),color='red')
plt.legend()

plt.ylim(0.3,0.7)

#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_deu[x],color='blue') for x in np.arange(0,np.size(linreg_plt_glu4_GM_deu))]
#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_deu[x],color='red') for x in np.arange(0,np.size(linreg_plt_glu4_WM_deu))]

#linreg_plt_glu4_GM_avg_deu=fit_subject_time_deu*np.mean(glu4_GM_deu_slopes)+np.mean(glu4_GM_deu_intercept)
#linreg_plt_glu4_WM_avg_deu=fit_subject_time_deu*np.mean(glu4_WM_deu_slopes)+np.mean(glu4_WM_deu_intercept)
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_avg_deu,color='blue',label='R²='+str(np.round(np.mean(glu4_GM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_GM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_GM_deu_slopes),5)))
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_avg_deu,color='red',label='R²='+str(np.round(np.mean(glu4_WM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_WM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_WM_deu_slopes),5)))


plt.savefig('../MRSI/glu4_mrsi_extra.pdf',format='pdf')


######## glu4_mrsi_linear DEX pooled

plt.figure('glu4_mrsi_dex')
plt.title('glu4_mrsi_dex')

[plt.plot(mrsi_subject_time_dex[i],glu4_GM_dex[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.shape(mrsi_subject_time_dex)[0])]

[plt.plot(mrsi_subject_time_dex[i],glu4_WM_dex[i],ls='none',marker='o',color='red') for i in np.arange(0,np.shape(mrsi_subject_time_dex)[0])]

plt.plot(fit_subject_time_dex,pooled_linreg_plt_glu4_GM_dex,label='GM: R²='+str(np.round(np.mean(pooled_glu4_GM_dex_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_GM_dex_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_GM_dex_slopes),5)),color='blue')
plt.plot(fit_subject_time_dex,pooled_linreg_plt_glu4_WM_dex,label='WM: R²='+str(np.round(np.mean(pooled_glu4_WM_dex_r2value),2))+' p='+str(np.round(np.mean(pooled_glu4_WM_dex_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glu4_WM_dex_slopes),5)),color='red')
plt.legend()

plt.ylim(0.3,0.7)

#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_deu[x],color='blue') for x in np.arange(0,np.size(linreg_plt_glu4_GM_deu))]
#[plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_deu[x],color='red') for x in np.arange(0,np.size(linreg_plt_glu4_WM_deu))]

#linreg_plt_glu4_GM_avg_deu=fit_subject_time_deu*np.mean(glu4_GM_deu_slopes)+np.mean(glu4_GM_deu_intercept)
#linreg_plt_glu4_WM_avg_deu=fit_subject_time_deu*np.mean(glu4_WM_deu_slopes)+np.mean(glu4_WM_deu_intercept)
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_GM_avg_deu,color='blue',label='R²='+str(np.round(np.mean(glu4_GM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_GM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_GM_deu_slopes),5)))
#plt.plot(fit_subject_time_deu,linreg_plt_glu4_WM_avg_deu,color='red',label='R²='+str(np.round(np.mean(glu4_WM_deu_r2value),2))+' p='+str(np.round(np.mean(glu4_WM_deu_pvalue),2))+ ' slope='+str(np.round(np.mean(glu4_WM_deu_slopes),5)))


plt.savefig('../MRSI/glu4_mrsi_dex.pdf',format='pdf')

#### glx_2H_linear pooled 

plt.figure('glx_2H_mrsi')
plt.title('glx_2H_mrsi')

[plt.plot(subject_time_2H[i],glx_GM_mrsi[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.size(subject_time_2H))]
[plt.plot(subject_time_2H[i],glx_WM_mrsi[i],ls='none',marker='o',color='red') for i in np.arange(0,np.size(subject_time_2H))]

plt.plot(fit_subject_time_extra,pooled_linreg_plt_glx_GM_mrsi,label='GM: R²='+str(np.round(np.mean(pooled_glx_GM_mrsi_r2value),2))+' p='+str(np.round(np.mean(pooled_glx_GM_mrsi_pvalue),15))+ ' slope='+str(np.mean(pooled_glx_GM_mrsi_slopes))+'+-'+str(np.std(pooled_glx_GM_mrsi_slopes)),color='blue')
plt.plot(fit_subject_time_extra,pooled_linreg_plt_glx_WM_mrsi,label='WM: R²='+str(np.round(np.mean(pooled_glx_WM_mrsi_r2value),2))+' p='+str(np.round(np.mean(pooled_glx_WM_mrsi_pvalue),15))+ ' slope='+str(np.mean(pooled_glx_WM_mrsi_slopes))+'+-'+str(np.std(pooled_glx_WM_mrsi_slopes)),color='red')
plt.legend()



plt.savefig('../MRSI/glx_mrsi_2H.pdf',format='pdf')

#### glc_2H_linear

plt.figure('glc_2H_mrsi')
plt.title('glc_2H_mrsi')

[plt.plot(subject_time_2H[i],glc_GM_mrsi[i],ls='none',marker='o',color='blue') for i in np.arange(0,np.size(subject_time_2H))]
[plt.plot(subject_time_2H[i],glc_WM_mrsi[i],ls='none',marker='o',color='red') for i in np.arange(0,np.size(subject_time_2H))]

plt.plot(fit_subject_time_extra,pooled_linreg_plt_glc_GM_mrsi,label='GM: R²='+str(np.round(np.mean(pooled_glc_GM_mrsi_r2value),2))+' p='+str(np.round(np.mean(pooled_glc_GM_mrsi_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glc_GM_mrsi_slopes),5)),color='blue')
plt.plot(fit_subject_time_extra,pooled_linreg_plt_glc_WM_mrsi,label='WM: R²='+str(np.round(np.mean(pooled_glc_WM_mrsi_r2value),2))+' p='+str(np.round(np.mean(pooled_glc_WM_mrsi_pvalue),2))+ ' slope='+str(np.round(np.mean(pooled_glc_WM_mrsi_slopes),5)),color='red')
plt.legend()


plt.savefig('../MRSI/glc_mrsi_2H.pdf',format='pdf')



