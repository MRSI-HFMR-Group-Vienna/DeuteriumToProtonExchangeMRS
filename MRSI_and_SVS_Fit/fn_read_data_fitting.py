#!/usr/bin/python3


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

glu4_deu.append(np.array([7.591597578,7.061428165,6.651216276,6.837578857,6.552679738,6.529116653,6.204588709]))
glu4_deu.append(np.array([7.28942269,6.971726386,6.8020124,6.312532426,6.532953641,6.459479903]))
glu4_deu.append(np.array([7.499749229,7.135764737,6.880245538,6.993925508,6.725891083,6.537119756,6.571536628,6.489144723]))
glu4_deu.append(np.array([7.40554402,6.834823994,6.599946665,6.650960746,6.254538829]))
glu4_deu.append(np.array([7.361203632,7.151363164,6.722239406,6.328788527,6.34767417,6.398035882,5.997240587]))


glu4_dex=[]

glu4_dex.append(np.array([8.087082454,7.904463642,7.858245301,7.800754193,8.004791261,7.95180926,7.899954536]))
glu4_dex.append(np.array([7.470822804,7.720892186,7.559389043,7.711514585,7.678172]))
glu4_dex.append(np.array([7.997229466,8.042075693,8.155755663,8.071277887,8.10778063,8.238147568]))
glu4_dex.append(np.array([7.669753575,7.559502167,7.718638633,7.450290867,7.792486273,7.680154651]))
glu4_dex.append(np.array([7.704162352,7.778735703,7.672652485,7.9131778,7.864862671,7.838604449,7.831252147]))


gln4_deu=[]

gln4_deu.append(np.array([3.3855869,3.223858452,3.149956049,3.166021789,2.962522418,3.170305986,2.764378294]))
gln4_deu.append(np.array([2.392553274,2.246640639,2.110041577,2.279755563,2.145226184,2.147295866]))
gln4_deu.append(np.array([2.921262354,2.868072644,2.87120145,2.828441094,2.672000768,2.687644801,2.588565928,2.590651799]))
gln4_deu.append(np.array([2.577273859,2.569834306,2.338145357,2.545390059,2.117084341]))
gln4_deu.append(np.array([2.560053715,2.703794436,2.563201322,2.040698555,1.920040286,2.128831552,1.821415266]))


gln4_dex=[]

gln4_dex.append(np.array([3.295029553,3.088737932,3.009828569,3.083101549,3.110156188,3.118047124,3.075210613]))
gln4_dex.append(np.array([2.362113709,2.40379194,2.482980578,2.101624769,2.51423925]))
gln4_dex.append(np.array([3.009911873,2.993224904,3.058929841,2.854514482,3.0870891,2.886845483]))
gln4_dex.append(np.array([2.727162174,2.573226247,2.863416272,2.64603378,2.815571322,2.529541727]))
gln4_dex.append(np.array([2.818032391,2.816982062,2.700395556,2.98818567,2.85689456,2.743459041,2.823284036]))


gaba_deu=[]


gaba_deu.append(np.array([1.619426571,1.343095847,1.310964367,1.439490286,1.379511524,1.300253874,1.167443759]))
gaba_deu.append(np.array([1.253192913,0.572267284,1.047259478,0.833047312,1.200416002,1.167301078]))
gaba_deu.append(np.array([1.503913001,1.312012868,1.342257997,1.444465677,1.180602994,1.127413283,1.095082282,1.019990926]))
gaba_deu.append(np.array([1.098928318,1.030909544,0.96289077,1.134000498,1.03941189]))
gaba_deu.append(np.array([0.829919053,1.185598647,1.188746254,0.524601171,0.893920395,1.060743568,0.76591771]))

gaba_dex=[]


gaba_dex.append(np.array([1.427132198,1.158840363,1.23324062,1.349350111,1.416986708,1.348222835,1.572550881]))
gaba_dex.append(np.array([1.305570568,1.052375318,1.133647868,1.185745656,1.239927355]))
gaba_dex.append(np.array([1.431950451,1.29741177,1.308884061,1.436122193,1.565446196,1.51329942]))
gaba_dex.append(np.array([0.902813412,1.047388371,1.102514075,0.980821484,1.068190524,1.027626327]))
gaba_dex.append(np.array([1.051379211,0.987309149,1.074486447,1.120700917,1.072385789,0.927440403,0.99571178]))

#### Time Points
subject_time_deu=[]

subject_time_deu.append(np.array([21,39,49,60,71,81,93]))
subject_time_deu.append(np.array([24,42,55,66,78,90]))
subject_time_deu.append(np.array([24,42,54,66,77,87,98,109]))
subject_time_deu.append(np.array([20,43,55,67,78]))
subject_time_deu.append(np.array([21,38,49,60,71,82,97]))


subject_time_dex=[]

subject_time_dex.append(np.array([23,42,53,64,75,86,97]))
subject_time_dex.append(np.array([19,37,48,60,70]))
subject_time_dex.append(np.array([24,44,54,66,77,88]))
subject_time_dex.append(np.array([22,39,50,61,72,83]))
subject_time_dex.append(np.array([19,37,48,59,72,83,95]))



###################### DATA POINTS MRSI for 5 Subjects glu4/tcr gln4/tcr ratios for GM/WM

glu4_GM_mrsi_deu=[]
glu4_GM_mrsi_dex=[]

#################### glu4/tcr numbers GM after given deu and dex (control)

glu4_GM_mrsi_deu.append(np.array([0.552410157253131,0.512011808473785,0.491329914427083,0.484089421862397,0.478732139924454,0.479995729168255,0.460140380235207]))
glu4_GM_mrsi_deu.append(np.array([0.513499789886539,0.499711087035031,0.470204227236244,0.481128879855948,0.474579497473524,0.472186150984452]))
glu4_GM_mrsi_deu.append(np.array([0.51779996506957,0.485622577877703,0.479469399172592,0.462508937082936,0.466406635421806,0.44481780167264,0.44618496427095,0.442449856531019]))
glu4_GM_mrsi_deu.append(np.array([0.510889450300549,0.473980534501069,0.465039907477476,0.455054811205847,0.437326586063629,0.437633592080099,0.437884631164355,0.434889700906753]))
glu4_GM_mrsi_deu.append(np.array([0.539985772297478,0.52708562264263,0.522318909447299,0.514248782540732,0.502581523772467,0.498119378577269,0.495225531787315,0.490649397809115,0.489530441510711,0.482381751394495,0.471052552597617,0.468255401095264,0.460649084183253,0.455208176762363]))

glu4_GM_mrsi_dex.append(np.array([0.530463090032873,0.521588323000863,0.51603144842364,0.51603144842364,0.528055450773293,0.530169842714972,0.526752345659544]))
glu4_GM_mrsi_dex.append(np.array([0.57980395783244,0.57270126802055,0.577864616517295,0.568044397763079,0.582190396442699]))
glu4_GM_mrsi_dex.append(np.array([0.582839414028772,0.643462793790043,0.578633701633008,0.637334705588822,0.609523669681663,0.580487358844237]))
glu4_GM_mrsi_dex.append(np.array([0.503259178266538,0.5077753448348,0.504539595756434,0.510097347354328,0.510204081632653,0.508957643440068,0.5072,0.510359609430817]))
glu4_GM_mrsi_dex.append(np.array([0.522980723881392,0.526959968531574,0.533286325244876,0.530381891908249,0.532018636905407,0.525660894031265,0.532814456159404,0.536276433061987]))



glu4_WM_mrsi_deu=[]
glu4_WM_mrsi_dex=[]


#################### glu4/tcr numbers WM after given deu and dex (control)
glu4_WM_mrsi_deu.append(np.array([0.51310270832656,0.483450514673496,0.463658781920244,0.452808603134653,0.450231819213436,0.450044252138853,0.441160683788487]))
glu4_WM_mrsi_deu.append(np.array([0.464895635673624,0.44534833575097,0.432814984294412,0.442445488866631,0.426016106110848,0.421797958949662]))
glu4_WM_mrsi_deu.append(np.array([0.51779996506957,0.485622577877703,0.479469399172592,0.462508937082936,0.466406635421806,0.44481780167264,0.44618496427095,0.442449856531019]))
glu4_WM_mrsi_deu.append(np.array([0.470976898210203,0.445782943054398,0.431234241341439,0.419211136890951,0.416184971098266,0.409884653578795,0.410188351028733,0.40027131039974]))
glu4_WM_mrsi_deu.append(np.array([0.505997695067119,0.488940666841691,0.485848512321767,0.477209625173531,0.467413357966448,0.459606056633456,0.463672477756985,0.454533317757789,0.46193819736522,0.446770682750943,0.432512007368906,0.429864629183623,0.426903185426755,0.418728644294416,]))



glu4_WM_mrsi_dex.append(np.array([0.531618922514694,0.52457831721057,0.521564510215645,0.521564510215645,0.527966241932815,0.529112297350732,0.525121769999174]))
glu4_WM_mrsi_dex.append(np.array([0.500055098153425,0.497994884587451,0.501205687085056,0.495978510163687,0.495902029924712]))
glu4_WM_mrsi_dex.append(np.array([0.525081753817824,0.580730515875953,0.518578609149329,0.555951754962534,0.54910541098335,0.563272727272727]))
glu4_WM_mrsi_dex.append(np.array([0.472649693318901,0.469542157398019,0.479521806508745,0.467002721425037,0.476009434840632,0.477035959195678,0.479811097992916,0.479506800086744]))
glu4_WM_mrsi_dex.append(np.array([0.481573394108213,0.485949049519287,0.483391342547224,0.49136857690884,0.486995866573679,0.485919602514297,0.489636567751881,0.507786892206314]))


### time points

subject_time_mrsi_deu=[]
subject_time_mrsi_dex=[]


#### time points
subject_time_mrsi_deu.append(np.array([17,34,45,56,67,77,89]))
subject_time_mrsi_deu.append(np.array([20,38,51,63,74,87]))
subject_time_mrsi_deu.append(np.array([19,37,50,61,73,84,94,105]))
subject_time_mrsi_deu.append(np.array([16,34,45,56,67,78,91,107]))
subject_time_mrsi_deu.append(np.array([5,10,15,20,25,30,35,40,45,50,55,60,73,78]))


subject_time_mrsi_dex.append(np.array([19,38,49,60,71,82,93]))
subject_time_mrsi_dex.append(np.array([15,33,44,55,66]))
subject_time_mrsi_dex.append(np.array([21,39,50,61,72,84]))
subject_time_mrsi_dex.append(np.array([17,35,47,57,68,79,90,101]))
subject_time_mrsi_dex.append(np.array([15,32,44,55,68,79,91,103]))

################# Exponential Fitting glu4 DEU

##Starting Values ### glu4/tcr ratios
m0=[0.5,0.5,0.5,0.5,0.5,0.5]  
tau=[20,20,20,20,20,20]
const=[0,0,0,0,0,0]

### SVS
fitting_glu4=[fit(exp_decrease,subject_time_deu[i],glu4_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_deu))]
m0=[fitting_glu4[i][0][0] for i in np.arange(0,np.size(glu4_deu))]
tau=[fitting_glu4[i][0][1] for i in np.arange(0,np.size(glu4_deu))]
const=[fitting_glu4[i][0][2] for i in np.arange(0,np.size(glu4_deu))]

fit_subject_time_deu=np.arange(0,125,1)
fit_subject_time_dex=np.arange(0,125,1)
plot_fit_glu4=[exp_decrease(fit_subject_time_deu,m0[i],tau[i],const[i]) for i in np.arange(0,np.size(m0))]

#### MRSI
alldata_GM_glu4=[]
alldata_WM_glu4=[]
time_points_all_mrsi=[]
for y in np.arange(0,5):
	[alldata_GM_glu4.append(glu4_GM_mrsi_dex[y][x]) for x in np.arange(0,np.size(glu4_GM_mrsi_dex[y]))]
	[alldata_WM_glu4.append(glu4_WM_mrsi_dex[y][x]) for x in np.arange(0,np.size(glu4_WM_mrsi_dex[y]))]
	[time_points_all_mrsi.append(subject_time_mrsi_dex[y][x]) for x in np.arange(0,np.size(subject_time_mrsi_dex[y]))]

##Exponential
fitting_GM_mrsi_glu4=[fit(exp_decrease,subject_time_mrsi_deu[i],glu4_GM_mrsi_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_GM_mrsi_deu))]
fitting_WM_mrsi_glu4=[fit(exp_decrease,subject_time_mrsi_deu[i],glu4_WM_mrsi_deu[i],[m0[i],tau[i],const[i]]) for i in np.arange(0,np.size(glu4_WM_mrsi_deu))]

m0_GM=[fitting_GM_mrsi_glu4[i][0][0] for i in np.arange(0,np.size(glu4_GM_mrsi_deu))]
tau_GM=[fitting_GM_mrsi_glu4[i][0][1] for i in np.arange(0,np.size(glu4_GM_mrsi_deu))]
tau_GM_sd=[fitting_GM_mrsi_glu4[i][1][1] for i in np.arange(0,np.size(glu4_GM_mrsi_deu))]

const_GM=[fitting_GM_mrsi_glu4[i][0][2] for i in np.arange(0,np.size(glu4_GM_mrsi_deu))]

fit_subject_time_mrsi_deu=np.arange(0,125,1)
plot_fit_mrsi_GM_glu4=[exp_decrease(fit_subject_time_mrsi_deu,m0_GM[i],tau_GM[i],const_GM[i]) for i in np.arange(0,np.size(m0_GM))]

m0_WM=[fitting_WM_mrsi_glu4[i][0][0] for i in np.arange(0,np.size(glu4_WM_mrsi_deu))]
tau_WM=[fitting_WM_mrsi_glu4[i][0][1] for i in np.arange(0,np.size(glu4_WM_mrsi_deu))]
tau_WM_sd=[fitting_WM_mrsi_glu4[i][1][1] for i in np.arange(0,np.size(glu4_WM_mrsi_deu))]

const_WM=[fitting_WM_mrsi_glu4[i][0][2] for i in np.arange(0,np.size(glu4_WM_mrsi_deu))]

plot_fit_mrsi_WM_glu4=[exp_decrease(fit_subject_time_mrsi_deu,m0_WM[i],tau_WM[i],const_WM[i]) for i in np.arange(0,np.size(m0_WM))]

### ALL PLOTS MRSI intersubject, individual exponential fits
#### glu4_mrsi_deu exponential

plt.figure('glu4_mrsi_deu_averaged')
plt.title('glu4_mrsi_deu_averaged')

plot_fit_all_mrsi_WM=exp_decrease(fit_subject_time_mrsi_deu,np.average(m0_WM),np.average(tau_WM),np.average(const_WM))
plot_fit_all_mrsi_GM=exp_decrease(fit_subject_time_mrsi_deu,np.average(m0_GM),np.average(tau_GM),np.average(const_GM))

plt.plot(fit_subject_time_mrsi_deu,plot_fit_all_mrsi_GM,color='blue',label='GM: tau = '+str(np.round(np.average(tau_GM),0)) + '±' + str(np.round(np.std(tau_GM),0)) +' min')
plt.plot(fit_subject_time_mrsi_deu,plot_fit_all_mrsi_WM,color='red',label='WM: tau = '+str(np.round(np.average(tau_WM),0)) + '±' + str(np.round(np.std(tau_WM),0)) +' min')

[plt.plot(subject_time_mrsi_deu[x],glu4_GM_mrsi_deu[x],ls='none',marker='o',color='blue') for x in np.arange(0,5)]
[plt.plot(subject_time_mrsi_deu[x],glu4_WM_mrsi_deu[x],ls='none',marker='o',color='red') for x in np.arange(0,5)]
plt.legend()

plt.ylim(0.39,0.57)
#plt.ylim(0.4,0.5)

plt.savefig('glu4_mrsi_deu_exp_averaged.pdf',format='pdf')




### Linear fitting

(linfit_all_GM_mrsi_glu4_slopes,linfit_all_GM_mrsi_glu4_intercept,linfit_all_GM_mrsi_glu4_rvalue,linfit_all_GM_mrsi_glu4_pvalue,linfit_all_GM_mrsi_glu4_stderr)=stats.linregress(time_points_all_mrsi,alldata_GM_glu4)
(linfit_all_WM_mrsi_glu4_slopes,linfit_all_WM_mrsi_glu4_intercept,linfit_all_WM_mrsi_glu4_rvalue,linfit_all_WM_mrsi_glu4_pvalue,linfit_all_WM_mrsi_glu4_stderr)=stats.linregress(time_points_all_mrsi,alldata_WM_glu4)

linfit_GM_mrsi_glu4_slopes=[stats.linregress(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x])[0] for x in np.arange(0,5)]
linfit_GM_mrsi_glu4_intercept=[stats.linregress(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x])[1] for x in np.arange(0,5)]
linfit_GM_mrsi_glu4_rvalue=[stats.linregress(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x])[2] for x in np.arange(0,5)]
linfit_GM_mrsi_glu4_pvalue=[stats.linregress(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x])[3] for x in np.arange(0,5)]
linfit_GM_mrsi_glu4_stderr=[stats.linregress(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x])[4] for x in np.arange(0,5)]

linfit_WM_mrsi_glu4_slopes=[stats.linregress(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x])[0] for x in np.arange(0,5)]
linfit_WM_mrsi_glu4_intercept=[stats.linregress(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x])[1] for x in np.arange(0,5)]
linfit_WM_mrsi_glu4_rvalue=[stats.linregress(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x])[2] for x in np.arange(0,5)]
linfit_WM_mrsi_glu4_pvalue=[stats.linregress(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x])[3] for x in np.arange(0,5)]
linfit_WM_mrsi_glu4_stderr=[stats.linregress(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x])[4] for x in np.arange(0,5)]


plt.figure('glu4_mrsi_dex_linear_pooled')
plt.title('glu4_mrsi_dex_linear_pooled')

linfit_all_plot_GM_mrsi_glu4_dex=linfit_all_GM_mrsi_glu4_slopes*fit_subject_time_mrsi_deu+linfit_all_GM_mrsi_glu4_intercept
linfit_all_plot_WM_mrsi_glu4_dex=linfit_all_WM_mrsi_glu4_slopes*fit_subject_time_mrsi_deu+linfit_all_WM_mrsi_glu4_intercept

plt.plot(fit_subject_time_mrsi_deu,linfit_all_plot_GM_mrsi_glu4_dex,c='blue',label='GM: R² = ' + str(round(linfit_all_GM_mrsi_glu4_rvalue**2,2))+' p = ' + str(round(linfit_all_GM_mrsi_glu4_pvalue,2)))
plt.plot(fit_subject_time_mrsi_deu,linfit_all_plot_WM_mrsi_glu4_dex,c='red',label='WM: R² = ' + str(round(linfit_all_WM_mrsi_glu4_rvalue**2,2))+' p = ' + str(round(linfit_all_WM_mrsi_glu4_pvalue,2)))
[plt.plot(subject_time_mrsi_dex[x],glu4_GM_mrsi_dex[x],c='blue',ls='none',marker='o') for x in np.arange(0,5)]
[plt.plot(subject_time_mrsi_dex[x],glu4_WM_mrsi_dex[x],c='red',ls='none',marker='o') for x in np.arange(0,5)]
plt.legend()
plt.ylim(0.35,0.7)
plt.savefig('glu4_mrsi_dex_linear_pooled.pdf',format='pdf')



########################################## Linear Fitting SVS 
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



########################################################### ALL PLOTS SVS
#### glu4_deu exponential
plt.figure('glu4_deu')
plt.title('glu4_deu')

[plt.plot(subject_time_deu[i],glu4_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' tau='+str(np.round(tau[i],2))+' min') for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,plot_fit_glu4[i]) for i in np.arange(0,np.size(m0))]
plt.ylim(5,8.5)
plt.text(50,7.5,'Intersubject COV [tau] = '+str(np.round(np.std(tau)/np.mean(tau)*100,2))+'%')
plt.legend()

plt.savefig('glu4_svs_deu.pdf',format='pdf')

#### glu4_dex_linear
plt.figure('glu4_dex')
plt.title('glu4_dex')

[plt.plot(subject_time_dex[i],glu4_dex[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(np.round(glu4_dex_r2value[i],2))+' p='+str(np.round(glu4_dex_pvalue[i],2))) for i in np.arange(0,np.size(subject_time_dex))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,linreg_plt_glu4_dex[i]) for i in np.arange(0,np.size(glu4_dex_r2value))]
plt.legend()
plt.ylim(5,8.5)

plt.savefig('glu4_svs_dex.pdf',format='pdf')

#### gln4_deu_linear
plt.figure('gln4_deu')
plt.title('gln4_deu')

[plt.plot(subject_time_deu[i],gln4_deu[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(gln4_deu_r2value[i])+' p='+str(gln4_deu_pvalue[i])) for i in np.arange(0,np.size(subject_time_deu))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_deu,linreg_plt_gln4_deu[i]) for i in np.arange(0,np.size(gln4_deu_r2value))]
plt.legend()
plt.ylim(1.2,4)

plt.savefig('gln4_svs_deu.pdf',format='pdf')


#### gln4_dex_linear

plt.figure('gln4_dex')
plt.title('gln4_dex')

[plt.plot(subject_time_dex[i],gln4_dex[i],ls='none',marker='o',label='sub '+str(i+1)+' R²='+str(gln4_dex_r2value[i])+' p='+str(gln4_dex_pvalue[i])) for i in np.arange(0,np.size(subject_time_dex))]
plt.gca().set_prop_cycle(None) 
[plt.plot(fit_subject_time_dex,linreg_plt_gln4_dex[i]) for i in np.arange(0,np.size(gln4_dex_r2value))]
plt.legend()
plt.ylim(1.2,4)

plt.savefig('gln4_svs_dex.pdf',format='pdf')


####### read in 4D MRSI data of one representative subject with 14 time points
data_mrsi_L2=[]
g=open('MRSI_input/Glu_GD_amp_mapBLOCK1.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK2.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK3.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK4.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK5.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK6.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK7.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK8.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK9.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK10.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK11.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK12.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK13.txt')
data_mrsi_L2.append(g.readlines())
g=open('MRSI_input/Glu_GD_amp_mapBLOCK14.txt')
data_mrsi_L2.append(g.readlines())
g.close()

### GM Mask
f=open('MRSI_input/GM_mask.txt')
GM_mask=f.readlines()
f.close()
##### Magnitude
f=open('MRSI_input/magnitude_rs.txt')
magnitude_rs=f.readlines()
f.close()

#### WM Mask
f=open('MRSI_input/WM_mask.txt')
WM_mask=f.readlines()
f.close()


####################### Format Data
time_points_number=np.arange(0,14)
voxel_number=np.arange(0,np.size(data_mrsi_L2[0]))
for t in time_points_number:
	for v in voxel_number:
		data_mrsi_L2[t][v]=np.array(data_mrsi_L2[t][v].split())
	

mrsi_4D_L2=np.zeros(shape=(36,36,23,14))

for t in time_points_number:
	for v in voxel_number:
		mrsi_4D_L2[int(data_mrsi_L2[t][v][0])][int(data_mrsi_L2[t][v][1])][int(data_mrsi_L2[t][v][2])][t]=float(data_mrsi_L2[t][v][3])

############# Create Masks	
for v in voxel_number:
	GM_mask[v]=np.array(GM_mask[v].split())
	WM_mask[v]=np.array(WM_mask[v].split())
	magnitude_rs[v]=np.array(magnitude_rs[v].split())

mrsi_GM_mask=np.zeros(shape=(36,36,23))
mrsi_WM_mask=np.zeros(shape=(36,36,23))
mrsi_magitude_rs_mask=np.zeros(shape=(36,36,23))


for v in voxel_number:
	mrsi_GM_mask[int(GM_mask[v][0])][int(GM_mask[v][1])][int(GM_mask[v][2])]=float(GM_mask[v][3])
	mrsi_WM_mask[int(WM_mask[v][0])][int(WM_mask[v][1])][int(WM_mask[v][2])]=float(WM_mask[v][3])
	mrsi_magitude_rs_mask[int(magnitude_rs[v][0])][int(magnitude_rs[v][1])][int(magnitude_rs[v][2])]=float(magnitude_rs[v][3])
	
	
np.save('mrsi_GM_mask.npy',mrsi_GM_mask)
np.save('mrsi_WM_mask.npy',mrsi_WM_mask)
np.save('mrsi_magitude_rs_mask.npy',mrsi_magitude_rs_mask)


mrsi_time_points=np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 73, 78],dtype=float)

################################### Exponential fit

m0=100
tau=5
const=1

mrsi_3D_tau_L2=np.zeros(shape=(36,36,23))
for x in np.arange(0,36):
	for y in np.arange(0,36):
		for z in np.arange(0,23):
			mrsi_3D_tau_L2[x][y][z]=fit(exp_decrease,mrsi_time_points,mrsi_4D_L2[x][y][z],[m0,tau,const])[0][1]



mask_L2=mrsi_3D_tau_L2<200
masked_mrsi_3D_tau_L2=mrsi_3D_tau_L2*mask_L2

np.save('masked_mrsi_3D_tau_L2.npy',masked_mrsi_3D_tau_L2)
np.save('mrsi_3D_tau_L2.npy',mrsi_3D_tau_L2)



#### SD
mrsi_3D_tau_sd_L2=np.zeros(shape=(36,36,23))
for x in np.arange(0,36):
	for y in np.arange(0,36):
		for z in np.arange(0,23):
			mrsi_3D_tau_sd_L2[x][y][z]=fit(exp_decrease,mrsi_time_points,mrsi_4D_L2[x][y][z],[m0,tau,const])[1][1]


np.save('mrsi_3D_tau_sd_L2.npy',mrsi_3D_tau_sd_L2)


############ Linear Fitting

mrsi_3D_slopes_L2=np.zeros(shape=(36,36,23))
mrsi_3D_intercept_L2=np.zeros(shape=(36,36,23))
mrsi_3D_rvalue_L2=np.zeros(shape=(36,36,23))
mrsi_3D_pvalue_L2=np.zeros(shape=(36,36,23))
mrsi_3D_stderr_L2=np.zeros(shape=(36,36,23))


for x in np.arange(0,36):
	for y in np.arange(0,36):
		for z in np.arange(0,23):
				(mrsi_3D_slopes_L2[x][y][z],mrsi_3D_intercept_L2[x][y][z],mrsi_3D_rvalue_L2[x][y][z],mrsi_3D_pvalue_L2[x][y][z],mrsi_3D_stderr_L2[x][y][z])=stats.linregress(mrsi_time_points,mrsi_4D_L2[x][y][z])
	


np.save('mrsi_3D_slopes_L2_map.npy',mrsi_3D_slopes_L2)
np.save('mrsi_3D_intercept_L2_map.npy',mrsi_3D_intercept_L2)
np.save('mrsi_3D_rvalue_L2_map.npy',mrsi_3D_rvalue_L2)
np.save('mrsi_3D_pvalue_L2_map.npy',mrsi_3D_pvalue_L2)
np.save('mrsi_3D_stderr_L2_map.npy',mrsi_3D_stderr_L2)




############################################### PLOT GM AND WM
#GM
slopes_GM_abs_L2=mrsi_3D_slopes_L2[(mrsi_3D_pvalue_L2<0.05)*(mrsi_3D_rvalue_L2<-0.8)*(mrsi_3D_slopes_L2*mrsi_GM_mask<0)]
#WM
slopes_WM_abs_L2=mrsi_3D_slopes_L2[(mrsi_3D_pvalue_L2<0.05)*(mrsi_3D_rvalue_L2<-0.8)*(mrsi_3D_slopes_L2*mrsi_WM_mask<0)]








