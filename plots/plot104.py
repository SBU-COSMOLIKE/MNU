import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

# GENERAL PLOT OPTIONS
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'

chaindir=os.getcwd()

analysissettings={'smooth_scale_1D':0.35,'smooth_scale_2D':0.35,'ignore_rows': u'0.5',
'range_confidence' : u'0.01'}

analysissettings2={'smooth_scale_1D':0.35,'smooth_scale_2D':0.35,'ignore_rows': u'0.0',
'range_confidence' : u'0.01'}

root_chains = (
  'EXAMPLE_MCMC37',
  'EXAMPLE_MCMC38',
  'EXAMPLE_MCMC39',
  'EXAMPLE_MCMC40',
  'EXAMPLE_MCMC41',
  'EXAMPLE_MCMC42',
  'EXAMPLE_MCMC43',
  'EXAMPLE_MCMC44',
  'EXAMPLE_MCMC45'
)

num_points_thin = 50000

# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[0],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP1')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[1],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP2')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[2],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP3')
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[3],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP4')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[4],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP5')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[5],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP6')
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[6],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP7')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[7],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP8')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + '/../chains/' + root_chains[8],settings=analysissettings)
samples.thin(factor = int(np.sum(samples.weights)/num_points_thin))
p = samples.getParams()
samples.addDerived(p.omegam*p.H0/100.,name='gamma',label='{\\Omega_m h}')
samples.addDerived(p.s8omegamp5/0.5477225575,name='SS8',label='{S_8}')
samples.addDerived(p.mnu,name='mnu2',label='{\\sum m_\\nu}', range=[0.06,0.7])
samples.saveAsText(chaindir + '/.VM_P104_TMP9')
# --------------------------------------------------------------------------------



#GET DIST PLOT SETUP
g = gplot.get_single_plotter(
  chain_dir=chaindir,
  analysis_settings=analysissettings2,
  width_inch=5.3
)

g.settings.axis_tick_x_rotation=65
g.settings.lw_contour = 1.2
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 13.0
g.settings.legend_fontsize = 9.0
g.settings.alpha_filled_add = 0.85
g.settings.lab_fontsize=15.5
g.legend_labels=False
g.settings.constrained_layout=True

g.plots_1d(
  [ 
    chaindir + '/.VM_P104_TMP1', 
    chaindir + '/.VM_P104_TMP4', 
    chaindir + '/.VM_P104_TMP7', 
    chaindir + '/.VM_P104_TMP3',
    chaindir + '/.VM_P104_TMP6',
    chaindir + '/.VM_P104_TMP9'
  ],
  params=[u'mnu2'],
  line_args=[
    {'lw': 1.0,'ls': 'solid', 'color':'lightcoral'},
    {'lw': 1.2,'ls': '--', 'color':'lightcoral'},
    {'lw': 1.4,'ls': '-.', 'color': 'lightcoral'},
    {'lw': 1.7,'ls': 'solid', 'color':'indigo'},
    {'lw': 1.9,'ls': '--', 'color':'indigo'},
    {'lw': 2.2,'ls': '-.', 'color': 'indigo'}
  ]
)

g.settings.tight_layout=True
g.finish_plot(
  legend_labels=[
    '$\\Lambda \\quad \\,\\,\\,$  P1 + ACTDR4 TTTEEE',
    '$\\Lambda \\quad \\,\\,\\,$  P2 + ACTDR4 TTTEEE',
    '$\\Lambda \\quad \\,\\,\\,$  P3 + ACTDR4 TTTEEE',
    '$w_0w_a$ P1 + SPT TTTEEE',
    '$w_0w_a$ P2 + SPT TTTEEE',
    '$w_0w_a$ P3 + SPT TTTEEE'
  ],
  no_extra_legend_space=True,
  legend_loc=(0.58,0.65)
)

ax = g.get_axes()
ax.set_xlim(0.06,0.7)

# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP3',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP3', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='indigo',alpha=0.25) 
# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP6',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP6', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='indigo',alpha=0.25) 
# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP9',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP9', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='indigo',alpha=0.25) 
# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP1',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP1', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='lightcoral',alpha=0.5) 
# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP4',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP4', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='lightcoral',alpha=0.5) 
# --------------------------------------------------------------------------------
samples        = loadMCSamples(chaindir + '/.VM_P104_TMP7',settings=analysissettings2)
stats          = samples.getMargeStats()
sampleAnalyser = gplot.MCSampleAnalysis(chaindir,settings=analysissettings2)
lims           = stats.parWithName('mnu2').limits
density        = sampleAnalyser.get_density(chaindir + '/.VM_P104_TMP7', param='mnu2')

tmp0 = []
tmp1 = []
for i in range(0,len(density.x)) :
  if((density.x[i] < lims[1].lower) or density.x[i] > lims[1].upper) :
    tmp0.append(density.x[i])
    tmp1.append(density.P[i])
plt.gca().fill_between(tmp0, 0, tmp1, facecolor='lightcoral',alpha=0.5) 
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------

g.export()

