import tkinter as tk

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
import matplotlib as mpl
#import matplotlib.pyplot as mpl_plot
#from matplotlib.figure import Figure
##from matplotlib.image import NonUniformImage

import numpy as np
import re


####################################
# Transport coefficient definitions
####################################


version='sims1'

#if (version == "duke"):
if (re.match("duke(-[a-z]+)?",version) is not None):

    shear_param_list={
    "eta_over_s_at_min":[0.,.3],
    "high_T_slope_in_GeV":[.0,2.],
    "high_T_curvature":[-1,1]
    }

    def eta_over_s( T_in_GeV, parameters):

        T_min_in_GeV=.154
        eta_over_s_at_min=parameters['eta_over_s_at_min']
        high_T_slope_in_GeV=parameters['high_T_slope_in_GeV']
        high_T_curvature=parameters['high_T_curvature']

        res=eta_over_s_at_min+high_T_slope_in_GeV*(T_in_GeV - T_min_in_GeV)*np.power(T_in_GeV/T_min_in_GeV,high_T_curvature)
            
        return res


    if (re.match("duke-jf",version) is not None):

        bulk_param_list={
        "zetas_max":[0.001, 0.3],
        "zetas_area_fourth":[0.0002**.25, (0.2*0.3)**.25],
        "zetas_t0":[.17,.2]
        }

        def zeta_over_s( T_in_GeV, parameters):

            zetas_max=np.maximum(parameters['zetas_max'],bulk_param_list['zetas_max'][0])
            zetas_area_fourth=parameters['zetas_area_fourth']
            zetas_t0=parameters['zetas_t0']
            zetas_width=2./3.141592*(zetas_area_fourth**4/zetas_max)

            res=zetas_max/(1+np.power((T_in_GeV - zetas_t0)/zetas_width,2))

            return res
    else:
        bulk_param_list={
        "T_peak_in_GeV":[0.15,.2],
        "zeta_over_s_at_peak":[.0,0.1],
        "width_in_GeV":[1e-5,.1]
        }

        def zeta_over_s( T_in_GeV, parameters):

            T_peak_in_GeV=parameters['T_peak_in_GeV']
            zeta_over_s_at_peak=parameters['zeta_over_s_at_peak']
            width_in_GeV=parameters['width_in_GeV']

            res=zeta_over_s_at_peak/(1+np.power((T_in_GeV - T_peak_in_GeV)/width_in_GeV,2))

            return res


if (version == "sims1"):

    shear_param_list={
    "T_kink_in_GeV":[.15,.22],
    "eta_over_s_at_kink":[.001,.25],
    "low_T_slope_in_GeV":[-1,2],
    "high_T_slope_in_GeV":[-1,1]
    }

    def eta_over_s( T_in_GeV, parameters):

        T_kink_in_GeV=parameters['T_kink_in_GeV']
        eta_over_s_at_kink=parameters['eta_over_s_at_kink']
        low_T_slope_in_GeV=parameters['low_T_slope_in_GeV']
        high_T_slope_in_GeV=parameters['high_T_slope_in_GeV']

        if (T_in_GeV<T_kink_in_GeV):
            res= eta_over_s_at_kink + low_T_slope_in_GeV*(T_kink_in_GeV - T_in_GeV)
        else:
            res=eta_over_s_at_kink + high_T_slope_in_GeV*(T_in_GeV - T_kink_in_GeV)

        res=np.maximum(0.001,res)

        return res

    bulk_param_list={
    "T_peak_in_GeV":[0.13,.3],
    "zeta_over_s_at_peak":[1e-5,0.3],
    "zeta_area_root":[(1e-5)**.25,.5**.25],
    "lambda":[-.8,0.8]
    }

    #\frac{\zeta}{s}(T)=\frac{1}{\pi\sigma_\zeta}\frac{A_{\zeta}}{1+\left(T-T_{\zeta,c}\right)^2\left[\sigma^2_{\zeta}\left(\lambda sign(T-T_{\zeta,c})+1\right)\right]^{-1}}
    #$T_{\zeta,c}\in [0.13,0.5]$ GeV; Jonah had this $T_{\zeta,c}\in [0.15,0.2]$ GeV
    #$\lambda\in[-1,1]$; Johah had this at $\lambda=0$
    #$\sigma_{\zeta}\in[0.03,0.15]$ GeV; Jonah had this $\sigma_{\zeta}\in[0.03,0.1]$ GeV
    #$A_{\zeta}\in[0,1]$; Jonah had this $A_{\zeta}\in[0,0.1]$


    def zeta_over_s( T_in_GeV, parameters):

        T_peak_in_GeV=parameters['T_peak_in_GeV']
        zeta_over_s_at_peak=np.maximum(parameters['zeta_over_s_at_peak'],0.001)
        zeta_over_s_area=(parameters['zeta_area_root'])**4.
        lambda_param=parameters['lambda']
        zeta_over_s_width=1./3.141592*(zeta_over_s_area/zeta_over_s_at_peak)

        res=zeta_over_s_at_peak/(1+np.power((T_in_GeV - T_peak_in_GeV)/zeta_over_s_width/(lambda_param*np.sign(T_in_GeV-T_peak_in_GeV)+1),2))

        return res


####################################
# Widget to play with parameters
####################################

root = tk.Tk()
root.wm_title("eta/s and zeta/s")

# Canvas for plots
fig = mpl.figure.Figure(figsize=(8, 4), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
#canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
canvas.get_tk_widget().grid(row=0,columnspan=2)

# Need to track the id of plotted lines so that they can be deleted
line_list=[None,None]

# List of axes
ax_list=[None,None]

# List of sliders
sliders_list=[{},{}]

# Function to update the plots when the sliders are varied
def eta_slider_update(val):
    slider_update("shear")

def zeta_slider_update(val):
    slider_update("bulk")

#
def slider_update(viscosity):

    # 
    if (viscosity == 'shear'):
        tid=0
        lparam_list=shear_param_list
        transport_fct=eta_over_s
    else:
        tid=1
        lparam_list=bulk_param_list
        transport_fct=zeta_over_s

    # Figure out which line to delete
    if (line_list[tid]):
        line_list[tid].remove()

    parameters={ key: float(sliders_list[tid][key].get()) for key, values in lparam_list.items() }

    T_range = np.arange(0, .6, .001)
    yaxis_res=[transport_fct(T,parameters) for T in T_range]

    lines, = ax_list[tid].plot(T_range, yaxis_res, 'b')
    line_list[tid]=lines
    canvas.draw()

# Loop over the two transport coefficients
transport_coeffs=['eta','zeta']
for transport_coeff in transport_coeffs:

    transport_id=0 if transport_coeff == 'eta' else 1
    tmp_ax=fig.add_subplot(1,2,transport_id+1)
    ax_list[transport_id] =tmp_ax

    axis_label=r'$\eta/s$' if transport_coeff == 'eta' else r'$\zeta/s$'
    tmp_ax.set_title(axis_label)

    tmp_ax.set_xlim(0.1,.6)
    if (transport_coeff == 'eta'):
        y_low_lim=0.0
        y_high_lim=0.5
    else:
        y_low_lim=0.0
        if (re.match("^duke-jf",version) is not None):
            y_high_lim=0.4
        if (re.match("^sims",version) is not None):
            y_high_lim=0.4
        else:
            y_high_lim=0.15

    tmp_ax.set_ylim(y_low_lim,y_high_lim)


    if (transport_coeff == 'eta'):
        param_list=shear_param_list
        slider_update_fct=eta_slider_update
        transport_fct=eta_over_s
    else:
        param_list=bulk_param_list
        slider_update_fct=zeta_slider_update
        transport_fct=zeta_over_s

  
    ###########
    # Sliders
    ###########

    # Create all the sliders
    for n, (param_name, (param_min, param_max)) in enumerate(param_list.items()):
        tmp_scale = tk.Scale(orient='horizontal', from_=param_min, to=param_max, resolution=(param_max-param_min)/50., command=slider_update_fct, label=param_name, length=200)
        sliders_list[transport_id][param_name]=tmp_scale
        tmp_scale.grid(row=n+1,column=transport_id)

    # Set sliders to mean value
    for param_name, (param_min, param_max) in param_list.items():
        sliders_list[transport_id][param_name].set((param_min+param_max)/2.0)



    ###################################################
    # Show prior given limits on the parameters
    ###################################################

    # Here's the idea: histogram \eta/s(T) and plot its density

    # temperature bins
    T_bins_edges = np.linspace(0.1, 0.6, 20)
    T_bins= (T_bins_edges[:-1] + T_bins_edges[1:]) / 2
    T_bins_size=len(T_bins)

    # transport bins
    transport_bins_edges=np.linspace(y_low_lim,y_high_lim,20)
    transport_bins= (transport_bins_edges[:-1] + transport_bins_edges[1:]) / 2
    transport_bins_size=len(transport_bins)

    # sample transport's parameter space
    samples_per_param=15
    meshgrid_input={param_name: np.linspace(param_min,param_max,samples_per_param) for (param_name, (param_min, param_max)) in param_list.items()}

    param_name_list=[ param_name for (param_name, (param_min, param_max)) in param_list.items() ]
    number_of_params=len(param_name_list)

    pre=[ meshgrid_input[name] for name  in param_name_list]

    param_sample=np.meshgrid(*pre)

    number_of_samples=np.power(samples_per_param,number_of_params)

    histo=np.zeros([T_bins_size,transport_bins_size])
    unit=1./number_of_samples
    #transport_for_histo=np.empty([T_array_size,number_of_samples])
    index_tuple=tuple(samples_per_param for j in range(number_of_params))
    for Ti, T in enumerate(T_bins): #range(T_bins_size):
        #tmp_array=np.empty([number_of_samples])
        for si in range(number_of_samples):
           
            # Convert sample iterator into a N-d index
            indices=np.unravel_index(si, index_tuple)
            tmp_dict=meshgrid_input={param_name: param_sample[n][indices] for n, (param_name, (param_min, param_max)) in enumerate(param_list.items())}
            tmp_transport=transport_fct(T,tmp_dict)

            transport_index=np.digitize(tmp_transport,transport_bins_edges)-1

            #if (T>.4):
            #    print(tmp_transport,transport_bins_edges,transport_index)

            # If the transport falls outsize of our binning, ignore it
            if (transport_index <0)or(transport_index>=transport_bins_size):
                continue

            histo[Ti][transport_index]+=unit

    cmap=mpl.cm.get_cmap("Purples")
    norm=mpl.colors.Normalize(vmin=1e-50, vmax=0.15)
    cmap.set_under(color='white')
    # This produces a binned 2D plot
    #tmp_ax.pcolormesh(T_bins_edges, transport_bins_edges, histo.T, cmap=cmap, vmin=1e-50, vmax=.25)
    # This produces a smooth 2D plot
    im = mpl.image.NonUniformImage(tmp_ax, interpolation='bilinear', cmap=cmap)
    im.set_norm(norm)
    im.set_data(T_bins, transport_bins, histo.T)
    tmp_ax.images.append(im)



####
# For validation 
####
#
#eta_res=np.loadtxt("eta_over_s.dat")
#zeta_res=np.loadtxt("zeta_over_s.dat")
#
#ax_list[0].plot(eta_res[:,0], eta_res[:,1], 'r')
#ax_list[1].plot(zeta_res[:,0], zeta_res[:,1], 'r')
#
#canvas.draw()


tk.mainloop()
