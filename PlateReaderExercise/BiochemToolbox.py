"""A library of code to support Chem 3530

This library contains a short collection of functions for use in *Python* 
notebooks. Import the library using `import BiochemToolbox as BT`. You can 
then access this documentation with the command help(BT).

help(BT) will also display the documentation for all the functions in 
this library.
"""

from scipy.stats import pearsonr
from scipy.optimize import curve_fit      ## import tools
import numpy as np                       
from matplotlib import pyplot as plt     
import pandas as pd
from sklearn.metrics import r2_score

import uncertainties as un               ## import Uncertainties tool
from uncertainties import unumpy as unp
from uncertainties import umath as um

def MM(S, Vmax, KM):
    """Michaelis-Menten calculation of initial rate from substrate conc.
    
    Arguments:
    ----------
    S: float or numpy array
        The substrate concentration(s)
    Vmax : float
        The Michaelis-Menten Vmax value
    KM : float
        The Michaelis-Menten Km value

    Returns:
    -------
    float or numpy array
        The calculated initial rate(s)
    """

    return Vmax * S / (S + KM)

def MM_curve_fit(x,y):
    """A function to fit x,y data to the Michaelis-Menten equation
    
    Arguments:
    ----------
    x, y : numpy arrays
        The x,y data for the Michaelis-menten experiment.
        x is the concentration of substrate
        y is the initial rate
        
    Returns:
    -------
    Vmax, KM : Uncertainty objects
        Returns a pair of uncertain numbers using Uncertainty objects.
        The first is Vmax and the second is KM
    """
    ### Perform the curve fit
    params, stats = curve_fit(MM, x, y)   ## two objects are returned

    ### Interpret the results
    v_max, KM = params   ### pull out the two values in the params object

    perr = np.sqrt(np.diag(stats)) ### convert covariance matrix to stdev
    stdev_v_max, stdev_KM = perr   ### pull out the two stdev values

    v_max_u = un.ufloat(v_max,stdev_v_max)  ### create uncertainty values
    KM_u = un.ufloat(KM,stdev_KM)

    return v_max_u, KM_u

def MM_Plot(x, y, x_lim = (None,None), y_lim = (None,None),
            title = "", x_label = "", y_label = "", file_name = "Plot.pdf" ):
    r"""A function to make a MM plot and return the MM parameters

    This function will take an x,y set of values for substrate concentration 
    and a set of style parameters and output a plot and save the plot as a pdf 
    file. It will return the curve fit parameters for Vmax and Km as
    Uncertainties objects (value with st_dev)

    Arguments:
    x (numpy array):  The values for concentration of substrate
    y (numpy array): The values for observed initial rate of reaction
    x_lim (tuple): (optional) The min and max values for the x axis. 
    y_lim (tuple): (optional) The min and max values for the y axis. 
    title (str): (optional) A string for adding a title to the plot. 
    x_label (str): (optional) A string for adding a label to the x-axis.
        Can be a raw string with LaTex math or a simple string.
        e.g. "Substrate" or r"[S] / $10^{-3} mole L^{-1}$" 
    y_label (str): (optional) A string for adding a y-axis label
        Can be a raw string with LaTex math or a simple string
        e.g. "Initial rate" or r"$\nu / 10^{-3}mole L^{-1}$" 
    file_name (str): (optional) The file name for the output plot.
    
    Returns:
    -------
    Vmax: Uncertainties, KM : Uncertainties 
    This function outputs a michaelis-menten plot with the line fit
    displayed on the data. It returnes two uncertainties objects, Vmax and KM.

    Vmax.n will give the nominal value
    Vmax.s will give the standard error
    """

    ####################################
    ### Perform the curve fit
    ####################################
    
    v_max, KM = MM_curve_fit(x,y)

    ################################
    ### make a list of x values from zero to the end of the line
    ################################
    
    x_fit = np.linspace(0, np.max(x), 100)
    
    ################################
    ### Feed that list into the function for the line fit
    ################################
    
    y_fit = MM(x_fit, v_max.n, KM.n)
    
    ######################
    ### Create an empty plot
    #####################
        
    plt.rcdefaults()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
    
    ######################
    ### Plot the data and the curve fit line
    #####################

    ax.scatter(x, y,
            marker='o',                 ### Plot the x and y data 
            color='white',              ### markers are this color
            edgecolors = 'black',       ### outline of markers is this color
            linewidths = 0.5,           ### outline of markers is this wide
            s=36,                       ### "s" is "size". sqrt(64) = 8 points wide
            zorder = 1                  ### everything is in layer 1 (the top layer in this case)
            )

    ax.plot(x_fit, y_fit,
            linestyle = '-',            ### use a line between points
            linewidth='0.5',            ### make the line thin
            color = 'black',            ### the line is black
            zorder = 0                  ### everything is in layer 0 (the bottom layer in this case)
            )

                                                                                                     
    ######################
    ### Apply the style parameters
    #####################

    ax.set(xlabel=x_label,
            ylabel=y_label,
            title = title,
            xlim=x_lim,
            ylim=y_lim
          )

    ######################
    ### Display and export the plot
    #####################
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(file_name)
    plt.show()

    return v_max, KM

def linear(x, slope, intercept):  ### Take x values, slope and intercept and return the y values
    """Linear line calculation from slope and intercept.
    
    Arguments:
    x: float or numpy array
        The x value or array of values
    slope (float)
    intercept (float)

    Returns:
    float or numpy array
        The calculated y values
    """
    return slope * x + intercept

def line_curve_fit(x,y):
    """A function to fit x,y data to a linear equation
    
    Arguments:
    ----------
    x, y : numpy arrays
        The x,y data .
        
    Returns:
    -------
    slope, intercept : Uncertainty objects
        Returns a pair of uncertain numbers using Uncertainty objects.
        The first is slope and the second is intercept
    """

    ####################################
    ### Perform the curve fit
    ####################################

    params, stats = curve_fit(linear, x, y)   ## two objects are returned
    
    ####################################
    ### Interpret the results
    ####################################
    
    slope, intercept = params   ### pull out the two values in the params object
    
    perr = np.sqrt(np.diag(stats))        ### convert covariance matrix to stdev values
    stdev_slope, stdev_intercept = perr   ### pull out the two stdev values
    
    slope_u = un.ufloat(slope,stdev_slope)               ## create a value with uncertainty built in
    intercept_u = un.ufloat(intercept,stdev_intercept)
    
    return slope_u, intercept_u

def Linear_Plot(x, y, x_lim = (None,None), y_lim = (None,None),           
            title = "", x_label = "", y_label = "", file_name = "Plot.pdf" ):
    r"""A function to make a linear plot and return the fit parameters

    This function will take an x,y set of values and a set of style 
    parameters and output a plot and save the plot as a pdf file. It will 
    return the linear fit parameters for slope and intercept as Uncertainties
    objects (value with st_dev)

    Arguments:
    x (numpy array):  The values for concentration of substrate
    y (numpy array): The values for observed initial rate of reaction
    x_lim (tuple): (optional) The min and max values for the x axis. 
    y_lim (tuple): (optional) The min and max values for the y axis. 
    title (str): (optional) A string for adding a title to the plot. 
    x_label (str): (optional) A string for adding a label to the x-axis.
        Can be a raw string with LaTex math or a simple string.
        e.g. "1/Substrate" or r"1/([S] / $10^{-3} mole L^{-1}$)" 
    y_label (str): (optional) A string for adding a y-axis label
        Can be a raw string with LaTex math or a simple string
        e.g. "1/Initial rate" or r"1/($\nu$ / $10^{-3} mole L^{-1}$)" 
    file_name (str): (optional) The file name for the output plot.
    
    Returns:
    slope (Uncertainties), intercept (Uncertainties): This function outputs 
        a linear plot with the line fit displayed on the data and returns 
        two uncertainties objects, slope and intercept. 
        slope.n will give the nominal value. 
        slope.s will give the standard error. 
    """
        
    ####################################
    ### Perform the curve fit
    ####################################
    
    slope, intercept = line_curve_fit(x,y)

    ################################
    ### make a list of x values from zero to the end of the line
    ################################
    
    x_fit = np.linspace(0, np.max(x), 100)  ## 100 points from 0 to the highest value on LB plot x-axis
    
    ################################
    ### Feed that list into the function for the line fit
    ################################
    
    y_fit = linear(x_fit, slope.n, intercept.n)
    
    ######################
    ### Create an empty plot
    #####################
        
    plt.rcdefaults()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
    
    ######################
    ### Plot the data and the curve fit line
    #####################

    ax.scatter(x, y,
            marker='o',                 ### Plot the x and y data 
            color='white',              ### markers are this color
            edgecolors = 'black',       ### outline of markers is this color
            linewidths = 0.5,           ### outline of markers is this wide
            s=64,                       ### "s" is "size". sqrt(64) = 8 points wide
            zorder = 1                  ### everything is in layer 1 (the top layer in this case)
            )

    ax.plot(x_fit, y_fit, 
            linestyle = '-',            ### use a line between points
            linewidth='0.5',            ### make the line thin
            color = 'black',            ### the line is black
            zorder = 0                  ### everything is in layer 0 (the bottom layer in this case)
            )
                                                                                     
    ######################
    ### Apply the style parameters
    #####################

    ax.set(xlabel=x_label,
            ylabel=y_label,
            title = title,
            xlim=x_lim,
            ylim=y_lim
          )

    ######################
    ### Display and export the plot
    #####################
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(file_name)
    plt.show()

    return slope, intercept

def get_integrated_MM_function():
    """Returns the integrated Michaelis-Meneten rate law and a function
    
    No parameters
    
    returns (eq, f):
        eq (SymPy equation object): 
            The integrated MM rate law
        f (function object): f(t, S0, KM, Vmax)
            t (float or array): time point or array of points
            S0 (float): Initial substrate concentration
            KM (float): KM value for enzyme
            Vmax (float): Vmax for enzyme
            returns s (float or array): substrate conc at time = t
    
    Example: 
    MMeq, s_conc = get_integrated_MM_function(t, S0, KM, Vmax)
    display(MMeq) # show the integrated rate equation
    s = s_conc(t_list, S0_value, KM_value, Vmax_value) # s at each t in t_list

    """
    import sympy as sym

    ##########################
    ### Set up equation 
    ##########################

    t = sym.symbols('t')           ### create x as a 'symbol', not a variable
    Vmax = sym.symbols('V_{max}')  ### create k as a 'symbol'
    St = sym.symbols('S_t')        ### create At as a 'symbol'
    S0 = sym.symbols('S_0')        ### create A0 as a 'symbol'
    KM = sym.symbols('K_M')

    xt = sym.Function('x_t')       ### create x as a 'function', not a variable

    lhs = sym.Derivative(-St, t)   ### Using Derivative function to get differential of A(t) w.r.t. t
    rhs = Vmax*(St/(KM+St))

    diffeq = sym.Eq(lhs, rhs)                   ### create a sympy equation
    diffeq = diffeq.subs({St: (S0 - xt(t))})    ### Substitute the term with S0-x

    ##########################
    ### Solve the differential equation
    ##########################

    res = sym.dsolve(diffeq, ics={xt(0): 0})    ### Solve the differential equation. Initial condition is x(t) = 0 when t = 0

    ##########################
    ### Clean up algebra 
    ##########################

    eq = res.subs(xt(t), S0-St)            ### substitute x for So - St
    eq = sym.simplify(eq)                  ### Simplify the result
    eq = sym.Eq(eq.lhs - S0, eq.rhs - S0)  ### Subtract S0 from both sides of the equation
    eq = sym.Eq(-eq.lhs, -eq.rhs)          ### take the negative of both sides of the equation

    ##########################
    ### Display the final form of equation
    ##########################

    #print("The integrated rate law for the MM equation")
#    display(eq)                         

    ##########################
    ### create function in terms of t, S0, KM and Vmax
    ##########################

    f = sym.lambdify([t, S0, KM, Vmax], eq.rhs)
    return eq,f

def calculate_e_NPA(pH = 7.0):
    ### parameters to get extinction coeff for NPA at give pH value
    e_NPAA = 18300  ### extinction coeff for NPA anion
    pKa_NPA = 7.15  ### pKa for p-nitrophenol

    ### Calculated Values from the above lists
    Ka = 10 ** -pKa_NPA   ### extinction coeff for NPA at given pH
    H = 10 ** -pH
    e_NPA = e_NPAA * (Ka / (H + Ka))   # abs of phenolate anion at pH value
    
    return e_NPA

def read_plate_setup(file_name):
    """Reads a file name for a plate plan and returns lists of the data

    Arguments:
    ----------
    file_name: string

    Returns:
    --------
    pandas dataframe
        The dataframe read in from the file. Some columns are different 
        lengths so it will be separated into lists (Pandas series) and 
        those are returned in the dictionary below.

    Example
    -------
    df = read_plate_setup{"platesetup.csv"}
    """

    df = pd.read_csv(file_name,
                     comment = "#",
                     skipinitialspace=True)
    return df

def plot_wells(data_file_path, plate_name_list, Column_list, Row_list,
               Fraction_time_span = 1,
               Make_Plots = False, 
               plot_file_path = "plots/analysis_plot_",
               Save_Data = False, 
               result_file_path = "data1/data/analysis_results",
               fancy = False, tiny_points = False, tiny_line = False):

    """Loads and plots data files. Can return line fit.

    Will make a plot of all the wells designated by columns and rows
    Will fit data to linear curve fit and collect slopes and intercept data
    Will save the plot as a pdf
    Will save the data as a csv file with same filename as pdf
    
    Arguments
    ---------
    data_file_name: string
        The data file name. Filenames will be this plus column and row
        e.g. "name_10_C.csv"
    Column_list, Row_list: Array like
        list of names for column, row to be plotted.
        Will usually be one column and some rows but can be as many as wanted
    Fraction_time_span: float
        The fraction of the time span to plt. Default is 1 for 100%
    Line_Fit: boolean
        If True then line fits will be performed for each well and data
        written to lists and put into the result dataframe
    Display_Plot, Display_Data: booleans
        If True then the plot will be created and displayed and the
        dataframe displayed, respectively.
    
    Returns
    -------
    result: pandas dataframe
        The results of line fits. Will be empty if Line_Fit = False
    """
    def linear_function(x, slope, intercept):
        return slope * x + intercept
    
    #print(Column_list)

    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    

    slope_list = []; slope_stderr_list = []
    int_list = []; int_stderr_list = []; rsq_list = [];
    plate_list = []; well_lane_list = []; well_row_list=[]

    for plate_name in plate_name_list:
        #print(plate_name)
        for lane_name in Column_list:
        #print(lane_name)

            for row_name in Row_list:
                in_file_name = data_file_path + plate_name \
                                + "_" + str(lane_name) \
                                + "_" + row_name + ".csv"
                df = pd.read_csv(in_file_name)
                points_used = int(Fraction_time_span * len(df["time"]))

                x = df["time"][0:points_used]
                y = df["abs"][0:points_used]


                if True:
                    param,cov = curve_fit(linear_function, x, y)
                    slope, intercept = un.correlated_values(param,cov)


                    rsq = r2_score(y, linear_function(x, *param))

                    slope_list.append(slope.n)
                    slope_stderr_list.append(slope.s)
                    int_list.append(intercept.n)
                    int_stderr_list.append(intercept.s)
                    rsq_list.append(rsq)
                    well_lane_list.append(str(lane_name))
                    well_row_list.append(row_name)
                    plate_list.append(plate_name)
                ### end of if:

                if Make_Plots:
                    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))

                    x_fit = np.linspace(0,np.max(x),10)
                    ax.plot(x_fit, linear_function(x_fit, *param),
                            linestyle = '-',
                            linewidth='0.5',
                            color = 'black',
                            zorder = 0)

                    if tiny_line:
                        ax.plot(x, y,
                            #marker=None,
                            color='black',
                            linewidth = 0.5,
                            ms=8,
                            zorder = 1)
                    if tiny_points:
                        ax.scatter(x, y,
                            marker='o',
                            color='black',
                            edgecolors = 'none',
                            linewidths = 0.5,
                            s=2,
                            zorder = 2)
                    else:
                        ax.scatter(x, y,
                            marker='o',
                            color='white',
                            edgecolors = 'black',
                            linewidths = 0.5,
                            s=32,
                            zorder = 2)
                        ax.scatter(x, y,
                            marker='o',
                            color='white',
                            edgecolors = None,
                            linewidths = 0.5,
                            s=64,
                            zorder = 1)
                    ### end of if tiny_points:
                ### end of for row_name:
            ###end of for column_name:
                    if np.min(y) < 0:
                        q = np.min(y) - 0.05*np.max(y)
                    else:
                        q = -0.05*np.max(y)

                    ax.set(xlabel= r"Time $/min$",
                           ylabel=r"$A_{405}$",
                           title = f"{plate_name}-{lane_name}-{row_name}",
                           xlim=[0, None],
                           ylim=[q, None]
                            )
                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.savefig(plot_file_path +"_"+ plate_name +"_"+ str(lane_name) \
                                +"_"+ row_name + ".pdf")     ### export the plot as this
                    plt.close()


    results = {"Plate": plate_list,
               "Column":well_lane_list,
               "Row":well_row_list,
               "slope":slope_list,
               "slope stderr":slope_stderr_list,
               "int": int_list,
               "int stderr":int_stderr_list,
               "RSQ": rsq_list}
    results = pd.DataFrame(results)

        #display(results)

    if Save_Data:
        results.to_csv(result_file_path + ".csv")
        print("Data saved as " + result_file_path + ".csv")

    return results

def dual_plot_w_residuals(filename, lane_name, row_name,
                          Fraction_time_span = 1,
                          plot_file = "plots/Cell_w_residuals_2_plot",
                          fancy = False):
    
    """Plot abs vs time data for a well with linear fit and residuals
    
    Arguments
    ---------
    
    filename: string
        The file root name. the well column and row will be added to 
        this name to get the file name needed.
    lane_name, row_name: strings
        The column and row label. Will be used to access the file
    Fraction_time_span: float
        The fraction of the time span to plt. Default is 1 for 100%
    plot_file: string
        The file name of the plot to be written as pdf. Will append the
        column and row label to give "file_12_A.pdf" for example

    Returns:
    --------
    
    Nul

    Does not return any objects but will output a figure with Two plots 
        Plot 1: Plot with time span chose by Fraction_time_span
        Plot 2: residuals for plot 
    """
    
    def linear_function(x, slope, int):
        return slope*x + int
    
    def linear_function_int0(x, slope):
        return slope*x
    
    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7,3))
    
    in_file_name = filename + "_" + lane_name + "_" + row_name + ".csv"
    df = pd.read_csv(in_file_name)
    
    points_used = int(Fraction_time_span * len(df["time"]))
    
    x = df["time"][0:points_used]
    y = df["abs"][0:points_used]
    
    param,cov = curve_fit(linear_function, x,y)
    slope, intercept = param
    rsq = r2_score(y, linear_function(x, *param))
    
    ##param,cov = curve_fit(linear_function_int0, x,y)
    ##[slope] = param
    
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr, int_stderr = perr
    print(f"slope = {slope:0.3g} +/- {slope_stderr:0.3g}, R2 = {rsq:0.3f}")
    ##print(f"slope = {slope:0.3g}")
    x_fit = np.linspace(0,np.max(x),10)
    
    ax[0].plot(x_fit, linear_function(x_fit, slope, intercept),
               linestyle = '-',
               linewidth='0.5',
               color = 'black',
               zorder = 0)
    ax[0].scatter(x, y,
                  marker='o',
                  color='lightgray',
                  edgecolors = 'black',
                  linewidths = 0.5,
                  s=8,
                  zorder = 2)
#    ax[0].scatter(x, y,
#                  marker='o',
#                  color='white',
#                  edgecolors = None,
#                  linewidths = 0.5,
#                  s=32,
#                  zorder = 1)
    ax[0].set(xlabel= r"Time $/min$",
              ylabel=r"$A_{405}$",
     #              title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[-0.05*np.max(y), None]
             )
    
    residuals = y - linear_function(x, slope, intercept)
    y = residuals
    ax[1].hlines(0, xmin = 0, xmax = np.max(x),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[1].plot(x, y,
#               linestyle = '-',
#               linewidth='3',
#               color = 'white',
#               zorder = 1)
#    ax[1].plot(x, y,
#               linestyle = '-',
#               linewidth='0.5',
#               color = 'black',
#               zorder = 1)
    ax[1].scatter(x, y,
                  marker='o',
                  color='lightgray',
                  edgecolors = 'black',
                  linewidths = 0.5,
                  s=8,
                  zorder = 3)
#    ax[1].scatter(x, y,
#                  marker='o',
#                  color='white',
#                  edgecolors = None,
#                  linewidths = 0.5,
#                  s=32,
#                  zorder = 2)
    ax[1].set(xlabel= r"Time  $/min$",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.02, +0.02])
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(plot_file+"_"+lane_name+"_"+row_name+".pdf")     ### export the plot as this
    plt.show()

    print("Plot saved as "+plot_file+"_"+lane_name+"_"+row_name+".pdf")

def plot_w_residuals(filename, lane_name, row_name,
                          Fraction_time_span = 1,
                          plot_file = "plots/plot_w_residuals.pdf",
                          fancy = False):
    
    """Plot abs vs time data for a well with linear fit and residuals
    
    Arguments
    ---------
    
    filename: string
        The file root name. the well column and row will be added to 
        this name to get the file name needed.
    lane_name, row_name: strings
        The column and row label. Will be used to access the file
    Fraction_time_span: float
        The fraction of the time span to plt. Default is 1 for 100%
    plot_file: string
        The file name of the plot to be written as pdf. Will append the
        column and row label to give "file_12_A.pdf" for example

    Returns:
    --------
    
    Nul

    Does not return any objects but will output a figure with Two plots 
        Plot 1: Plot with time span chose by Fraction_time_span
        Plot 2: residuals for plot 
    """
    
    def linear_function(x, slope, int):
        return slope*x + int
    
    def linear_function_int0(x, slope):
        return slope*x
    
    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=[4,5], height_ratios=[1, 4])  
    
    in_file_name = filename + "_" + lane_name + "_" + row_name + ".csv"
    df = pd.read_csv(in_file_name)
    
    points_used = int(Fraction_time_span * len(df["time"]))
    
    x = df["time"][0:points_used]
    y = df["abs"][0:points_used]
    
    param,cov = curve_fit(linear_function, x,y)
    slope, intercept = param
    rsq = r2_score(y, linear_function(x, *param))
    
    ##param,cov = curve_fit(linear_function_int0, x,y)
    ##[slope] = param
    
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr, int_stderr = perr
    print(f"slope = {slope:0.3g} +/- {slope_stderr:0.3g}, R2 = {rsq:0.3f}")
    ##print(f"slope = {slope:0.3g}")
    x_fit = np.linspace(0,np.max(x),10)
    
    ax[1].plot(x_fit, linear_function(x_fit, slope, intercept),
               linestyle = '-',
               linewidth='0.5',
               color = 'black',
               zorder = 0)
    ax[1].scatter(x, y,
                  marker='o',
                  color='lightgray',
                  edgecolors = 'black',
                  linewidths = 0.5,
                  s=8,
                  zorder = 2)
#    ax[1].scatter(x, y,
#                  marker='o',
#                  color='white',
#                  edgecolors = None,
#                  linewidths = 0.5,
#                  s=32,
#                  zorder = 1)

    if np.min(y) < 0:
        q = 1.1 * np.min(y)
    else:
        q = -0.05*np.max(y)
    ax[1].set(xlabel= r"Time $/min$",
            ylabel=r"$A_{405}$",
     #             title = "Lane # "+lane_name,
            xlim=[-0.05*np.max(x), None],           
            ylim=[q, None]
            )
    
    residuals = y - linear_function(x, slope, intercept)
    y = residuals
    ax[0].hlines(0, xmin = 0, xmax = np.max(x),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[0].plot(x, y,
#               linestyle = '-',
#               linewidth='3',
#               color = 'white',
#               zorder = 1)
#    ax[0].plot(x, y,
#               linestyle = '-',
#               linewidth='0.5',
#               color = 'black',
#               zorder = 1)
    ax[0].scatter(x, y,
                  marker='o',
                  color='lightgray',
                  edgecolors = 'black',
                  linewidths = 0.5,
                  s=8,
                  zorder = 3)
#    ax[0].scatter(x, y,
#                  marker='o',
#                  color='white',
#                  edgecolors = None,
#                  linewidths = 0.5,
#                  s=32,
#                  zorder = 2)
    ax[0].set(xlabel= r"",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.01, +0.01]
               )
    ax[0].set_xticks([])


    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(plot_file+"_"+lane_name+"_"+row_name+".pdf")     ### export the plot as this
    plt.show()

    print("Plot saved as "+plot_file+"_"+lane_name+"_"+row_name+".pdf")

def plot_four_w_residuals(filename, lane_name, row_name,
                          Fraction_time_span_medium = 0.2,
                          Fraction_time_span_short = 0.05,
                          plot_file = "plots/Cell_w_residuals_4_plot",
                          fancy = False):

    """DEPRECATED - Use plot_six_w_residuals() instead
    Plot 2x2 plot grid. abs vs time data for a well with different time spans
    
    Arguments
    ---------
    
    filename: string
        The file root name. the well column and row will be added to 
        this name to get the file name needed.
    lane_name, row_name: strings
        The column and row label. Will be used to access the file
    Fraction_time_span_medium: float
        The fraction of the time span to plot is second plot. Default is 0.2 for 20%
    Fraction_time_span_short: float
        The fraction of the time span to plot is third plot. Default is 0.05 for 5%
    plot_file: string
        The file name of the plot to be written as pdf. Will append the
        column and row label to give "file_12_A.pdf" for example

    Returns:
    -------

    Nul

    Does not return any objects but will output a figure with four plots 
        Plot 1: full plot
        Plot 2: plot with time fraction medium
        Plot 3: Plot with time fraction short
        Plot 4: residuals fopr plot with short time span

    Will print slope and stderr for slope for each of plot 2 and plot 3.
    """

    def linear_function(x, slope, int):
        return slope*x + int
    
    def linear_function_int0(x, slope):
        return slope*x
    
    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7,6))
    
    in_file_name = filename + "_" + lane_name + "_" + row_name + ".csv"
    df = pd.read_csv(in_file_name)
    
    points_used = int(Fraction_time_span_short * len(df["time"]))
    points_used2 = int(Fraction_time_span_medium * len(df["time"]))
    
    x_all = df["time"]
    y_all = df["abs"]
    
    x = x_all[0:points_used]
    y = y_all[0:points_used]
    
    x2 = x_all[0:points_used2]
    y2 = y_all[0:points_used2]
    
    #########################################################################

    param, cov = curve_fit(linear_function, x_all, y_all)
    slope_all, intercept_all = param
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr_all, int_stderr_all = perr
    print(f"All slope = {slope_all:0.3g} +/- {slope_stderr_all:0.3g}")
    
    x_fit = np.linspace(0,np.max(x_all),10)
    
    ax[0][0].plot(x_fit, linear_function(x_fit, slope_all, intercept_all),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)

    ax[0][0].scatter(x_all, y_all,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[0][0].scatter(x_all, y_all,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 1)
    if np.min(y_all) < 0:
        q = np.min(y_all) - 0.05*np.max(y_all)
    else:
        q = -0.05*np.max(y_all)

    ax[0][0].set(xlabel= r"Long Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )
    
    #########################################################################
    
    param, cov = curve_fit(linear_function, x2, y2)
    slope2, intercept2 = param
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr2, int_stderr2 = perr
    print(f"Med slope = {slope2:0.3g} +/- {slope_stderr2:0.3g}")
    
    x_fit = np.linspace(0,np.max(x2),10)
    
    ax[0][1].plot(x_fit, linear_function(x_fit, slope2, intercept2),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)
    ax[0][1].scatter(x2, y2,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[0][1].scatter(x2, y2,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=16,
#            zorder = 1)

    if np.min(y2) < 0:
        q = np.min(y2) - 0.05*np.max(y2)
    else:
        q = -0.05*np.max(y2)

    ax[0][1].set(xlabel= r"Med Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )
    
    #########################################################################
    
    param,cov = curve_fit(linear_function, x,y)
    slope, intercept = param
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr, int_stderr = perr
    print(f"short slope = {slope:0.3g} +/- {slope_stderr:0.3g}")
    
    x_fit = np.linspace(0,np.max(x),10)
    
    ax[1][0].plot(x_fit, linear_function(x_fit, slope, intercept),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)
    ax[1][0].scatter(x, y,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[1][0].scatter(x, y,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 1)

    if np.min(y) < 0:
        q = np.min(y) - 0.05*np.max(y)
    else:
        q = -0.05*np.max(y)

    ax[1][0].set(xlabel= r"Short Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )
    
    #########################################################################
    
    #x=x2; y=y2; slope=slope2; intercept=intercept2
    
    residuals = y - linear_function(x, slope, intercept)
    y = residuals
    ax[1][1].hlines(0, xmin = 0, xmax = np.max(x),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[1][1].plot(x, y,
#            linestyle = '-',
#            linewidth='3',
#            color = 'white',
#            zorder = 1)
#    ax[1][1].plot(x, y,
#            linestyle = '-',
#            linewidth='0.5',
#            color = 'black',
#            zorder = 1)
    ax[1][1].scatter(x, y,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8, 
            zorder = 3)
#    ax[1][1].scatter(x, y,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 2)
    ax[1][1].set(xlabel= r"Short Time  $/min$",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.01, +0.01]
               )
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(plot_file+"_"+lane_name+"_"+row_name+".pdf")     ### export the plot as this
    plt.show()

    print("Plot saved as "+plot_file+"_"+lane_name+"_"+row_name+".pdf")

def plot_six_w_residuals(filename, lane_name, row_name,
                          Fraction_time_span_medium = 0.2,
                          Fraction_time_span_short = 0.05,
                          plot_file = "plots/Cell_w_residuals_6_plot",
                          fancy = False):

    """Plot 3x2 plot grid. The first column is abs vs time data 
    and the second column is differemtial plots. Each row is a different time
    span. Top row is full time span, the remain two rows are fractions of
    that span.
    
    Arguments
    ---------
    
    filename: string
        The file root name. the well column and row will be added to 
        this name to get the file name needed.
    lane_name, row_name: strings
        The column and row label. Will be used to access the file
    Fraction_time_span_medium: float
        The fraction of the time span to plot is second plot. Default is 0.2 for 20%
    Fraction_time_span_short: float
        The fraction of the time span to plot is third plot. Default is 0.05 for 5%
    plot_file: string
        The file name of the plot to be written as pdf. Will append the
        column and row label to give "file_12_A.pdf" for example

    Returns:
    -------

    Nul

    Does not return any objects but will output a figure with six plots 
        Plot 1: full plot
        Plot 2: residuals for full plot 
        Plot 2: plot with time fraction medium
        Plot 4: residuals for medium time span plot 
        Plot 5: Plot with time fraction short
        Plot 6: residuals for plot with short time span

    Will print slope and stderr for slope for each of abs vs time plots.
    """
#    from sklearn.metrics import r2_score

    def linear_function(x, slope, int):
        return slope*x + int
    
    def linear_function_int0(x, slope):
        return slope*x
    
    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.close()
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(7,9))
    
    in_file_name = filename + "_" + lane_name + "_" + row_name + ".csv"
    df = pd.read_csv(in_file_name)
    
    points_used = int(Fraction_time_span_short * len(df["time"]))
    points_used2 = int(Fraction_time_span_medium * len(df["time"]))
    
    x_all = df["time"]
    y_all = df["abs"]
    
    x = x_all[0:points_used]
    y = y_all[0:points_used]
    
    x2 = x_all[0:points_used2]
    y2 = y_all[0:points_used2]
    
    #########################################################################

    param, cov = curve_fit(linear_function, x_all, y_all)
    slope_all, intercept_all = param
    rsq = r2_score(y_all, linear_function(x_all, *param))

    perr = np.sqrt(np.diag(cov))
    slope_stderr_all, int_stderr_all = perr
    print(f"All slope = {slope_all:0.3g} +/- {slope_stderr_all:0.3g}, r2 = {rsq:0.3f}")
    
    x_fit = np.linspace(0,np.max(x_all),10)
    
    ax[0][0].plot(x_fit, linear_function(x_fit, slope_all, intercept_all),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)

    ax[0][0].scatter(x_all, y_all,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[0][0].scatter(x_all, y_all,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 1)
    if np.min(y_all) < 0:
        q = np.min(y_all) - 0.05*np.max(y_all)
    else:
        q = -0.05*np.max(y_all)

    ax[0][0].set(xlabel= r"Long Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )


    #########################################################################
    
    #x=x2; y=y2; slope=slope2; intercept=intercept2
    
    residuals = y_all - linear_function(x_all, slope_all, intercept_all)

    ax[0][1].hlines(0, xmin = 0, xmax = np.max(x_all),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[0][1].plot(x_all, residuals,
#            linestyle = '-',
#            linewidth='3',
#            color = 'white',
#            zorder = 1)
#    ax[0][1].plot(x_all, residuals,
#            linestyle = '-',
#            linewidth='0.5',
#            color = 'black',
#            zorder = 1)
    ax[0][1].scatter(x_all, residuals,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8, 
            zorder = 3)
#    ax[0][1].scatter(x_all, residuals,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 2)
    ax[0][1].set(xlabel= r"Long Time  $/min$",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.01, +0.01]
               )

    #########################################################################

    param, cov = curve_fit(linear_function, x2, y2)
    slope2, intercept2 = param
    rsq = r2_score(y2, linear_function(x2, *param))

    perr = np.sqrt(np.diag(cov))
    slope_stderr2, int_stderr2 = perr
    print(f"Med slope = {slope2:0.3g} +/- {slope_stderr2:0.3g}, r2 = {rsq:0.3f}")

    x_fit = np.linspace(0,np.max(x2),10)

    ax[1][0].plot(x_fit, linear_function(x_fit, slope2, intercept2),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)
    ax[1][0].scatter(x2, y2,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[1][0].scatter(x2, y2,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=16,
#            zorder = 1)

    if np.min(y2) < 0:
        q = np.min(y2) - 0.05*np.max(y2)
    else:
        q = -0.05*np.max(y2)

    ax[1][0].set(xlabel= r"Med Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )

    #########################################################################
    
    #x=x2; y=y2; slope=slope2; intercept=intercept2
    
    residuals = y2 - linear_function(x2, slope2, intercept2)
    

    ax[1][1].hlines(0, xmin = 0, xmax = np.max(x2),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[1][1].plot(x2, residuals,
#            linestyle = '-',
#            linewidth='3',
#            color = 'white',
#            zorder = 1)
#    ax[1][1].plot(x2, residuals,
#            linestyle = '-',
#            linewidth='0.5',
#            color = 'black',
#            zorder = 1)
    ax[1][1].scatter(x2, residuals,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8, 
            zorder = 3)
#    ax[1][1].scatter(x2, residuals,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 2)
    ax[1][1].set(xlabel= r"Medium Time  $/min$",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.01, +0.01]
               )
    #########################################################################

    param,cov = curve_fit(linear_function, x,y)
    slope, intercept = param
    rsq = r2_score(y, linear_function(x, *param))
    
    perr = np.sqrt(np.diag(cov))
    slope_stderr, int_stderr = perr
    print(f"short slope = {slope:0.3g} +/- {slope_stderr:0.3g}, r2 = {rsq:0.3f}")
    
    x_fit = np.linspace(0,np.max(x),10)
    
    ax[2][0].plot(x_fit, linear_function(x_fit, slope, intercept),
            linestyle = '-',
            linewidth='0.5',
            color = 'black',
            zorder = 0)
    ax[2][0].scatter(x, y,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8,
            zorder = 2)
#    ax[2][0].scatter(x, y,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 1)

    if np.min(y) < 0:
        q = np.min(y) - 0.05*np.max(y)
    else:
        q = -0.05*np.max(y)

    ax[2][0].set(xlabel= r"Short Time $/min$",
              ylabel=r"$A_{405}$",
     #         title = "Lane # "+lane_name,
              xlim=[-0.05*np.max(x), None],
              ylim=[q, None]
             )
    
    #########################################################################
    
    #x=x2; y=y2; slope=slope2; intercept=intercept2
    
    residuals = y - linear_function(x, slope, intercept)

    ax[2][1].hlines(0, xmin = 0, xmax = np.max(x),
                 colors='black', linestyles='solid',
                 linewidths = 0.5, zorder = 0)
#    ax[2][1].plot(x, residuals,
#            linestyle = '-',
#            linewidth='3',
#            color = 'white',
#            zorder = 1)
#    ax[2][1].plot(x, residuals,
#            linestyle = '-',
#            linewidth='0.5',
#            color = 'black',
#            zorder = 1)
    ax[2][1].scatter(x, residuals,
            marker='o',
            color='lightgray',
            edgecolors = 'black',
            linewidths = 0.5,
            s=8, 
            zorder = 3)
#    ax[2][1].scatter(x, residuals,
#            marker='o',
#            color='white',
#            edgecolors = None,
#            linewidths = 0.5,
#            s=32,
#            zorder = 2)
    ax[2][1].set(xlabel= r"Short Time  $/min$",
              ylabel="Residuals",
     #         title = "Lane # "+lane_name,
     #         xlim=[None, None],
              ylim=[-.01, +0.01]
               )
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(plot_file+"_"+lane_name+"_"+row_name+".pdf")     ### export the plot as this
    plt.show()

    print("Plot saved as "+plot_file+"_"+lane_name+"_"+row_name+".pdf")

def contact_sheet(data_root_name,
                  columns = ("1","2","3","4","5","6",
                             "7","8","9","10","11","12"),
                  enzymes = ("","","","","","","","","","","",""),
                  rows = ("A","B","C","D","E","F","G","H"),
                  fancy = False):
    """Plots 12 plots, one for each column of plate
    
    Will produce a 'contact sheet' of plots, one for each plat column
    
    Arguments:
    ---------
    data_root_name: string
        The file root name. Column and row labels will be added.
    columns, row: lists
        These are the lists of names for the columns and rows. Will be appended
        to data_root_name to access the individual data files.
    enzymes: list
        The names of the enzymes (or other note) in each column.
        Will be printed on each plot for corresponding column
    Returns: 
    -------

    Nul

    Will output the 3 column X 4 row plot of all data sets in the plate
    and save the plot as a pdf file.
    """
    plt.ioff()      ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet

    fig, ax = plt.subplots(nrows=4,
                           ncols=3,
                           figsize=(7,10),
                        #   sharex=True,
                        #   sharey=True
                           )
    n = 0  ### set counter
    coldata = zip(columns, enzymes)
    for lane_name, enzyme in coldata:

        plot_row = n // 3   ### use counter to get column and row
        plot_col = n % 3

        ax[plot_row][plot_col].set(
                    xlabel = None,
                    ylabel = None,
            #        title = "Lane # "+lane_name,
            #        xlim = [None, None],
                    ylim = [-.1, 4.1])


        ### Column 0 gets y-axis label and ticks.
        ax[plot_row][0].set(ylabel= r"$A_{405}$")
        ax[plot_row][0].set_yticks([0,1,2,3])
        ### Columns 1 & 2 get no y-axis ticks.
        ax[plot_row][1].set_yticks([])
        ax[plot_row][2].set_yticks([])
    
        ### Row 3 (bottom) gets x-axis label and ticks.
        ax[3][plot_col].set(xlabel= r"$t\;/\;min$")
        ax[0][plot_col].set_xticks([])
        ### Other rows get no x-axis ticks.
        ax[1][plot_col].set_xticks([])
        ax[2][plot_col].set_xticks([])
    
        ax[plot_row][plot_col].text(0, 3.4, "Column: "+str(lane_name))
        ax[plot_row][plot_col].text(0, 3.0, enzyme)

    
        for row_name in rows:
            in_file_name = data_root_name + "_" \
                + str(lane_name) + "_" \
                + row_name + ".csv"
            df = pd.read_csv(in_file_name)
    
            x = df["time"]
            y = df["abs"]
            
            ax[plot_row][plot_col].plot(x, y,
                                        linestyle = '-',
                                        linewidth='0.3',
                                        color = 'black',
                                        zorder = 0)
            #ax[plot_row][plot_col].scatter(x, y,
            #                               marker='o',
            #                               color='black',
            #                               edgecolors = None,
            #                               linewidths = 0.5,
            #                               s=1,
            #                               zorder = 2)
            #ax[plot_row][plot_col].scatter(x, y,
            #                               marker='o',
            #                               color='white',
            #                               edgecolors = None,
            #                               linewidths = 0.5,
            #                               s=4,
            #                               zorder = 1)
        ### end of for row_name in row_name_list:

        n += 1    ### increment counter
    ### end of for lane_name in lane_name_list:

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    ### Last characters of data_root_name is expected to be the plate number
    save_name = "data2/plots/plot_contact_sheet_"+data_root_name[-2:]+".pdf"
    plt.savefig(save_name) ### export the plot as this
    #plt.show()                 ### display the plot in this notebook
    plt.close()
    print(f"Plot saved as {save_name}")

def make_plate_data_for_setup(file_name):
    '''takes a plate setup file name and uses the data to construct
    a plate plan.

    input: file_name: path to csv file. Must contain columns for...
        "Column": Corresponds to column in plate.
        "enzyme": text string of enzyme name
        "E_conc": float of enzyme concentration in nanomolar
        "Row": Row designations on plate/ Usually A to H
        "S_Conc": Substrate concentration in each row in millimolar
        "kcat": kcat value for each enzyme in per second
        "KM": KM value for each enzyme in millimolar
    
    returns: dataframe containing 

        "row_name": same as "Row"
        "S_conc": Substrate concentration converted to molar
        "lane_name": same as the "Column" 
        "E_conc": Enzyme conc converted to molar
        "kcat": kcat converted to per minute
        "Vmax": calculated from kcat and E_conc. In molar per minute
        "KM": KM converted to molar
        "Enzyme": Same as "enzyme" above
    '''
    #Load and clean data file
    df = pd.read_csv(file_name,
                     comment = "#",
                     skipinitialspace=True)
   # df.dropna(axis = 0, inplace = True)   # remove any rows with NaN values

    # Make dictionary of lists for constructing the final dataframe
    dictionary = {}
    dictionary["row_name"] = df["Row"]
    dictionary["S_conc"] = df["S_Conc"] * 1E-3  ### convert millimolar to molar
    dictionary["lane_name"] = df["Column"]
    dictionary["E_conc"] = df["E_Conc"] * 1E-9  ### convert nanomolar to molar
    dictionary["kcat"] = df["kcat"] * 60        ### convert /s to /min
    dictionary["Vmax"] = dictionary["kcat"] * dictionary["E_conc"] ## Vmax from kcat and [E]
    dictionary["KM"] = df["KM"] * 1E-3          ### convert mM to M
    dictionary["Enzyme"] = df["Enzyme"]

    df2 = pd.DataFrame.from_dict(dictionary)
    return df2

def make_data_files_plates(setupdata, file_name, pH = 7.0):
    '''Will take the dataframe of setupdata and output a series of csv files 
    that represent the plate reader data for the experiment
    
    '''

    eq, f = get_integrated_MM_function()  # eq is sympy equation object, f is the corresponding function

    time_start = 0.5        ### start at 30 second (0.5 minutes)
    time_end = 60           ### The end time (minutes)
    n_points = 360          ### number of points - increase if needed

    voltage_error = 0.001   ### parameters to define output range and error
    random_error = 0.0001
    max_value = 4

    dt = time_end / n_points          ### time step, delta t
    t_line = np.arange(time_start,    ### time vector (list of time points)
                    time_end + dt,
                    dt)

    ### Note: Lane names, enzyme conc list, KM list and Vmax list must all be
    ### same length or this will fail. Row names and row concentration lists
    ### must also have equal lengths.

    parameters = zip(setupdata["lane_name"],
                    setupdata["Enzyme"],
                    setupdata["Vmax"],
                    setupdata["KM"])

    for p in parameters:
        lane_name, E_name, Vmax_value, KM_value = p   ### unpack kcat, KM and [E]

        row_df = setupdata[["row_name","S_conc"]].dropna()
        print(row_df)
        row_info = zip(row_df["row_name"], row_df["S_conc"])


        for row in row_info:
            row_name, S0_value = row      ### unpack row name and substrate conc

            plate_df = pd.DataFrame([])   ### start with empty dataframe
            #print(row_name)

            ### Calculate product from enzyme reaction 
            product_E = S0_value - f(t_line, S0_value, KM_value, Vmax_value)
            product_E = np.real(product_E)  ### complex numbers fixed

            ### Calculate product from uncatalyzed reaction
            e_NPA = calculate_e_NPA(pH)

            product_NPA = 0*(S0_value - S0_value * np.exp(-1E-5 * t_line))   #### REMOVE 0 
            product = product_E + product_NPA
            absorbance = product * e_NPA   ### result in absorbance units

            ### Add voltage error 
            fraction_transmittance  = 1 / (10 ** absorbance)                      
            fraction_transmittance = np.random.normal(fraction_transmittance,
                                                    voltage_error, 
                                                    len(fraction_transmittance))

            ### Replace negative transmittance values with a very low value to avoid math errors with np.log10
            fraction_transmittance[fraction_transmittance <= 0] = 1E-6

            ### Calculate absorbance
            absorbance = -np.log10(fraction_transmittance)

            ### Add random error ro absorbance
            absorbance = np.random.normal(absorbance,
                                        random_error,
                                        len(absorbance))

            ### Clean generated data set.
            absorbance[absorbance > max_value] = max_value   ### cap values to maximum absorbance
            absorbance = np.nan_to_num(absorbance,  ### replace any NaN with max value
                            copy = True,
                            nan = max_value)

            ### insert the two data arrays into the dataframe
            plate_df["time"] = t_line
            plate_df["abs"] = absorbance

            ### Write data out to file using lane_name and row_name
            out_file_name = file_name + "_" + str(lane_name) + "_" + str(row_name) + ".csv"
            
            plate_df.to_csv(out_file_name, float_format='%10.4g')

def collect_lanes(data_file_path, plate_name_list, plate_plan_path,
                Fraction_time_span = 1, 
                pH = 7,
                result_file_path = "data1/results/analysis_results_",
                make_plots = True,
                plot_file_path = "data1/plots/",
                fancy = True):

    """Determines slope of initial rate in each cell given time fraction.

    Will fit data to linear curve fit and collect slopes and intercept data for
    all wells according to a list of plate file names. Slopes calculated for
    initial slope given a fraction of the total time span.
    Will save the data as a separate file for each lane.
    
    Arguments
    ---------
    data_file_name: string
        The data file name. Filenames will be this plus column and row
        e.g. "name_10_C.csv"
    Column_list, Row_list: Array like
        list of names for column, row to be plotted.
        Will usually be one column and some rows but can be as many as wanted
    Fraction_time_span: float
        The fraction of the time span to plt. Default is 1 for 100%
    
    Returns
    -------
    result: nul
    Will output .csv files
    """

    def linear_function(x, slope, intercept):
        return slope * x + intercept

    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet

    e_NPA = calculate_e_NPA(pH = pH)

    for plate_name in plate_name_list:

        plate_df = read_plate_setup(plate_plan_path + plate_name + ".csv")

        column_list = plate_df["Column"]
        enzyme_names = plate_df["Enzyme"]
        enzyme_concs = plate_df["E_Conc"]
        row_list = plate_df["Row"].dropna()
        substrate_concs = plate_df["S_Conc"].dropna()

        for lane_name, enzyme_name, enzyme_conc in zip(column_list, enzyme_names, enzyme_concs):

            slope_list = []; slope_stderr_list = []
            int_list = []; int_stderr_list = []; rsq_list = []
            plate_list = []; well_lane_list = []; well_row_list=[]
            enzyme_name_list = []; enzyme_conc_list = []
            substrate_conc_list = []

            for row_name, substrate_conc in zip(row_list, substrate_concs):
                ### get abs vs time data from file
                in_file_name = data_file_path + plate_name \
                                + "_" + str(lane_name) \
                                + "_" + row_name + ".csv"
                df = pd.read_csv(in_file_name)

                ### Trim data list to time fraction
                points_used = int(Fraction_time_span * len(df["time"]))
                x = df["time"][0:points_used]
                y = df["abs"][0:points_used]

                ### Curve fit to get slope in abs/time and intercept in abs
                param, cov = curve_fit(linear_function, x, y)

                rsq = r2_score(y, linear_function(x, *param))
                slope, intercept = un.correlated_values(param, cov)

                ### Convert slope from abs/time to conc/time
                slope = slope / e_NPA    # convert abs per min to molar per min        ############################
                slope = slope * 1E6  # convert molar per min to micromolar per minute  ############################

                intercept = intercept / e_NPA    # convert abs to molar
                intercept = intercept * 1E6  # convert molar micromolar

                slope_list.append(slope.nominal_value)
                slope_stderr_list.append(slope.std_dev)
                int_list.append(intercept.nominal_value)
                int_stderr_list.append(intercept.std_dev)
                rsq_list.append(rsq)
                well_lane_list.append(str(lane_name))
                well_row_list.append(row_name)
                plate_list.append(plate_name)
                enzyme_name_list.append(enzyme_name)
                enzyme_conc_list.append(enzyme_conc)
                substrate_conc_list.append(substrate_conc)

            results = {"Plate": plate_list,
                       "Column":well_lane_list,
                       "Enzyme":enzyme_name_list,
                       "E_Conc":enzyme_conc_list,
                       "Row":well_row_list,
                       "S_Conc":substrate_conc_list,
                       "slope":slope_list,                # micromolar per minute
                       "slope stderr":slope_stderr_list,  # micromolar per minute
                       "int": int_list,                   # micromolar
                       "int stderr":int_stderr_list,      # micromolar
                       "RSQ": rsq_list}
            
            result = pd.DataFrame(results)        
            result.to_csv(result_file_path +"_"+ plate_name +"_"+ str(lane_name) + ".csv")
            print("Data saved as " + result_file_path +"_"+ plate_name +"_"+ str(lane_name) + ".csv")
        print(f"{plate_name} complete")

def plot_wells2(data_file_path, plate_name_list, Column_list, Row_list,
                plate_plan_path = "data2/plateplans/",
                Fraction_time_span = 1,
                Make_Plots = False, 
                plot_file_path = "plots/analysis_plot_",
                Save_Data = False, 
                result_file_path = "data1/data/analysis_results",
                fancy = False, tiny_points = False, tiny_line = False):

    """Loads and plots data files. Can return line fit.

    Will make a plot of all the wells designated by columns and rows
    Will fit data to linear curve fit and collect slopes and intercept data
    Will save the plot as a pdf
    Will save the data as a csv file with same filename as pdf
    
    Arguments
    ---------
    data_file_name: string
        The data file name. Filenames will be this plus column and row
        e.g. "name_10_C.csv"
    Column_list, Row_list: Array like
        list of names for column, row to be plotted.
        Will usually be one column and some rows but can be as many as wanted
    Fraction_time_span: float
        The fraction of the time span to plt. Default is 1 for 100%
    Line_Fit: boolean
        If True then line fits will be performed for each well and data
        written to lists and put into the result dataframe
    Display_Plot, Display_Data: booleans
        If True then the plot will be created and displayed and the
        dataframe displayed, respectively.
    
    Returns
    -------
    result: pandas dataframe
        The results of line fits. Will be empty if Line_Fit = False
    """

    def linear_function(x, slope, intercept):
        return slope * x + intercept
    
    #print(Column_list)

    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    

    slope_list = []; slope_stderr_list = []
    int_list = []; int_stderr_list = []; rsq_list = [];
    plate_list = []; well_lane_list = []; well_row_list=[];

    for plate_name in plate_name_list:
        #print(plate_name)

        plateplan_df = read_plate_setup(plate_plan_path + plate_name + ".csv")

        mask = plateplan_df["Column"].isin(Column_list)
        column_df = plateplan_df[mask][["Column", "Enzyme", "E_Conc"]]

        mask = plateplan_df["Row"].isin(Row_list)
        row_df = plateplan_df[mask][["Row", "S_Conc"]]

        column_name_list = column_df["Column"]
        enzyme_names = column_df["Enzyme"]
        enzyme_concs = column_df["E_Conc"]
        row_list = row_df["Row"]
        substrate_concs = row_df["S_Conc"]

        for lane_name, e_name, e_conc in zip(column_name_list, enzyme_names, enzyme_concs):
        #print(lane_name)

            slope_list = []; slope_stderr_list = []
            int_list = []; int_stderr_list = []; rsq_list = [];
            plate_list = []; well_lane_list = []; well_row_list=[]
            E_conc_list = []; S_conc_list = []; E_name_list=[]


            for row_name, s_conc in zip(row_list, substrate_concs):
                in_file_name = data_file_path + plate_name \
                                + "_" + str(lane_name) \
                                + "_" + row_name + ".csv"
                df = pd.read_csv(in_file_name)

                points_used = int(Fraction_time_span * len(df["time"]))
                x = df["time"][0:points_used]
                y = df["abs"][0:points_used]

                param,cov = curve_fit(linear_function, x, y)
                slope, intercept = un.correlated_values(param, cov)

                rsq = r2_score(y, linear_function(x, *param))

                slope_list.append(slope.n)
                slope_stderr_list.append(slope.s)
                int_list.append(intercept.n)
                int_stderr_list.append(intercept.s)
                rsq_list.append(rsq)
                well_lane_list.append(str(lane_name))
                well_row_list.append(row_name)
                plate_list.append(plate_name)
                E_conc_list.append(e_conc)
                E_name_list.append(e_name)
                S_conc_list.append(s_conc)

                if Make_Plots:
                    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))

                    x_fit = np.linspace(0,np.max(x),10)
                    ax.plot(x_fit, linear_function(x_fit, *param),
                            linestyle = '-',
                            linewidth='0.5',
                            color = 'black',
                            zorder = 0)

                    if tiny_line:
                        ax.plot(x, y,
                            #marker=None,
                            color='black',
                            linewidth = 0.5,
                            ms=8,
                            zorder = 1)
                    if tiny_points:
                        ax.scatter(x, y,
                            marker='o',
                            color='black',
                            edgecolors = 'none',
                            linewidths = 0.5,
                            s=2,
                            zorder = 2)
                    else:
                        ax.scatter(x, y,
                            marker='o',
                            color='white',
                            edgecolors = 'black',
                            linewidths = 0.5,
                            s=32,
                            zorder = 2)
                        ax.scatter(x, y,
                            marker='o',
                            color='white',
                            edgecolors = None,
                            linewidths = 0.5,
                            s=64,
                            zorder = 1)
                    ### end of if tiny_points:
                    if np.min(y) < 0:
                        q = np.min(y) - 0.05*np.max(y)
                    else:
                        q = -0.05*np.max(y)

                    ax.set(xlabel= r"Time $/min$",
                           ylabel=r"$A_{405}$",
                           #title = "Lane # "+lane_name,
                           xlim=[0, None],
                           ylim=[q, None]
                            )
                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.savefig(plot_file_path +"_"+ plate_name +"_"+ str(lane_name) \
                                +"_"+ row_name + ".pdf")     ### export the plot as this
                    plt.close()
                ### end of if plots
            ### end of for row_name:


            results = {"Plate": plate_list,
                       "Column": well_lane_list,
                       "E_name": E_name_list,
                       "E_conc": E_conc_list,
                       "Row": well_row_list,
                       "S_conc": S_conc_list,
                       "slope": slope_list,
                       "slope stderr": slope_stderr_list,
                       "int": int_list,
                       "int stderr": int_stderr_list,
                       "RSQ": rsq_list
                       }
            results = pd.DataFrame(results)

                #display(results)

            if Save_Data:
                results.to_csv(f"{result_file_path}_{plate_name}_{lane_name}.csv")
                print(f"Data saved as {result_file_path}_{plate_name}_{lane_name}.csv")

def plot_lanes_MM(result_file_path, plate_name_list, Column_list, 
                    final_file_path = "data1/kinetic_values",
                    Make_Plots = False,
                    plot_file_path = "plots/lanes_plot_",
                    fancy = True,
                    pH = 7.0
                    ):
    
    """Loads and plots data files. Can return line fit.

    Will make a plot of all the lanes given for the plates given
    Will fit data to MM curve fit and collect KM and kcat data
    Will save the plots as a pdf if selected
    Will save the result as a csv file with same filename as pdf
    
    Arguments
    ---------
    result_file_path: string
        The path to result files. Filenames will be this plus column and row
    plate_name_list: Array like
        The names of the plates used to create the result files
        e.g. "name_10_C.csv"
    Column_list: Array like
        list of names for lanes in each plate to be used
    Make_plots: boolean
        If True then curve fits will be performed for each file and data
        written to lists and put into the result dataframe
    plot_file_path: string
        path and root name for plot files. Plat number and lane label will
        be appended to this string.
    fancy, tiny_points, tiny_line: booleans
        Style control.
    
    Returns
    -------
    result: pandas dataframe
        The results of curve fits. Will be empty if Line_Fit = False
    """
    def linear_function(x, slope, intercept):
        return slope * x + intercept
    
    #print(Column_list)

    plt.ioff()           ### switch off interactive display of plots. plt.show() needed to display a plot now
    plt.rcdefaults()     ### resets the plot defaults so we always start in the same place
    if fancy:
        plt.style.use("../styles/tufte.mplstyle")     ### Then add a fancy style sheet
    
    e_NPA = calculate_e_NPA(pH = pH)
    print(f"eNPA = {e_NPA:0.2f}")
    KM_list = []; KM_U_list = []
    kcat_list = []; kcat_U_list = []
    rsq_list = []
    plate_list = []; lane_list = []
    kcat_KM_list = []; kcat_KM_U_list = []
    enzyme_list = []; E_conc_list = []

    for plate_name in plate_name_list:
        #print(plate_name)
        for lane_name in Column_list:
        #print(lane_name)

            in_file_name = result_file_path + plate_name \
                            + "_" + str(lane_name) + ".csv"
            df = pd.read_csv(in_file_name)
            y = df["slope"]
            y_u = df["slope stderr"]
            x = df["S_conc"]
            
            enzyme = df["E_name"][0]
            E_conc = df["E_conc"][0]

            (param, cov, *other) = curve_fit(MM, x, y,
                                   p0=None, 
                                   sigma=y_u, 
                                   absolute_sigma=True, 
                                   bounds=(0, np.inf)
                                )
            rsq = r2_score(y, MM(x, *param))

            Vmax_u, KM_u = un.correlated_values(param, cov)
            print(Vmax_u, KM_u)

            Vmax_u_molar_s = Vmax_u / e_NPA     # from abs per min to molar per minute
            E_conc_molar = E_conc / 1E9       # from nanomolar to molar
            kcat_u = (Vmax_u_molar_s / 60) / E_conc_molar   # Convert Vmax to molar per second and get kcat in per second 
            KM_u_molar = KM_u /1000    # From millimolar to molar
            kcat_KM_u = kcat_u / KM_u_molar

            KM_list.append(KM_u.n)              # millimolar
            KM_U_list.append(KM_u.s)
            kcat_list.append(kcat_u.n)          # per second
            kcat_U_list.append(kcat_u.s)
            rsq_list.append(rsq)
            plate_list.append(plate_name)
            lane_list.append(lane_name)
            kcat_KM_list.append(kcat_KM_u.n)   # per molar per second
            kcat_KM_U_list.append(kcat_KM_u.s)
            enzyme_list.append(enzyme)
            E_conc_list.append(E_conc)          # nanomolar

            ### end of if:
            if Make_Plots:
                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
                x_fit = np.linspace(0,np.max(x),100)
                ax.plot(x_fit, MM(x_fit, *param),
                        linestyle = '-',
                        linewidth='0.5',
                        color = 'black',
                        zorder = 0)
                ax.scatter(x, y,
                    marker='o',
                    color='white',
                    edgecolors = 'black',
                    linewidths = 0.5,
                    s=16,
                    zorder = 2)
                ax.errorbar(x, y, yerr=y_u * 3, xerr=None, fmt="None", ecolor="red", elinewidth = 0.5, zorder=1)

                ax.set(xlabel= r"[S] $/mM$",
                       ylabel=r"rate $\mu M/min$",
                       title = f"{enzyme}, {plate_name}, Lane {lane_name}",
                       xlim=[-0.05*np.max(x), np.max(x)*1.05],
                       ylim=[-0.05*np.max(y), np.max(y)*1.05]
                        )
                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                plt.savefig(f"{plot_file_path}_{plate_name}_{lane_name}.pdf")     ### export the plot as this
                plt.close()
        print(f"{plate_name} complete")

    results = {"Plate": plate_list,
               "Column":lane_list,
               "Enzyme": enzyme_list,
               "E_conc": E_conc_list,
               "KM (mM)":KM_list,
               "KM stderr":KM_U_list,
               "kcat (/s)": kcat_list,
               "kcat stderr":kcat_U_list,
               "kcat/KM (M/s)": kcat_KM_list,
               "kcat/KM stderr":kcat_KM_U_list,
               "RSQ": rsq_list}
    
    results = pd.DataFrame(results)

    #display(results)

    results.to_csv(result_file_path + ".csv", float_format='%10.4g')
    print("Data saved as " + result_file_path + ".csv")

    return results

