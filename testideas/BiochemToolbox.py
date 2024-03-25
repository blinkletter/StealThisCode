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

    return(Vmax * S / (S + KM))
        
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

    return(v_max_u, KM_u)     

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
    return(slope * x + intercept)

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
    
    return(slope_u, intercept_u)

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
    return(eq,f)

def read_plate_setup(file_name, pH = 7.0):
    """Reads a file name for a plate plan and returns lists of the data

    Arguments:
    ----------
    file_name: string
    pH: float

    Returns:
    ----------
    pandas dataframe
        The dataframe read in from the file. Some columns are different 
        lengths so it will be separated into lists (Pandas series) and 
        those are returned in the dictionary below.
    dictionary 
        {"row_name_list",
         "S_conc_list",
         "lane_name_list",
         "E_conc_list",
         "E_Name_list"}
    float
        The extinction coefficent of nitrophenolate at the given pH

    Example
    -------
    df, plate, e405 = read_plate_setup{"platesetup.csv", pH = 7.0}
    """
    
    df = pd.read_csv(file_name,
                     comment = "#",
                     skipinitialspace=True)

    row_name_list = df["Row"].dropna()
    S_conc_list = df["S_Conc"].dropna() * 1E-3  ### convert millimolar to molar

    lane_name_list = df["Column"]
    E_conc_list = df["E_Conc"] * 1E-9  ### convert nanomolar to molar
    E_Name_list = df["Enzyme"]

    ### parameters to get extinction coeff for NPA at give pH value
    pH = 7.0
    e_NPAA = 18300  ### extinction coeff for NPA anion
    pKa_NPA = 7.15 ### pKa for p-nitrophenol

    ### Calculated Values from the above lists
    Ka = 10 ** -pKa_NPA   ### extinction coeff for NPA at given pH
    H = 10 ** -pH
    e_NPA = e_NPAA * (Ka / (H + Ka))

    return df, {"row_name_list": row_name_list,
            "S_conc_list": S_conc_list,
            "lane_name_list": lane_name_list,
            "E_conc_list": E_conc_list,
            "E_Name_list": E_Name_list,}, e_NPA

    #display(df)

def plot_lanes(ax, data_file_name, Column_list, Row_list, 
               Fraction_time_span = 1, Line_Fit = True):
    """Loads and plots data files. Can return line fit.
    
    Arguments
    ---------
    ax: pyplot axes object
        axes created and sent to this function for adding the points
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
    
    Returns
    -------
    ax: pyplot axes object
        The plot axes will new elements added
    result: pandas dataframe
        The results of line fits. Will be empty if Line_Fit = False
    """
    def linear_function(x, slope, intercept):
        return slope * x + intercept
    
    #print(Column_list)

    slope_list = []; slope_stderr_list = []
    int_list = []; int_stderr_list = []; rsq_list = [];
    well_lane_list = []; well_row_list=[]

    for lane_name in Column_list:
        #print(lane_name)

        for row_name in Row_list:
            in_file_name = data_file_name \
                            + "_" + str(lane_name) \
                            + "_" + row_name + ".csv"
            df = pd.read_csv(in_file_name)
            points_used = int(Fraction_time_span * len(df["time"]))
    
            x = df["time"][0:points_used] 
            y = df["abs"][0:points_used]
    
            
            if Line_Fit:
                param,cov = curve_fit(linear_function, x,y)
                slope, intercept = param

                perr = np.sqrt(np.diag(cov))
                slope_stderr, int_stderr = perr
                r, p = pearsonr(x, y)
                rsq = r ** 2

                slope_list.append(slope)
                slope_stderr_list.append(slope_stderr)
                int_list.append(intercept)
                int_stderr_list.append(int_stderr)
                rsq_list.append(rsq)
                well_lane_list.append(lane_name)
                well_row_list.append(row_name)

                x_fit = np.linspace(0,np.max(x),10)
                ax.plot(x_fit, linear_function(x_fit, slope, intercept), 
                        linestyle = '-', 
                        linewidth='0.5', 
                        color = 'black', 
                        zorder = 0)
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
    
    results = {"Column":well_lane_list,
               "Row":well_row_list,
               "slope":slope_list,
               "slope stderr":slope_stderr_list,
               "int": int_list,
               "int stderr":int_stderr_list,
               "RSQ": rsq_list}
    results = pd.DataFrame(results) 
     
    return(ax, results)
