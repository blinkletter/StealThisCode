# The Erying Equation and Experimental Error
The famous Erying equationis used to create a linear equation that relates the free energy of activation or the free energy in an equilibrium with temperature. Because the enthalpy term is independant of temperature and the entropy term involves temperature, we can the separate these effects and determine the enthalpy and the entropy difference in the chemical change between reactants and transtion state or products. But, do you numbers make sense? Unless you are calculating the standard deviations from your plot you may be unaware of a low quality relationship.

## Enzyme Kinetics Example

The data sets for this experiment comes from "Linear Eyring Plots Conceal a Change in the Rate-Limiting Step in an Enzyme Reaction", Teresa F. G. Machado, Tracey M. Gloster, and Rafael G. da Silva, *Biochemistry*, **2018**, *57*, 6757-6761.https://doi.org/10.1021/acs.biochem.8b01099.

Machado et al. were presenting an example where they hypothesized that there was a hidden effect due to changes in the rate-determining step. The plots are very linear so we will use the data at face value. I'm not really convinced that there was a change in rds.

The reaction is the reduction of acetoacetate (3-oxobutanoate) or 3-Oxopentanoate by NADH with the action of *(R)-3-hydroxybutyrate dehydrogenase* (EC 1.1.1.30). Two sources of the enzyme were used. The common bacteria *Acinetobacter baumannii* and the cold-favouring (psychrophilic) bacteria *Psychrobacter arcticus*. The Michaelis-Menten turnover number ($k_{cat}$, the first-order rate constant of the enzyme-NADH-substrate complex to give product) was measured at temperatures for $-10$ to $67\ ^\circ C$. The Eyring plot for the *P. arcticus* curves slightly as temperature increases. The enzyme comes from an organism adapted got cold conditions. The plot for the *A. baumannii* enzyme were linear so I will only use that data in this example.


## The Eyring Equation

The Erying equation is...

$$ k = \frac{\kappa k_B}{h} T e^{-\Delta G^\ddagger}$$

...and we can include that $\Delta G^\ddagger = \Delta H^\ddagger - T \Delta S^\ddagger$ and we obtain...

$$ k = \frac{\kappa k_B}{h} T e^{\frac{-\Delta H^\ddagger}{RT}}e^{\frac{\Delta S^\ddagger}{R}}$$

This can be written in a linear form...

$$\ln \frac{k}{T} = \frac{-\Delta H^\ddagger}{R}\cdot\frac{1}{T} + \ln\frac{\kappa k_B}{h} + \frac{\Delta S^\ddagger}{R}$$

...where $k$ is the observed rate constant and $T$ is the absolute temperature.

From this, we can find $\Delta H^\ddagger$ and $\Delta S^\ddagger$ from the slope and intercept by plotting $\ln \frac{k}{T}$ vs. $\frac{1}{T}$

## Know Your Error

Or rather, estimate your error correctly. It's easy to get a standard deviation for the set of results from a repeated experiment. How do you get the standard deviation for parameters derived from a line fit? Microsoft Excel can do this. But to fit an arbitrary equation, like a curve fit of the Erying equation, you will need to use the Solver ad-on, which is difficult to automate. Python can do all of this easily using many possible options. 

Using the tools described in the pages that follow, we will be able to fit daya to linear or non-linear equation models, obtain estimates of errors for the fitted parameters, use experimental error determined for data points to weight the curve fits and to propagate error through complex equations. Error propagation is often ignored in data analysis because it is complicated and opaque. We will use Python as a 'black box' that handles errors correctly. It will still be opaque, but it wont be complicated.

## Activities

The following pages are interactive *Python* notebooks that are designed to run in Google *Colab*. Use them to explore python as a data anaysis tool.

1. [Short and Sweet](temp)

This notebook will plot a 4-point Erying plot and determine the $\Delta H^\ddagger$ and $\Delta S^\ddagger$ for the reaction. The standard deviations for the parameters will be examined and used to calculate the confidence interval for the rate constant at a given temperature. The ```scipy.stats.linregress``` tool will be demonstrated along with the ```matplotlib.pyplot``` library for plotting. Along the way we will learn the correct way to use the ```uncertainties``` module in our calculations.

2. [A Better Fit](temp)

The ```lmfit``` library provides a better set of tools for handling data and interfacing with the ```uncertainties``` package. It has many built-in tools for analyzing the curve fit. We will be using ```lmfit``` rather than ```scipy.optimize.curvefit``` in the rest of this exercise.

3. [Read it and Don't Weep](temp)

This notebook will read in data from a text file. Creating data files for your experimental data will allow you to use a notebook to analyze a new set of data by changing the file name and nothing else. This will enable your data analysis notebook to be more versatile. More tools in ```lmfit``` for presenting data will be explored.

4. [Style Matters](temp)

Journals often have a specific style for plots. In a thesis you should use the same style for all plots. In this notebook we will explore styling plots using many option available in the ```matplotlib.pyplot``` library. Once you get a style you like in a notebook, you never need to change. I haven't changed my style since 1985, that why I look so good.

5. [By Your Bootstrap](temp)

The confidence intervals produced by ```lmfit``` follow typical rules for error propagation and present a range that is symmetrical above and below any value. This is the common $a \pm b$ way of presenting error. However, the confidence range may not be symmetrical (this is especially true in log scales.) One way to express this is to use more sophisticated error propagation tools in $python$ such as ```soerp``` or ```mcerp``` that can handle second-order effects in error propoagation that become important in complicated systems. You have a whole lifetime ahead of you to explore these tools so we will not start here.

In this notebook we will explore a robust but inefficient method for determining the confidence interval for a data point and generate a confidence band on a plot that reflects the "real world" error in your data. This method is called "bootstrapping."  

## Steal This Code

All of the notebooks presented here are for you to explore and use. Feel free to take any code and repurpose it for your own work. We are not programers. We pick and hack at code we find and hopefully get it working for our own needs. Once we build a tool that works we can use it and start forgetting how to program. The real benefit is that others can follow exactly what we did in our data analysis buy examining the *Python* code. No one will wonder why they cannot exactly obtain the results that you gleaned from your data. The cod ewill tell them exactly how it was done. 

Good documentation is the backbone of science. Repeatable experiments require complete information. A statement like "analyzed using standard methods" does not reveal that you used linear error propagation (the ```Uncertainties```  library) when perhaps a Monte-Carlo method for error propagation (the ```mcerp``` library) would have been more appropriate. Your readers will just be guessing at what you might have done. Provide you code in supplemental material of a paper or in the appendices of your thesis. Then others can evaluate your analysis and properly praise your work or find your error.