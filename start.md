
# Steal This Code

This book contains workbooks that outline various methods of data analysis and plotting with *Python*. It is a collection of *Python* notebooks that I have written over the years for various purposes. The major purpose of this is to provide a library of my own code and documnetation for myself so that I can more easily re-use the work that I have done and not have to relearn and re-invent wheels when I need to solve a problem that I had already solved a few years prior. This record of examples is meant to be mined for useful bits. Steal some of it or steal all of it, just don';t spend time writing code that someone else already has written.

## Programming as Theft

We will use *Python* only as a **tool**. All code will be written for you in the first sets of exercises. You will need only change a few values to adapt the code to your needs. Once you have experienced changing and **reusing** code that was written for you you will be asked to **steal** it and **adapt** it to solve the final problems in this course. You will never need to change much, but will hopefully gain the ability to find useful code, steal it and **use it**. Steal it from me; steal it from the internet; give your code to friends. Always **give credit** where credit is due.

## Things to Steal

**All the code** that you will need in the first half of our exercises will be **available** in the information pages for the lab or tutorial in question. This book is intended to provide some basic **examples** of *Python* code for performing **calculations and visualizing data**. You will see that the programs that you are using during our tutorials were all stolen from this book and changed to suit you. Examine, play with, and **steal the code** in the interactive notebooks linked below.

### Math and Calculations

>**Title**: [Basic Math with *Python*](1_Math/1_1_Calculator.ipynb) <br>
>**Tools**: *Python* math operators, *NumPy*  <br>
>**Skills**: Basic math operations and math functions using the tools of *NumPy*. Math with *NumPy* arrays, formatting numbers in f-strings.   

>**Title**: [Using *Uncertainty*](1_Math/1_2_uncertainties.ipynb) <br>
>**Tools**: *Uncertainties*, *NumPy*  <br>
>**Skills**: Math operations with uncertain values, using the tools of *NumPy* and their versions provided by the *Uncertainties* package.   

### Plotting with MatPlotLib

>**Title**: [Plotting Data](2_MatPlotLib/2_1_Plotting.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*  <br>
>**Skills**: Making simple x,y data plots. Styling plots  

>**Title**: [Linear Fits](2_MatPlotLib/2_2_Curve_Fits.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*, *SciPy*  <br>
>**Skills**: Linear parameter fits to x,y data. Styling plots   

>**Title**: [Curve Fits](2_MatPlotLib/2_3_Curve_Fits.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*, *SciPy*  <br>
>**Skills**: Curve fits to x,y data. Styling plots. Some more advanced styling.   

>**Title**: [Advanced Plotting Example](2_MatPlotLib/2_4_Fancy_1.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*, *SciPy*  <br>
>**Skills**: All of the above plus advanced plotting using axes objects. Two plots in one figure.   

>**Title**: [Notes on Plotting](2_MatPlotLib/2_5_Notes.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*, *SciPy*, *Pandas*  <br>
>**Skills**: More examples of plots and curve fits. A residual plot is demonstrated. An example of a curve fit with the function plotted over the data is presented. 

>**Title**: [Example: Order of Reaction](2_MatPlotLib/2_6_Plotting_Reaction_Kinetics_to_Determine_Order.ipynb) <br>
>**Tools**: *MatPlotLib*, *NumPy*, *SciPy*, *Pandas*  <br>
>**Skills**: More examples of plots and curve fits. Using plots and line fits to determine the order of a reaction. 


### Calculus with SymPy

>**Title**: [Symbolic Math](3_SymPy/K01_solving_with_sympy.ipynb) <br>
>**Tools**: *SymPy*, *NumPy* <br>
>**Skills**: Symbolic algebra, solving for unknowns automatically.  

>**Title**: [Integrating Rate Laws](3_SymPy/K02_Integrating_Rate_Law_1.ipynb) <br>
>**Tools**: *SymPy*, *NumPy*, *MatPlotLib* <br>
>**Skills**: Calculus with *SymPy*. Integrating 1<sup>st</sup> and 2<sup>nd</sup> order rate laws. Plotting algebra expressions with *SymPy*. Creating and plotting functions made from *SymPy* expressions.

>**Title**: [Integrating the Michaelis-Menten Rate Law](3_SymPy/K03_Integrating_MM.ipynb) <br>
>**Tools**: *SymPy*, *NumPy*, *SciPy*, *LMFit*, *MatPlotLib* <br>
>**Skills**: Advanced calculus without having a clue about calculus with *SymPy*. Curve fitting data against functions made from *SymPy* expressions using `sympy.curve_fit` and `lmfit.Model`.  

### Numerical Methods with *SciPy*

>**Title**: [Numeric Integration](4_SciPy/K04_NumericIntegration.ipynb) <br>
>**Tools**: *SymPy*, *NumPy*, *SciPy*, *MatPlotLib* <br>
>**Skills**: Integrating a differential equation using numerical methods. Comparing the results to the analytical symbolic integration from the previous notebook. Curve fitting the data against a function that performs numerical integration.  

>**Title**: [Integrating Systems of Equations](4_SciPy/K05_NumericIntegration_SystemEq.ipynb) <br>
>**Tools**: *SymPy*, *NumPy*, *SciPy*, *MatPlotLib* <br>
>**Skills**: Going beyond the MM equation and its simplifications and assumptions to integrate the system of rate laws that fully describes the mechanism of Michaelis-Menten model for enzyme catalysis.  Comparing the results to those from integrating the Michaelis-Menten equation.

### PhysOrg Examples

This set of five notebooks explores some examples of the methods presented above. These notebooks supported the introduction to *Python for my Physical Organic Class. You can read a description of each exercise by going to the **[Intro to PhysOrg Examples](6_PhysOrgExamples/6_0_Introduction.md)** page.

---
These is a **Juptyter-book** that was built from a set of **interactive *Python* Jupyter notebooks**. The original notebook for any given chapter can be obtained using the **download link** at the top of the page.

---
The Ouroboros with benzene image was accessed at https://commons.wikimedia.org/wiki/File:Ouroboros-benzene.svg on June 24, 2022. It is licensed under the Creative Commons license CC BY-SA 3.0.