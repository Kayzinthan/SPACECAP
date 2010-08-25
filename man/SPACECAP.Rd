\name{SPACECAP}
\alias{SPACECAP}                                                      

\title{A Program to Estimate Animal Abundance and Density using Bayesian Spatially-Explicit Capture-Recapture Models}

\description{
\pkg{SPACECAP} is a user-friendly software package for estimating animal densities using closed model capture-recapture 
sampling specifically based on photographic captures. It implements a new generation of spatially explicit capture-recapture 
models developed in a recent paper, Royle et al.(2009), cited below. Spatially explicit capture-recapture models implemented in 
\pkg{SPACECAP} directly estimate animal density by explicitly using the information on capture histories in combination 
with spatial locations of captures under a unified Bayesian modeling framework. This approach offers advantage such as: 
substantially dealing with problems posed by individual heterogeneity in capture probabilities in conventional 
capture-recapture analyses. It also offers non-asymptotic inferences which are more appropriate for small samples of 
capture data typical of photo-capture studies.
}

\usage{
SPACECAP()
}

\details{

\bold{PROGRAM SPACECAP Version 1.0.1}

\bold{Introduction}

\pkg{SPACECAP}, is a user-friendly software package for estimating animal densities using closed model capture-recapture 
sampling specifically based on photographic captures. It implements a new generation of spatially explicit 
capture-recapture models developed in a recent paper titled 
"Bayesian Inference in Camera Trapping Studies for a Class of Spatially Explicit Capture Models" by J. A  Royle, 
K. U. Karanth, A. M. Gopalaswamy, and N. S. Kumar, published in \emph{Ecology},  \bold{90(11)}, 2009, pp. 3233-3244.

The conventional approach to the analysis of animal density from camera trap surveys is to apply closed 
capture-recapture model analyses, and, then convert resulting estimates of abundances to densities using a wide range 
of essentially ad hoc methods. For example, a standard practice in camera trapping studies has been to use auxiliary 
capture location information to estimate the mean or maximum distance moved by the study species as a basis to apply a 
"buffer width" around the trap array or, to use various other heuristic "adjustments" (Wilson and Anderson 1985a, b, 
Karanth and Nichols 1998, Parmenter et al. 2003, Trolle and Kery 2003) to estimate the effectively sampled area. While 
these approaches appear to work adequately in practice, little had been known about the range of conditions under which 
they work well. This is because most real world study situations involve study area of odd shapes and sizes and difficult 
terrain that makes setting camera traps challenging and conditions assumed by ad hoc approaches may not apply. 

Spatially explicit capture-recapture models implemented in \pkg{SPACECAP} package directly estimate animal density by 
explicitly using the information on capture histories in combination with spatial locations of captures under a unified 
Bayesian modeling framework. This approach offers advantage such as: substantially dealing with problems posed by 
individual heterogeneity in capture probabilities in conventional capture-recapture analyses. It also offers 
non-asymptotic inferences which are more appropriate for small samples of capture data typical of photo-capture studies. 
Further technical details about the models and analyses can be obtained from the paper by Royle et al. (2009), cited above. We also note another spatially explicit modeling and estimation approach based on conventional likelihood based inference that is available (Borchers and Efford 2008).

\bold{The Spatial Capture-Recapture Model and its Parameters}

The model considered in the current version of \pkg{SPACECAP} applies to binary observations y(i,j,k) for individual "i",  
trap "j" and sample occasion (e.g., night) "k". The model is a type of binary regression model, similar to logistic 
regression, in which 

    \deqn{y(i,j,k) ~ Bernoulli(p(i,j,k))}

Here p(I, j, k) is the probability of detecting an individual "i", at trap "j" and sample occasion "k". The probability 
p(i,j,k) is then related to suitable covariates of interest by applying a suitable transformation. In \pkg{SPACECAP} we make 
use of the complementary log-log link transformation. Thus, the simplest possible model is expressed 

    \deqn{cloglog(p(i,j,k)) = b0} 
        
where b0 is a parameter to be estimated. By taking the inverse of the cloglog transformation 
    
    \deqn{p(i,j,k) = 1-exp(-exp(b0))} 
    
This particular function arises by considering the binary observations to be formally reduced from a Poisson encounter 
frequency model (see Royle et al 2009 for details). That is, p(i,j,k) is the Pr(at least 1 encounter in a trap) under 
the Poisson model. The most general model that \pkg{SPACECAP} presently allows is a model in which 
    
    \deqn{cloglog(p(i,j,k)) = b0 + b1*x(i,j,k) + b2*dist(s(i),u(j))^2}
    
In this expression x(i,j,k) is an indicator of previous capture of individual i in trap j. Thus "b1" is a measure of 
the behavioral response (see below for more discussion of this; also see Royle et al 2009 for more details).  The parameter 
b2 is a regression coefficient on the effect of distance between individual activity center s(i) and the location of trap j, 
u(j). b2 is constrained to be < 0 implying that the probability of encountering an individual in a trap decreases as the 
distance between the individual and the trap increases. In the model, s(i) are the collection of animal locations within a 
prescribed region S (this is provided by the user - see below). S is referred to subsequently as the state-space of the 
random variables s(i). 

In \pkg{SPACECAP}, the parameter N is the population size of individuals - the number of activity centers located in S 
(see further discussion of this below).  

In \pkg{SPACECAP}, N is reported under the name Nsuper. \cr
Density, D = N/||S|| where ||S|| is the area of the state-space. The units for Density, D is animals/100 km2. 

\pkg{SPACECAP} also reports several derived parameters: 

lam0 = exp(b0), is the intercept in terms of expected encounter frequency; this may be thought of as the expected 
       encounter rate of an individual "i" in trap location "j" at sampling occasion "k", whose home-range centre 
       is exactly at the trap location.  

sigma = sqrt(1/b2), is the scale parameter of a bivariate normal encounter function (or an exponential, but this 
        option isn't provided in \pkg{SPACECAP} yet). This may also be viewed as a "range parameter" of an animal. 
        For a highly mobile animal, this value will tend to be large (eg:- tigers will have a higher sigma value 
        in comparison to civets)

beta = b1, is the regression coefficient that measures the behavioral response

psi = the ratio of the number of animals actually present within S to the maximum allowable number 
      (set by the user during data augmentation - see below). 

\bold{Installing PROGRAM SPACECAP Version 1.0.1}

STEP 1: Download latest version of program \R (R Development Core Team)

\pkg{SPACECAP} works within the \R programming environment, so your very first step will be to connect to the internet, 
go to the website http://www.r-project.org, download and install the latest version of \R (R 2.9.2 or higher)  from the 
nearest CRAN mirror. 


STEP 2: Download the package \pkg{SPACECAP} to your computer.

When you launch \R, go to Packages->Install package(s), once again select your nearest CRAN mirror and select package 
\pkg{SPACECAP} for installation. 


STEP 3: OPEN the package \pkg{SPACECAP}

In the \R environment, go to Packages->Load package, select \pkg{SPACECAP} and load the package. 


STEP 4: Launching \pkg{SPACECAP}

In the \R environment, at the prompt ">", type the command SPACECAP(). 
This will launch the Graphic-User Interface of SPACECAP and you are now set to begin the Bayesian Spatially-Explicit 
Capture-Recapture (SECR) analysis of your camera trap survey data.
                                                                  
\bold{SECR Analysis using \pkg{SPACECAP}}

Running an SECR Analysis in \pkg{SPACECAP} essentially involves four simple steps: 

1. Setting up the input files \cr
2. Selecting the appropriate model combination \cr
3. Selecting the Markov chain Monte Carlo (MCMC) settings \cr
4. Hitting the "RUN" button. \cr


\bold{STEP 1: SETTING UP THE INPUT FILES FOR ANALYSIS}

\pkg{SPACECAP} requires you to first create three input files and store these files at a suitable location on your computer. 
These files are:

1. Animal Capture Details File (Animal ID no., Trap Location no., Sampling Occasion no.) 

2. Trap Deployment details File (Trap Spatial Location, Deployment Activity, Sampling Occasion no.)

3. State-space details File (describing the Potential Animal Home Range Center Details in terms of their spatial 
location and habitat suitability indicator for these home range centers)

These three raw data files can most easily be created using spreadsheet applications like Microsoft EXCEL, 
OpenOffice or other software you are comfortable with. However, eventually, all files must be saved in an ASCII comma 
separated format \emph{(.csv)}, because \pkg{SPACECAP} can only read these types of input files. 

Upon launching \pkg{SPACECAP}, you will notice on the input data panel, three separate buttons to load the three input 
data files. Pressing of these buttons will help you locate the corresponding input files using your desktop browser 
(like Windows Explorer or Finder). 

\bold{INPUT FILE 1: Animal Capture Details} 

The input file containing individual animal capture histories and locations consists of a 3-column table, each column 
representing the Location Number, the Animal Identity Number and the Sampling Occasion number, strictly in that order. 
Please note that these are all "number" fields. Please not enter labels containing alpha-numeric characters such as 
"CAT-123", "PLACE-100" or "January 2009" etc. for these fields, because \pkg{SPACECAP} will not recognize them. Use simple 
integer numbers. Each unique individual captured during camera trap sampling should be given a unique identification 
number, ranging from 1 to n, where n is the total number of unique individuals caught during the camera trap survey. 

Duration of the overall survey is determined by species biology in order to meet the assumption of demographic closure. 
The duration of sampling occasions (or periods) will, in turn, be based on how many such occasions are there in the 
survey duration. 

Because of shortage or camera traps or logistical reasons in many camera trap surveys, the study area is quite often 
divided into blocks or sub-units and camera trapping is conducted sequentially in these (Karanth and Nichols 2002). 
With data from conventional capture-recapture analyses the resulting capture history matrix is artificially constructed 
to ensure that assumption that camera trap survey simultaneously covers the entire area on each sampling occasion 
(See Karanth and Nichols 2002 page 132-133).

However, since spatial locations of traps are explicitly incorporated into capture-recapture modeling in \pkg{SPACECAP} this 
artificial construction of sampling occasions is not required, which greatly enhances flexibility in survey and analyses. 
Each sampling occasion must have a unique identity number, ranging from 1: T, where T is the total number of sampling 
occasions. We note that each sampling occasion need not cover the entire survey area. 

Each camera trap location must be given a unique identification number, ranging from 1: J, where J is the total number 
of camera trap locations used in the survey. 

Assuming a camera trapping survey was conducted from 10 Jan 2009 to 30 Jan 2009, and, we treat each day as a sampling 
occasion, we end up with 20 sampling occasions. Let us assume there were 16 camera trap locations used in this survey and 
for logistical reasons the study area was partitioned into 4 blocks and each block contained 4 camera trap locations. 
Let us assume further that camera traps were deployed for 5 successive days in each block covering each of the four 
blocks successively in 20 days.  Assuming only 6 animals were photo-captured and identified, the \bold{INPUT FILE 1} for 
\pkg{SPACECAP} would look like:

\tabular{lll}{
LOC_ID    \tab  ANIMAL_ID \tab  SO  \cr
1    \tab  4 \tab  17  \cr
5    \tab  5 \tab  13  \cr
7    \tab  6 \tab  17  \cr
8    \tab  6 \tab  16  \cr
9    \tab  3 \tab  20  \cr
11   \tab  1 \tab  14  \cr
15   \tab  2 \tab  12  \cr
}

The first data row tells us that Animal ID no 4 was captured at Location ID 1 on the 17th sampling occasion. 
In your spreadsheet application, you can build up the animal capture details data as per the above format 
(please make sure that you INCLUDE the header row with the exact titles on column headings as shown above). 
The file must eventually be saved as a \emph{.csv} file (for example, "captures.csv") and saved in the working directory.

\bold{INPUT FILE 2: Trap Deployment Details} 

In many camera trap sample surveys of animals, all camera trap stations in the study area 
may not be operational simultaneously for logistical reasons (for example: limited number of cameras or manpower). 
Therefore, the trap deployment details input file provides \pkg{SPACECAP} with the information on the dates when each camera 
trap location was active and operational during the survey. Some trap-nights or trap-days of capture data may be "lost" 
as a result of camera trap failure, theft, vandalism or animal-damage. This type of trap activity/passivity information 
can also be effectively fed into and used in \pkg{SPACECAP}. The trap deployment details file records both these types of 
information, thus accurately accounting for trapping effort. 

The trap deployment data are stored in a two dimensional matrix of camera trap locations and sampling occasions in a 
binary, 0/1 format, where 0 indicates that a particular camera trap station was NOT operational on a particular sampling 
occasion, and 1 indicates that it was operational. The trap location is denoted in 3 columns in the table representing the 
Trap Location ID no., the spatial location expressed in X and Y-coordinates (in Universal Transverse Mercator UTM 
projection system in GIS). It is important that these coordinates are represented in the UTM projection system, because 
it is used for all distance measurements and computations in \pkg{SPACECAP}. 

We illustrate trap deployment data entry using a Table for the same example of animal capture details described above. 
Please recall that the camera trapping survey was conducted with a camera trap array of 16 trap locations, deployed in 
4 blocks with 4 camera trap sites in each block. The camera trap survey was conducted over 20 sampling occasions, during 
which each trapping block was sampled over 5 sampling occasions. The resulting TRAP DEPLOYMENT DATA file would look like 
the one below:   

\tabular{lllllllllllllllllllllll}{
LOC_ID \tab X_Coord \tab Y_Coord \tab 1 \tab 2 \tab 3 \tab 4 \tab 5 \tab 6 \tab 7 \tab 8 \tab 9 \tab 10 \tab 11 \tab 12 \tab 13 \tab 14 \tab 15 \tab 16 \tab 17 \tab 18 \tab 19 \tab 20 \cr
1 \tab 619303 \tab 1325966 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
2 \tab 624151 \tab 1325013 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
3 \tab 624722 \tab 1323864 \tab 1 \tab 1 \tab 1 \tab 0 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
4 \tab 621806 \tab 1322453 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
5 \tab 622451 \tab 1320137 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
6 \tab 622599 \tab 1317937 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
7 \tab 623179 \tab 1315941 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
8 \tab 625156 \tab 1315587 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
9 \tab 626022 \tab 1314224 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
10 \tab 627568 \tab 1315494 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
11 \tab 619604 \tab 1324739 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 0 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
12 \tab 621478 \tab 1324515 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \cr
13 \tab 623317 \tab 1323989 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \cr
14 \tab 624406 \tab 1321603 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \cr
15 \tab 624482 \tab 1320577 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \cr
16 \tab 629229 \tab 1319793 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \cr
}

The table shows that the camera trap sites (Loc 1-4) were operational during sampling occasions 1-5 but were not 
operational on the remaining sampling occasions 6-20. Additionally, an odd '0' corresponding to Loc 3 and sampling 
occasion 3 indicates that a trap night was "lost" here. 

The trap deployment details file should be constructed exactly as above in your spreadsheet application. 
It must then be converted to a comma separated ASCII file (.csv file) with an appropriate filename (e.g. traps .csv) 
in the working directory accessed by \pkg{SPACECAP}. 

\bold{INPUT FILE 3: Potential Home-Range Centers}

In \pkg{SPACECAP} analyses, the surveyed area containing the camera trap array combined with an extended area surrounding it, 
known as the "state-space" of the underlying point process, say S, which is represented by a large number of equally 
spaced points in the form of a very fine mesh. These points are visualized as representing all possible potential activity 
centers (or home range centers) of all the animals in the animal population being surveyed. This fine grid or mesh of 
points can be easily generated using a GIS software (like ArcView, MAPINFO etc.) as briefly described below. We view this 
as an approximation to an underlying continuous state-space which, in practice, would normally be difficult to 
characterize for computational purposes except in very basic situations where regular polygons might be reasonable. 
While estimates of population size, N, will be sensitive to the size and extent of the state-space, the estimated density 
D=N/||S|| is invariant as the extent of the state-space increases. Thus, S should be chosen to be sufficiently large so as 
to ensure stability of the density estimate. Conceptually, this occurs (under the models fitted in \pkg{SPACECAP}) by choosing S 
to buffer the trap array by 2 or 3 times the encounter probability scale parameter. 

First a rectangle is formed by connecting the outermost camera trap locations using the GIS software of your choice. This 
rectangle is called the "Minimum area rectangle". A buffer distance (which is sufficiently large to ensure that no 
individual animal outside of the buffered region has any probability of being photo-captured by the camera traps in the 
array during the survey) is added to the rectangle around encompassing the trap array. Thereafter, using GIS numerous 
equally spaced points representing home range centers are generated for this extended area. 
In practice, some of these potential home range centers in the mesh may end up in habitats known to be entirely unsuitable 
for the study species (say in the middle of a village, for tiger data). \pkg{SPACECAP} appropriately deals with this problem 
because "Grid Cells" input file clearly specifies which of these potential home range/activity centers lie within suitable 
habitat and which do not. 

The potential home-range centers data file essentially consists of a 3 column table. The first two columns are the X and Y 
coordinates (both in UTM projection system) of all the potential activity centers (the equally spaced points generated 
from your GIS software) and the third column is a habitat suitability indicator column, indicated with 1s if the potential 
activity centers lies within suitable species habitat or with 0s otherwise. The format of this file is as given below 
(this example is a small subset of the actual data set, purely for illustration):

Example of potential home-range centers input file format:

\tabular{lll}{
X_COORD    \tab  Y_COORD \tab  HABITAT  \cr
611734      \tab  1299581 \tab  1  \cr
612301.6    \tab  1299583 \tab  1  \cr
612869.2    \tab  1299585 \tab  0  \cr
613436.9    \tab  1299587 \tab  0  \cr
614004.5    \tab  1299589 \tab  0  \cr
614572.2    \tab  1299591 \tab  0  \cr
615139.8    \tab  1299593 \tab  0  \cr
615707.5    \tab  1299595 \tab  1  \cr
616275.1    \tab  1299597 \tab  1  \cr
616842.8    \tab  1299600 \tab  0  \cr
617410.4    \tab  1299602 \tab  1  \cr
617978.1    \tab  1299604 \tab  1  \cr
618545.7    \tab  1299606 \tab  1  \cr
}

Such a table must be created for all the potential activity centers (this will be a very large table) and saved as a 
comma delimited ASCII file (.csv file) with an appropriate name, for example, centers.csv, and saved in the working 
directory.

In the input data panel, there is a text box for specifying the area of potential home-range centers. This area may be 
imagined to be composed of a point at the centre of a square, which we call "pixel". Please enter the pixel size area 
used by the GIS software in this box. The "fineness of the mesh" determined by the pixel size used for spacing of 
potential home range centers is dictated by species biology (e.g. a few hundred meters for tigers but only a few meters 
for civets). Caution should be exercised here to not specify a state-space that is too fine because the MCMC algorithm 
run time increases linearly with the size of the state-space grid. Regarding this as a discrete approximation to some 
underlying continuous state-space probably justifies a reasonably coarse state-space grid.  Regardless of the desired 
state-space dimension, we recommend carrying out a trial run with a coarse grid to evaluate the performance.


Please load these files using on screen 'buttons' provided in the input data panel of \pkg{SPACECAP}. Please specify the pixel 
size area of a potential home-range center (as created in input file 3 above) in square kilometers and then click on "OK". 
Please check the frame at the bottom for status or error messages. In case you wish to edit your selection, please click 
on the "Edit" button and start the selection all over again. 

Proceed to STEP 2, the Model Definition frame.

\bold{STEP 2: SELECTING THE APPROPRIATE MODEL COMBINATION FOR ANALYSIS}

The Model Definition panel of \pkg{SPACECAP} consists of a set of options to select an appropriate model combination for the 
Spatial-Capture Recapture Analysis. These are simple radio buttons indicating each model choice. Some of these model 
options that are "grayed out" are expected to be made available in future developments of \pkg{SPACECAP}. 

The model choices are:

1.	\bold{Trap response present OR Trap response absent}
Selecting "Trap response present " option runs the behavioral response option (equivalent to Model "Mb"). Select the 
"Trap response absent" option if you decide otherwise.  This model implements a local or "trap-specific" behavioral 
response under which the probability of encounter in a trap increases subsequent to initial capture in that trap. This is 
in contrast to the conventional "global" behavioral response which parameterizes a constant increase in encounter 
probability (on the logit scale, usually) that is not trap specific. 

2.   \bold{Spatial Capture-Recapture OR Non-spatial Capture-Recapture}
Select "Spatial Capture-Recapture" for running a spatially explicit capture-recapture analysis, or 
"Non-spatial Capture-Recapture" for running a conventional capture-recapture analysis (this is equivalent to the Null 
Model "Mo" in non-spatial CR analysis)
   
3.	\bold{Half Normal OR Negative Exponential}
Currently \pkg{SPACECAP} analyses SECR models with only the Half-Normal detection function.

4.	\bold{Bernoulli (binary)  OR Poisson encounter process}
       Currently the analysis is run with the Bernoulli encounter model in which the probability of success is derived 
       as the probability of a positive response under a Poisson encounter rate model. This motivates use of the 
       complementary log-log link which relates encounter probability to distance and other covariates. After the model 
       definition is complete, please click on "OK". Please check the frame at the bottom for status or error messages. 
       In case you wish to edit your selections, please click on the "Edit" button, change your model definition and click 
       on "OK" again. Proceed to the MCMC simulation settings frame.

\bold{STEP 3: Setting the Markov-Chain Monte Carlo (MCMC) parameters (for advanced users)}

\pkg{SPACECAP}, uses the Markov-Chain Monte Carlo simulation algorithm written in Program \R (R Development Core Team ) 
to estimate the parameters of the Spatially-explicit Capture Recapture models of Royle et al 2009). The relevant settings 
can be set in the MCMC simulation settings panel of \pkg{SPACECAP}. 

\bold{No of iterations} - This defines the number of MCMC iterations for the analysis (if you aren't sure what this means please 
set this to a value of about 50,000)

\bold{Burn-in} - This defines the number of initial values to discard during the MCMC analysis. Setting this at about 1000 
usually works well in our experience.  Note that some evaluation of whether this is sufficient should be carried out 
using conventional methods (a topic we will address in subsequent releases).

\bold{Thinning} - This defines the thinning rate. Only iteration numbers defined by the thinning rate are stored during the 
analysis (if you aren't sure of what this means, you may set this at a value of 1 -that is, no thinning. If you are 
setting it to 1, please make sure that the number of iterations specified is less than 65,000. This is because Microsoft 
Excel, which in most versions has a maximum limit of 65000 rows, will not be able to open our output file).    

\bold{Data augmentation} - Since we are uncertain about the total number of animals, which is likely to be larger than the 
minimum number caught during your camera trapping survey, you will need to "augment" this value by a certain amount. 
Ideally, you would like this to be a very large number relative to the number you have caught, but setting it up to 
be very high will cause the analysis to run for a very long time. As a rule of thumb, you could set this to a value 
of about 5-10 times the number of animals you have photo-captured during the survey. Data augmentation (DA) is a 
computational device that enables a convenient Bayesian analysis of capture-recapture models where N is unknown. 
In the context of SECR models, N is the population of individuals having their activity centers on the prescribed 
state-space. The basic idea of DA is to provide an upper bound on N, say M, which is equal to N plus the number of 
augmented individuals.  Technically, M is the upper limit of a uniform(0,M) prior for N, which is a customary 
"noninformative" prior for N in this context. As a practical matter, data augmentation creates a list of 
pseudo-individuals that are always available for the MCMC algorithm to "use" if necessary. That is, these 
pseudo-individuals leave and enter the population depending on the current values of the model parameters.  
See Royle, Dorazio and Link (2007) for some general context and Royle et al. (2009) for details in the context of spatial 
capture-recapture models. 

After the MCMC simulation values have been specified, please click on "OK". Please check the frame at the bottom for 
status or error messages. In case you wish to edit these settings, please click on the "Edit" button, edit these values 
and click on "OK" again. You are now all set to start the analysis. 

\bold{STEP 4: Running the analysis}

The last step will simply involve activating the RUN option in the top menu bar. This will start performing the 
analysis and you will see a progress bar indicating the status of the analysis. Currently, an analysis involving 
50000 iterations takes about 14 hours on a fast computer.  We are working hard to make the algorithms run faster - 
we promise! 

\bold{Results}

The posterior density estimates along with standard errors appear as a table in the output panel upon the completion of 
the analysis. This table also reports estimates of parameters lam0, sigma, psi and beta. If the analysis was run with 
trap response present, the estimates of "p1" and "p2" are also reported.

Additionally, all the results are written into a comma separated file (called \cr param_val_<timestamp>.csv) and is saved 
into the current working directory. All the summary statistics are written into a file called 
summary_stats_<timestamp>.csv, which is also saved to the current directory. And the posterior density graphs of all 
parameters are all stored as jpeg files (.jpg) in the current working directory. 

\bold{Further Developments in progress}

We plan to add a number of additional utility modules for various applications, which we will be incorporating in the 
subsequent versions to be released in next few months. Some of the features we are currently working on are:

1. Include the Poisson encounter process model. \cr
2. Include the Negative-Exponential distance function model. \cr
3. Provide with the Hazard-Rate distance function model. \cr
4. Generate FIT statistics for model selection using the Bayesian p-value. \cr

\bold{Suggested Citation:}

Pallavi Singh, Arjun M. Gopalaswamy, Andrew J. Royle, N. Samba Kumar and K. Ullas Karanth  (2010). 
SPACECAP: A Program to Estimate Animal Abundance and Density using Bayesian Spatially-Explicit Capture-Recapture Models. 
Wildlife Conservation Society - India Program, Centre for Wildlife Studies, Bangalure, India. Version 1.0.1

}

\author{
Pallavi Singh, Arjun M. Gopalaswamy, Andrew J. Royle, N. Samba Kumar and K. Ullas Karanth 
with contributions from Sumanta Mukherjee, Vatsaraj & Dipti Bharadwaj
}

\references{

Borchers, D. L., and M. G. Efford. 2008. Spatially explicit maximum likelihood methods for capture-recapture studies. 
  \emph{Biometrics} 64:377-385.

Karanth, K. U., and J. D. Nichols. 1998. Estimation of tiger densities in India using photographic 
  captures and recaptures. \emph{Ecology} \bold{79}:2852-2862.

Parmenter, R., et al. 2003. Small-mammal density estimation: a field comparison of grid-based
  vs. web-based density estimators. \emph{Ecological Monographs} \bold{73}:1-26.

R Development Core Team.2005. R: A language and environment for statistical computing. R Foundation for Statistical 
  Computing,  Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org.

Royle, J. A., R. M. Dorazio, and W. A. Link. 2007. Analysis of multinomial models with unknown index using data 
  augmentation. \emph{Journal of Computational and Graphical Statistics} \bold{16}:67-85. 

Royle, J. A., K. U. Karanth, A. M. Gopalaswamy and N. S. Kumar. 2009. Bayesian inference in camera trapping studies for a class of spatial capture-recapture models. \emph{Ecology} \bold{90(11)}, 3233-3244.

Trolle, M., and M. Kery. 2003. Estimation of ocelot density in the Pantanal using capture-
            recapture analysis of camera trapping data. \emph{Journal of Mammalogy} \bold{84}:607-614.

Wilson, K. R., and D. R. Anderson. 1985a.  Evaluation of a nested grid approach for estimating density. 
  \emph{Journal of Wildlife Management} \bold{49}:675-678.

Wilson, K. R., and D. R. Anderson. 1985b.  Evaluation of two density estimators of small mammal population size. 
  \emph{Journal of Mammalogy} \bold{66}:13-21. 

}

\examples{
SPACECAP()
}
