The scripts in this folder are designed to work if 
(a) all of them are in the same folder. 
(b) R's working directory is set to that folder. 
(c) The folder has two folders called ARTR and PSSP that hold the
    necessary data files (supplied along with the scripts). 

If those are not all true, the scripts will crash. The 'main programs' are:

(1) CreateTestData.R does what it says: it creates a series of artificial data sets
with known competition kernel, based on the real ARTR and PSSP data sets in the
subfolders of the same name.  

(2) FitMonoDistanceSpline_SIM.R does a fit to an artificial data set. Comments in the
script explain how to choose which dataset is analyzed and assign a name to the folder
where the results get stored. 

(3) FitMonoDistanceSpline_REAL.R does a fit to a real data set. Again, comments in the
script explain how to choose the data set, and assign a name to the results folder.  

All three of these should run 'out of the box' and do an analysis like the ones
reported in the paper. However, the scripts use a very short list of 
smoothing parameter values so that the run-time is not too long. A real analysis
should use a longer list (e.g., log10 lambda = seq(-6,6,by=1) or even more).   

LoadData.R, DistanceKernelSubs.R, and DistanceSplineSubs.R are "utility" files
that define functions called by the three 'main programs'. The purpose of each 
function is explained be comments in each of the files. 

The analysis of real data with random year effects mentioned in the paper is
implemented by sourcing DistanceSplineSubs_T.R after you source
DistanceSplineSubs.R. See line 17 of FitMonoDistanceSpline_REAL.R, where 
you can change the value of doYear to make this happen. This is 
'experimental' -- it should work, but it has not been thoroughly tested 
on artificial data. 

