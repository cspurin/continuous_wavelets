# continuous_wavelets
Here we take the differential pressure from flow in porous media experiments, and use continuous wavelet transformation (CWT) to gain further insight into the pressure data. The fluctuations in the pressure data provide key insights to the underlying pore-scale dynamics. 

The paper associated with work is: https://doi.org/10.1029/2023GL104473

# Guide to getting the code running 
Recommended steps for running the notebook.

1. Download and install vscode: https://code.visualstudio.com/
2. Download and install python: https://www.python.org/downloads/ (download version 3.10.11)
3. Install the python extension in vscode. In vscode, click the extensions icon in the sidebar and search for "python"
4. Open a new folder (link it the folder containing the contents of this repository, wherever that is on your computer)
5. Create a virtual enviornment. In vscode click View>Terminal and enter the code below: <br> 
   python -m venv myenv <br> 
   Then for Mac enter: source myenv/bin/activate <br> 
   or for Windows: myenv\Scripts\activate 
7. Install the necessary packages by typing this in the terminal:
   pip install -r requirements.txt
   NB must be in the folder where the requirements.txt file is. Alternatively, you can install the packages individually e.g. pip install numpy.
   Now open the notebook (wavelet_plotting_GRL_paper.ipynb) and begin!

# Guide to code inputs 
The variables for the wavelet analysis are 
1. dt. This is the time step between pressure data points in seconds. 
2. mirror = True. This mirrors the pressure data to minimise any edge effects. Can also be set to 'False".
3. wf = 'dog' # or 'morlet'. This is the wavelet. Options are 'dog', which is the Derivative-of-Gaussian wavelet. The Morlet wavelet can also be selected.
4. om0 = 6. This is the derivative of the wavelet. 
5. dj=0.1. This is the scale step. The power is rectified by scale as done in Liu et al. (2007). See manuscript for more details. 
6. mirrormethod = 2 # 1 or 2 (mirroring method reduces end effects, so type of mirroring depends on the pressure signal and overall trend). Mirror method 1 mirrors the data in the x axes only, mirror method 2 mirrors the data in the x and y axes. 
7. filter. If you have a lot of data, you can filter it to make the code run faster. Filter can = 'None' for no filtering. Filter = 1 for a moving average filter, and filter = 2 for downsampling of the data set.
8. infile. This is the file where your pressure data is saved. It needs to be a .dat file containing just the pressure measurements (code assumes input of kPa). 

Update these for your data set and then you are good to go! 

# Guide to code outputs 
## The first output is the pressure data that is used in the CWT 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/5e300e9b-edfd-4d25-a40d-58f58de95666)

## The second output is the input pressure data, and the pressure data used in the CWT
Note, if the filter is set to 'None' then these pressure plots will be the same. The mirrored pressure signal is also given, to see how well edge effects will be managed in the continuous wavelet transformation. 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/696a0e3e-36bb-49fe-bd03-f4c4e7fc72d0)

## The third output if the rectified power spectra 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/8603a9af-fa22-4682-bcea-9ce7dac6f2d2)

## The fourth output is the amplitude of fluctuations versus their frequency
This is the equivalent of a Fourier analysis plot. Red noise is indicated on the plot as this type of noise is common for flow in porous media. See https://doi.org/10.1029/2022WR031947 for further details on this. 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/831b4825-1d9b-400f-a293-ffc84eab23d2)

## The fifth output is the amplitude of fluctuations versus their frequency for different sections of the signal
The sections are shown in the pressure plot and can be adjusted within the code. 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/20455d8a-9607-4879-83be-031ea30474a3)
