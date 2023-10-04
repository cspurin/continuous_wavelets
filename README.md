# continuous_wavelets
Here we take the differential pressure from flow in porous media experiments, and use continuous wavelet transformation (CWT) to gain further insight into the pressure data. 

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
1. dt
2. mirror = True
3. om0 = 6
4. dj=0.1
5. wf = 'dog' # or 'morlet'
6. mirrormethod = 2 # 1 or 2 (mirroring method reduces end effects, so type of mirroring depends on the pressure signal and overall trend)
7. filter
8. infile

Update these for your data set and then you are good to go! 

# Guide to code outputs 
## The first output is the pressure data that is used in the CWT 
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/5e300e9b-edfd-4d25-a40d-58f58de95666)
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/696a0e3e-36bb-49fe-bd03-f4c4e7fc72d0)
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/8603a9af-fa22-4682-bcea-9ce7dac6f2d2)
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/831b4825-1d9b-400f-a293-ffc84eab23d2)
![image](https://github.com/cspurin/continuous_wavelets/assets/108369280/20455d8a-9607-4879-83be-031ea30474a3)
