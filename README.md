# CalDyn: a tool for an automatic calcium dynamics characterization

We give here a quick description on how to use the software, for further details on how everything are calculated, please refer to:
"CalDyn: An advanced high-throughput automated algorithm for monitoring intracellular calcium dynamics and electrophysiological signals in patient-derived cells".... to be cited


A Graphical User Interface (gui) is launched running the file Principal_V6_HA_*.m in matlab.


There exist different versions of this file which has been finely tuned for different signals. But any of them can be launched and the user is invited to create its own starting fileonce conditions have been chosen for a given experiment to make the analysis even faster.

In this file, the user will have to set the path to the excel file to be analyzed (folder variable). All results will be saved in the very same folder, in subfolder having the same name as the excel file.

See below for a list of the parameters that can be set in the gui either at start or also directly using the interface.

The GUI should look like:
![Gui_Figure](./Figures/Figure_GUI.png)


In [1], the signal is plotted in light grey, and in blue, the filtered signal is represented. This latter signal is the  one on which everything will be calculated, so it should be faithful to the original signal but with less noise, as depicted here (The filtering could be adjusted in box [2], the higher the value, the higher the filtering).
In [1], all the results from the peak detection algorithms are represented. In black, one will find the maxima and minima associated with each peak. In cyan, is the starting point of the peak and in green its finishing point in magenta are points taking along the decrease in signal corresponding to 0.8, 0.7 0.5 0.3 and 0.05 times the height of the initial signal.

The detection of the peak is made using the derivative of the signal searching for a couple of minimum/maximum above a threshold (corresponding to the rise and fall associated to each peak). The derivative of the signal is represented in red of box [3], the signal undergoes a new filtering step which can be adjusted using the smooth_length box [4].
The most important parameter to be set is the threshold  in box [5] (value between 0 and 1), it will set the dotted line represented in [3], above which the derivative maxima/minima will be searched. If this parameter is set to low, we will mainly detect the noise in the signal, if set to high, we will miss some peaks.


There is a possibility to manually remove peaks using [7] and clicking near the peak which needs to be erased.

It is also possible to remove the base line using [6] (The base line can be fitted either using a constant or a polynomial of degree 1 or 2). See below the result on the 
initial signal:
![remove_baseline](./Figures/remove_base_line.png)
If the box Ref=base line is checked, all the minima of the signals will be set to zero.

There is also an option for peaks classification (the panel corresponding to this option can be made apparent clicking on button [9]). Three thresholds will need to be set to classify the peaks between small medium and high, and also to  detect multipeaks (peaks which will not go down to the base line before rising again). If the classification option is set to 'on', then in the excel file where the results ares saved, there will be now 

Saved Parameters can be set using button [8], a new checkbox menu will appear on which the user will be able to choose the parameters he will want to save:
![parameters_to_be_saved](./Figures/parameters.png)


Parameters for the analysis ca








