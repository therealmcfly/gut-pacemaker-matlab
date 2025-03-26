# Detection Algorithm Model

This folder encapsulates the model that has been developed for the duration of the project. It is composed of many files, 
with the key script being: **detection_script**, that calls on many other functions to model functionality of the detection 
algorithm.

## List of files
* NEO_transform.m - mathematical function that operates on a signal - used for magnitude change detection
* artifact_detect.m - function to detect pacing artifacts
* artifact_remove.m - function to remove the detected artifacts
* detection_function.m - combined detection algorithm functionality encapsulated within a function - provided input parameters (used mainly in the application)
* detection_script.m - overall combined detection algorithm model  
* exp_10_output.mat - low-resolution experimental data
* exp_16_output.mat - low-resolution experimental data
* fir_51.mat - high-pass filter coefficients
* importfile.m - function to import relevant files  
* moving_average_1s_window.m - moving average filter function
* pig37exp2.mat - high-resolution experimental data
* pig41exp2.mat - high-resolution experimental data
* pigApp.mlapp - makeshift application to provide a more interactive demonstration of the algorithm
* pig_img.png - image for the app0
* visualise_signal.m - function to plot signals

## Running the model
No setup is required - just ensure that the workspace is assigned correctly.

The main model is run through the script **detection_script**. This can be run as any normal MATLAB script, and 
parameters can be changed as desired.

### A list of changeable parameters are as follows:

#### To test "Low-Resolution" signals:
* changing between exp16 .mat and exp10 .mat - line 16
* channel selection (1 to 256) - line 18

#### To test "High-Resolution" signals:
* changing between pig41 and pig37 - line 27
* channel selection (1 to 40 for pig41 and 1 to 30 for pig37)

#### Other functionality changes
* artifact detection threshold - line 33
* buffer size - line 34
* any lines containing **visualise signal** for plotting and viewing purposes

*Note that you must comment out the "Bad signal" block if intending to use the "Good signal" block and vice versa.

# Application
A makeshift application was created for the exhibition day to have a more interactive demonstration of the 
detection algorithm. This is found in **pigApp.m** and is run like any normal MATLAB app. 

*Note that the application is the sole user of the **detection_function.m** file which is a decomposed version of **detection_script.m**

From here you can select the desired experiment, and iterate over the window just as the function does.
It also allows to change the artifact detection threshold.
