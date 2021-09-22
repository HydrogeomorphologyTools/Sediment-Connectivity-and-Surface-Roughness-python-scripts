# Sediment-Connectivity-and-Surface-Roughness-python-scripts <br/>
Python scripts that can be used singularly or integrated into batch/workflow processing for computing Sediment Connectivity as in Cavalli et al. 2013 and Surface Roughness as in Cavalli and Marchi 2008 <br/>

Welcome to the **Sediment-Connectivity-and-Surface-Roughness-python-scripts** wiki!

**Authors:** Stefano Crema, Alessandro Sarretta, Giorgia Macchi and Marco Cavalli - CNR IRPI 2021

***

# DEPENDANCIES

In order to run, the provided functions **need to have TauDEM installed**, you can download and install the latest complete windows installer from [https://hydrology.usu.edu/taudem/taudem5/downloads.html](https://hydrology.usu.edu/taudem/taudem5/downloads.html)

For Linux users, TauDEM is provided also as C++ source code and a compilation can be tried (with some manual tuning) bearing in mind that it makes use of MPI, for windows users the complete windows installer will take care of installing the dependencies (GDAL, MPI/HPCPack, Visual Studio libraries).<br/>

**Mandatory Python imports are** <br/>

osgeo --> gdal<br/>
os<br/>
sys<br/>
time<br/>
numpy<br/>
scipy --> signal<br/>

***

# DESCRIPTION

Two main functions are made available within the repository:

1. * _**cavalli_roughness.py**_ computes surface roughness as in Cavalli and Marchi, (2008) i.e. as the standard deviation of residual topography
The function needs an input DTM and a parameter indicating the desired moving window extent and provides as output a raster of the surface roughness together with a raster of the roughness-derived weighting factor to be used for sediment connectivity calculation as reported in Cavalli et al., (2013).
Roughness is derived with a robust approach by first detrending the original surface (original-averaged surface over desired distance) and then looking at the detrended variability as a measure of surface roughness. Robustness, implications, and applications of such an approach are reported also in Crema et al., (2020)

2. * _**SedInConnect_target.py**_ calculates the index of sediment connectivity as described in Cavalli et al., (2013) and Crema and Cavalli, (2018).
The function needs an input pit-filled DTM, a parameter indicating the cell size, an input weighting factor (that can be derived with the _cavalli_roughness_ function), and a target shapefile that we encourage to be polygonal (if it's a line or polyline we strongly encourage to make a 1-cell buffer of it and use such a polygon as target shapefile so as to force all drainage lines intersecting the desired target). The function provides as output the index of connectivity in respect to the desired target. It also clips automatically the extent to the basin draining into the target only.

***

# USAGE

Both functions can be modified and tuned based on the user's needs. As provided, the two functions work in the same way that is to say they accept at the beginning of the script the path to the input files and parameters.\
After the function definition, at the end of the file, the function itself is called with the input/output files and parameters set at the beginning and produces the desired outcomes.\
This way it can be used for a single run, but it can be also called in a loop (or thread) to run in batch mode or integrated into an existing workflow.


***

# **References**

* Cavalli, M., Marchi, L., 2008. Characterisation of the surface morphology of an alpine alluvial fan using airborne LiDAR. Nat Hazards Earth Syst Sci 8, 323–333. https://doi.org/10.5194/nhess-8-323-2008
* Cavalli, M., Trevisani, S., Comiti, F., Marchi, L., 2013. Geomorphometric assessment of spatial sediment connectivity in small Alpine catchments. Geomorphology, Sediment sources, source-to-sink fluxes and sedimentary budgets 188, 31–41. https://doi.org/10.1016/j.geomorph.2012.05.007
* Crema, S., Cavalli, M., 2018. SedInConnect: a stand-alone, free and open source tool for the assessment of sediment connectivity. Comput. Geosci. 111, 39–45. https://doi.org/10.1016/j.cageo.2017.10.009
* Crema, S., Llena, M., Calsamiglia, A., Estrany, J., Marchi, L., Vericat, D., Cavalli, M., 2020. Can inpainting improve digital terrain analysis? Comparing techniques for void filling, surface reconstruction and geomorphometric analyses. Earth Surf. Process. Landf. 45, 736–755. https://doi.org/10.1002/esp.4739


# USAGE EXAMPLE _**cavalli_roughness.py**_
##################################################################################### <br/>
  ----------------------------------USAGE------------------------------------- <br/>
############# Simply set the following input and output files and parameters (examples provided) ######### <br/>
############# Remember to use forward slash in paths ######################################### <br/>

dtm_r_txt = 'D:/Research/SedInConnect_python/dtmfel.tif'  # input DTM <br/>
f_s_txt = 5  # input moving window size in pixel units <br/>
ri_out_txt = 'D:/Research/SedInConnect_python/ri.tif'  # output roughness index <br/>
w_out_txt = 'D:/Research/SedInConnect_python/w.tif'  # output <br/>

##################################################################################### <br/>

And the function runs with such input/output settings by calling internally: rw_cavalli(dtm_r_txt, f_s_txt, ri_out_txt, w_out_txt)<br/>

You can use this piece of code for single or batch processing and/or in combination within a programming workflow. <br/>
