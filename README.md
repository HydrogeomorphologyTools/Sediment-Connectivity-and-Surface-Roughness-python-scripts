# Sediment-Connectivity-and-Surface-Roughness-python-scripts <br/>
Python scripts that can be used singularly or integrated into batch/workflow processing for computing Sediment Connectivity as in Cavalli et al. 2013 and Surface Roughness as in Cavalli and Marchi 2008 <br/>

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

Necessary imports are:<br/>
osgeo --> gdal<br/>
os<br/>
sys<br/>
time<br/>
numpy<br/>
scipy --> signal<br/>
