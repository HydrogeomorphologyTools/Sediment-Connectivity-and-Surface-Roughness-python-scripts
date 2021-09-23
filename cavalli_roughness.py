## Imports
## [suggested --> pip install or, in case of Windows OSs --> using unofficial windows binaries @ https://www.lfd.uci.edu/~gohlke/pythonlibs/
## unofficial binaries are optimal for gdal for example]

import osgeo
from osgeo import gdal
import os
import sys
import time
import numpy
from scipy import signal

"""
This script computes Surface Roughness as expressed in Cavalli et al. (2008)
For the Guidelines on the application usage please visit the groupware session of www.sedalp.eu
authors: Stefano Crema, Giorgia Macchi, Alessandro Sarretta & Marco Cavalli,
last edited: September 2021
Copyright (C) 2014-2021  Stefano Crema, Alessandro Sarretta, Giorgia Macchi & Marco Cavalli, @ CNR-IRPI, corso stati uniti, 4, 35127, Padova (Italy)
mail: stefano.crema@irpi.cnr.it; marco.cavalli@irpi.cnr.it

This software has been developed in the framework of the INTERREG V-A Italiy-Austria 2014-2020 ITAT3032 SedInOut Project


###############################################################################
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
For a copy of the GNU-GPL v2 license please visit: http://www.gnu.org/licenses/gpl-2.0.txt
"""

"""
####################################################################################
----------------------------------    USAGE    -------------------------------------

############# Here you set input and output files and parameters ###################
############# Remember to use forward slash in paths ###############################
"""
dtm_r_txt = 'D:/Research/SedInConnect_python/dtmfel.tif'  # input DTM
f_s_txt = 5  # input moving window size
ri_out_txt = 'D:/Research/SedInConnect_python/ri.tif'  # output roughness index
w_out_txt = 'D:/Research/SedInConnect_python/w.tif'  # output
"""
#####################################################################################
-------------------------------------------------------------------------------------

#####################################################################################
"""

"""
--------------------------------------
###########################################################################        

CAVALLI ROUGHNESS FUNCTION

###########################################################################
--------------------------------------
"""


# noinspection SpellCheckingInspection
def rw_cavalli(dtm_r, f_s, ri_out, w_out):
    start = time.time()  # take the time
    filename = str(dtm_r)
    tif = osgeo.gdal.Open(filename)  # opening the file (1 band, bit not specified, but usually 32)
    #
    # folder path
    dir_pat = os.path.dirname(os.path.realpath(filename))  # path of the selected input file
    dir_path = dir_pat.encode('ascii', 'ignore')
    #
    #
    # Opening Messages
    if tif is None:
        print(time.strftime("%d/%m/%Y %H:%M:%S    "),
              "couldn't open input DTM dataset, for Weighting factor computatin")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%Y %H:%M:%S    "), "opening DTM for weighting factor computation was successful!")
    #
    # cols, rows, bands
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands = tif.RasterCount
    #
    # dtm as an array                #################################       PAY ATTENTION!     ##########################
    i_ar = tif.ReadAsArray()
    i_ar = i_ar.astype(float)
    band = tif.GetRasterBand(1)  # bands retrieving
    geoinf = tif.GetGeoTransform()  # upper left coordinates, cellsize
    proj = tif.GetProjection()
    tipo = i_ar.dtype
    del tif
    #
    # mask_NoData
    ndv = numpy.min(i_ar)
    i_ar[i_ar == ndv] = numpy.nan
    #
    # kernel size
    size_filter = int(f_s)
    #
    #
    i_ar_p = numpy.ones((rows, cols),
                        dtype=tipo)  # matrice di numeri uno e zero sui nodata da usare per far la media poi
    i_ar_p[numpy.isnan(i_ar) == 1] = 0
    i_ar_d = numpy.ones((rows, cols), dtype=tipo)  # build up array for summing correctly elevation data later
    i_ar_d[numpy.isnan(i_ar) == 1] = 0
    i_ar_d[numpy.isnan(i_ar) == 0] = i_ar[numpy.isnan(i_ar) == 0]
    #
    #
    ker = numpy.ones((size_filter, size_filter))  # creating weights for sum and then an averaging filter
    # C = scipy.signal.convolve2d(I_ar,ker, 'same')#borders kept the same but not realistic values, it divides by n^2 even if it finds only 2 elements scipy version
    d_num_el = signal.convolve2d(i_ar_p, ker, 'same')  # i get the number of elements per window
    d_num_el[numpy.isnan(i_ar) == 1] = numpy.nan

    e_sum = signal.convolve2d(i_ar_d, ker, 'same')  # sum of values in the DTM
    e_sum[numpy.isnan(i_ar) == 1] = numpy.nan

    m = e_sum / d_num_el  # average =sum/number of elements
    dtm_r = numpy.float64(i_ar) - numpy.float64(
        m)  # residual DTM, now I need to compute the Std dev of this averaged DTM (so I need again mean values of this and then squared values to fasten the computation)
    #
    """Use float64 to have more precision in decimals"""
    #
    dtm_r[numpy.isnan(i_ar) == 1] = 0
    # mean of residual topography
    e_r = signal.convolve2d(dtm_r, ker, 'same')  # sum of values in the DTM
    e_r[numpy.isnan(i_ar) == 1] = numpy.nan
    m_r = e_r / d_num_el  # average =sum/number of elements
    # Now playing for standard deviation following: http://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
    #
    # DTM_SQ = numpy.float64(I_ar_d) * numpy.float64(I_ar_d) #build up array for summing correctly elevation data later, this way I squared the values otherwise I loose the decimal precision if array is starting with 4 digits+ decimal + sign
    #
    dtm_r_sq = numpy.float64(dtm_r) * numpy.float64(dtm_r)  # residuals squared
    #
    m_r_sq = numpy.power(m_r, 2)
    e_r_sq = signal.convolve2d(dtm_r_sq, ker, 'same')  # sum of values in the DTM
    e_r_sq[numpy.isnan(i_ar) == 1] = numpy.nan
    #
    m_r_sq = e_r_sq / d_num_el  # E[X^2] average of squared residual DTM values
    mrsq = numpy.power(m_r, 2)  # square of average residual DTM values
    ri = numpy.sqrt(m_r_sq - mrsq)  # E[X^2] - E[X]^2       WC = Weight Cavalli
    #
    """Pay attention to the number we're treating, if we need to compute std dev via convolve and squared trick, for a DTM we need to use numpy.float64(DTM), otherwise the number of digits is not enough"""

    del (i_ar_p, i_ar_d, d_num_el, ker, e_sum, m, dtm_r, e_r, m_r, dtm_r_sq, m_r_sq, e_r_sq, mrsq)
    #
    #
    sur_rough = gdal.GetDriverByName('GTiff').Create(str(ri_out), i_ar.shape[1], i_ar.shape[0], 1,
                                                     gdal.GDT_Float32)  # shape sono righe e colonne, 1 è il numero di bande
    sur_rough.SetGeoTransform(geoinf)  # conferisco coordinate e proiezione del raster in input
    sur_rough.SetProjection(proj)
    sur_rough.GetRasterBand(1).WriteArray(ri, 0,
                                          0)  # scrivo effettivamente il raster, prima avevo solo allocato la memoria
    sur_rough.GetRasterBand(1).GetStatistics(0, 1)  # calculate statistics for visualization
    sur_rough = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del sur_rough
    # mean_t[y,x]=numpy.average(isnan([I_ar_t[y,x], I_ar_t[y-1,x], I_ar_t[y+1,x], I_ar_t[y,x-1], I_ar_t[y,x+1], I_ar_t[y-1,x-1], I_ar_t[y+1,x+1], I_ar_t[y+1,x-1], I_ar_t[y-1,x+1]])==False)
    #
    ## compute the weighting factor as 1-(R/(R+0.001))
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), 'computing weighting factor ... ...')
    #
    maxr = float(ri[numpy.isnan(ri) == False].max())  # max element in surface roughness
    # Mr = maxr + 0.001#increase a bit, in this release, old style of standardization

    weig_fac = 1.0 - (ri / maxr)  # 1-... to avoid 1-1=0 at denominator in downslope component
    weig_fac[numpy.where(
        weig_fac < 0.001)] = 0.001  # all values in the range 0-0-001 are shifted to 0.001, the rest remains NoData
    #
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), 'Weighting factor calculated!')
    #
    wf = gdal.GetDriverByName('GTiff').Create(str(w_out), i_ar.shape[1], i_ar.shape[0], 1,
                                              gdal.GDT_Float32)  # shape sono righe e colonne, 1 è il numero di bande
    wf.SetGeoTransform(geoinf)  # conferisco coordinate e proiezione del raster in input
    wf.SetProjection(proj)
    wf.GetRasterBand(1).WriteArray(weig_fac, 0,
                                   0)  # scrivo effettivamente il raster, prima avevo solo allocato la memoria
    wf.GetRasterBand(1).GetStatistics(0, 1)  # calculate statistics for visualization
    wf = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del wf
    del i_ar
    del ri
    del weig_fac
    print(time.strftime("%d/%m/%Y %H:%M:%S    "), 'RI and W files saved!')
    tot_time = time.time() - start
    print('Execution Time = ', str(tot_time), ' seconds \nCongratulations and enjoy!')


## Execution
rw_cavalli(dtm_r_txt, f_s_txt, ri_out_txt, w_out_txt)
#
## Finish line!!!!
