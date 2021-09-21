#!/usr/bin/python
# -*- coding: utf-8 -*-

import osgeo
from osgeo import gdal
from osgeo import ogr
import os
import sys
import time
import numpy
import math

numpy.seterr(invalid='ignore')  # Caution, setting the error management in numpy

"""
This script computes the index of connectivity as expressed in Cavalli et al. (2013) and Crema and cavalli (2018)

authors: Stefano Crema, Giorgia Macchi, Alessandro Sarretta & Marco Cavalli,
last edited: September 2021
Copyright (C) 2014-2021  Stefano Crema, Alessandro Sarretta, Giorgia Macchi & Marco Cavalli, @ CNR-IRPI, corso stati uniti, 4, 35127, Padova (Italy)
mail: stefano.crema@irpi.cnr.it; marco.cavalli@irpi.cnr.it

This software has been developed in the framework of the Interreg SedInOut Project

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
cell_s_txt = 2.5  # input raster cell size
tgt_file_txt = 'D:/Research/SedInConnect_python/target.shp'  # polygonal target shapefile
w_txt = 'D:/Research/SedInConnect_python/w.tif'  # weighting factor raster
ic_out_txt = 'D:/Research/SedInConnect_python/ic_tg.tif'  # output

"""
--------------------------------------
###########################################################################        

CAVALLI INDEX OF CONNECTIVITY IN RESPECT TO SELECTED TARGET

###########################################################################
--------------------------------------
"""

"""
--------------------------------------
###########################################################################       

CONNECTIVITY TO TARGETS

###########################################################################
--------------------------------------
"""


def CavalliConnectivitytg(dtm_f, c_s, tg_f, w_f, out_ic_tg):
    #
    start_0 = time.time()
    #
    filename = dtm_f.replace('\\', '/')
    filename = str(filename)
    tif = osgeo.gdal.Open(filename)  # opening the file (1 band, bit not specified, but usually 32)
    #
    # folder path
    dir_path = os.path.dirname(os.path.realpath(filename))  # path of the selected input file
    dir_path = dir_path.replace('\\', '/')
    #
    # Opening Messages
    if tif is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input DTM dataset (Connectivity to targets)")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening DTM (Connectivity to targets) was successful!")
    #
    # cols, rows, bands
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands = tif.RasterCount
    #
    # dtm as an array                #################################       ATTENTION     ##########################
    tif_ar = tif.ReadAsArray()  # con GDAL old 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
    tif_ar = tif_ar.astype(float)
    band = tif.GetRasterBand(1)  # bands retrieving
    geoinf = tif.GetGeoTransform()  # upper left coordinates, cellsize
    proj = tif.GetProjection()
    tipo = tif_ar.dtype
    del tif
    #
    # cell_size
    cell_s = numpy.float32(str(c_s))
    #
    # create constant array to transform into raster
    const_ar = (tif_ar >= 0) * cell_s
    c1 = numpy.where(const_ar == 0)
    const_ar[c1] = -1  # Replace -1 for nodata
    #
    # target file selection
    #
    no_data_value = -9999
    file_tar = tg_f.replace('\\', '/')
    file_tar = str(file_tar)
    temp_tar = gdal.GetDriverByName('GTiff').Create((dir_path + "/targets.tif"), const_ar.shape[1], const_ar.shape[0],
                                                    1, gdal.GDT_Float32)
    temp_tar.SetGeoTransform(geoinf)  # Set coordinates and projection of input raster
    temp_tar.SetProjection(proj)
    band = temp_tar.GetRasterBand(1)
    band.SetNoDataValue(no_data_value)

    targ_ds = ogr.Open(file_tar)
    targ_layer = targ_ds.GetLayer()

    gdal.RasterizeLayer(temp_tar, [1], targ_layer, None, None, burn_values=[10])  # ,options= ['ALL_TOUCHED=TRUE'])
    temp_tar = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del temp_tar

    # tg=poly2ras.rasterize(file_tar, dir_path+"/targets.tif", prototype=filename, burn_values=[10], options=["MERGE_ALG=ADD"])#assign value to variable or arbitrary one
    tg = osgeo.gdal.Open(dir_path + "/targets.tif")
    target_ar = tg.GetRasterBand(1).ReadAsArray()
    target_ar = target_ar.astype(float)
    # target_ar[numpy.where(tif_ar == tif_ar.min())] = numpy.NaN #where input has nodata, also the target must have, useful in case of sinks that clip the target too, no need for connectivity to the outlet
    del tg
    #
    #
    # create constant array to transform into Fl length raster
    fl_ar = (tif_ar >= 0) * 1  # build constant raster based on true/false array
    c1 = numpy.where(fl_ar == 0)  # 1 data, -1 NoData
    fl_ar[c1] = -1
    #
    # d8 flow directions & outputs
    # os.system("mpiexec -n 8 D8Flowdir -p "+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#windows case
    os.system(("mpiexec -n 8 D8Flowdir -p ".lower()) + filename[0:-4] + "p.tif" + " -sd8 " + filename[
                                                                                             0:-4] + "sd8.tif -fel " + filename)  # unix case
    # Dinf flow directions & outputs
    # os.system("mpiexec -n 8 DinfFlowdir -ang "+filename[0:-4]+"angt.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#windows case
    os.system(("mpiexec -n 8 DinfFlowdir -ang ".lower()) + filename[0:-4] + "angt.tif" + " -slp " + filename[
                                                                                                    0:-4] + "slp.tif -fel " + filename)  # unix case
    # removing unnecessary Dinf slope tif file
    os.remove(filename[0:-4] + "slp.tif")
    #
    """
    ---------------------------------------------
    WORKING ON D_DOWN COMPONENT
    ---------------------------------------------
    """
    # opening d8 flow directions
    tif_fdir8 = osgeo.gdal.Open(
        filename[0:-4] + "p.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_fdir8 is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input fdir_8 dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening fdir_8 was successful!")
    # Array
    tif_fdir8_ar = tif_fdir8.ReadAsArray()
    tif_fdir8_ar = tif_fdir8_ar.astype(numpy.float)
    ndv = numpy.min(tif_fdir8_ar)
    tif_fdir8_ar[tif_fdir8_ar == ndv] = numpy.NaN
    del tif_fdir8
    # Now_setting_NoData_where I have targets==10
    tif_fdir8_ar[(target_ar == 10)] = -1000  # I will then consider all negative values as NoData
    #
    # opening Dinf flow directions
    tif_fdirinf = osgeo.gdal.Open(
        filename[0:-4] + "angt.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_fdirinf is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input dirinf_tg dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening dirinf_tg was successful!")
    # Array
    tif_fdirinf_ar = tif_fdirinf.ReadAsArray()
    ndv = numpy.min(tif_fdirinf_ar)
    tif_fdirinf_ar[tif_fdirinf_ar == ndv] = -9999
    del tif_fdirinf  # otherwise it is in use and will not be deleted
    # Now_setting_NoData_where I have targets==10
    tif_fdirinf_ar[(target_ar == 10)] = -1000  # I will then consider all negative values as NoData
    # removing  Dinf  tif file, and rewrite it, the one with NoData
    # os.remove(filename[0:-4]+"angt.tif")
    # creating the new Dinf
    fd_inf_ds = gdal.GetDriverByName('GTiff').Create((filename[0:-4] + 'ang.tif'), const_ar.shape[1], const_ar.shape[0],
                                                     1, gdal.GDT_Float32)
    fd_inf_ds.SetGeoTransform(geoinf)
    fd_inf_ds.SetProjection(proj)
    fd_inf_ds.GetRasterBand(1).SetNoDataValue(-9999)
    fd_inf_ds.GetRasterBand(1).WriteArray(tif_fdirinf_ar, 0,
                                          0)
    fd_inf_ds = None
    del fd_inf_ds
    #
    # Slope d8 file modifications
    tif_sd8 = osgeo.gdal.Open(
        filename[0:-4] + "sd8.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_sd8 is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input sd8 dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening sd8 was successful!")
    # Array
    tif_sd8_ar = tif_sd8.ReadAsArray()
    del tif_sd8
    #
    # imposing upper and lower limits to slope, no data here are -1
    tif_sd8_ar[(tif_sd8_ar >= 0) & (tif_sd8_ar < 0.005)] = 0.005
    tif_sd8_ar[(tif_sd8_ar > 1)] = 1
    tif_sd8_ar[(tif_sd8_ar < 0)] = -1
    #
    # Create a reclassified slope d8 tiff array from numpy using gdal(better after removing it if it exists)
    dst_d8sl_ds = gdal.GetDriverByName('GTiff').Create((filename[0:-4] + 's.tif'), const_ar.shape[1], const_ar.shape[0],
                                                       1,
                                                       gdal.GDT_Float32)  # raster shape is rows and cols order, 1 is number of bands
    dst_d8sl_ds.SetGeoTransform(geoinf)
    dst_d8sl_ds.SetProjection(proj)
    dst_d8sl_ds.GetRasterBand(1).WriteArray(tif_sd8_ar, 0, 0)  # write the raster, previously only memory was allocated
    dst_d8sl_ds = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_d8sl_ds
    #
    # Selecting weighting factor
    filewgt = w_f.replace('\\', '/')
    filewgt = str(filewgt)
    tif_wgt = osgeo.gdal.Open(filewgt)  # opening the file (1 band, bit not specified, but usually 32)
    #
    # Opening Messages
    if tif_wgt is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input Weight dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening Weight was successful!")

    # wgt factor matrix
    tif_wgt_ar = tif_wgt.ReadAsArray()
    del tif_wgt
    #
    # Computing 1/(W*S)
    ws_1 = 1 / (tif_wgt_ar * tif_sd8_ar)
    #
    #
    # Calculating flow accumulation in order to retrieve after the coordinates of the outlet
    # os.system("mpiexec -n 8 AreaD8 -p "+filename[0:-4]+"p.tif"+ " -ad8 "+filename[0:-4]+"ad8.tif -nc")#windows case
    os.system(("mpiexec -n 8 AreaD8 -p ".lower()) + filename[0:-4] + "p.tif" + " -ad8 " + filename[
                                                                                          0:-4] + "ad8.tif -nc")  # unix case

    #
    # Trying to do a weighted flow length, USING FLOW DIR MATRIX, that removes the border cells
    # adding top row and left column with zeros
    #
    start = time.time()  # for computational time
    #
    # zero matrix bigger than F_dir8, to avoid border indexing problems
    fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0]) + 1, (tif_fdir8_ar.shape[1]) + 1), dtype=numpy.float32)
    fd8[1:fd8.shape[0], 1:fd8.shape[1]] = fd8[1:fd8.shape[0], 1:fd8.shape[1]] + tif_fdir8_ar
    # adding bottom row and right y axis with zeros
    fdir8 = numpy.zeros(shape=((fd8.shape[0]) + 1, (fd8.shape[1]) + 1), dtype=numpy.int16)
    fdir8[:fdir8.shape[0] - 1, :fdir8.shape[1] - 1] = fd8
    #
    # zero matrix bigger than weight, to avoid border indexing problems, and have same indexing as Fdir8
    wg = numpy.zeros(shape=((tif_wgt_ar.shape[0]) + 1, (tif_wgt_ar.shape[1]) + 1), dtype=numpy.float32)
    wg[1:wg.shape[0], 1:wg.shape[1]] = wg[1:fd8.shape[0], 1:wg.shape[1]] + ws_1  # the weight to weight tha flow length
    # adding bottom row and right y axis with zeros
    wgt = numpy.zeros(shape=((wg.shape[0]) + 1, (wg.shape[1]) + 1), dtype=numpy.float32)
    wgt[:wgt.shape[0] - 1, :wgt.shape[1] - 1] = wg
    #
    # Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
    w_fl = numpy.zeros(shape=((wgt.shape[0]), (wgt.shape[1])), dtype=numpy.float32)
    w_fl = w_fl - 1  # to give -1 to NoData after the while loop calculation
    #
    # Let's go for the search and algo-rhytm for the weighted-Flow-Length
    nd = numpy.where(
        fdir8 < 0)  # fast coordinates all the NoData values, starting from them to go forward and compute flow length
    #
    y = nd[0]  # rows, NoData indexes
    x = nd[1]  # columns, NoData indexes pay attention not to invert values !!!!!!!!!!!!!!
    #
    # initializing lists for outlet and moving cell coordinates, in function of their position
    yc1 = []
    yc2 = []
    yc3 = []
    yc4 = []
    yc5 = []
    yc6 = []
    yc7 = []
    yc8 = []
    xc1 = []
    xc2 = []
    xc3 = []
    xc4 = []
    xc5 = []
    xc6 = []
    xc7 = []
    xc8 = []
    #
    #   Flow Directions TauDEM
    #   4   3   2
    #   5   -   1
    #   6   7   8
    #
    #   Draining in Direction Matrix
    #   8   7   6
    #   1   -   5
    #   2   3   4
    #
    #
    i1 = fdir8[y, x - 1]  # Searching for NoData with cells draining into them, 8 directions
    d1 = numpy.where(i1 == 1)  # l
    yc1.extend(y[d1])  # coordinates satisfying the conditions
    xc1.extend(x[d1])
    w_fl[yc1, xc1] = 0  # initialize flow length at cells draining to NoData
    #
    i2 = fdir8[y + 1, x - 1]  # Searching for NoData with cells draining into them, 8 directions
    d2 = numpy.where(i2 == 2)  # lrad2
    yc2.extend(y[d2])  # coordinates satisfying the conditions
    xc2.extend(x[d2])
    w_fl[yc2, xc2] = 0  # initialize flow length at cells draining to NoData
    #
    i3 = fdir8[y + 1, x]  # Searching for NoData with cells draining into them, 8 directions
    d3 = numpy.where(i3 == 3)  # l
    yc3.extend(y[d3])  # coordinates satisfying the conditions
    xc3.extend(x[d3])
    w_fl[yc3, xc3] = 0  # initialize flow length at cells draining to NoData
    #
    i4 = fdir8[y + 1, x + 1]  # Searching for NoData with cells draining into them, 8 directions
    d4 = numpy.where(i4 == 4)  # lrad2
    yc4.extend(y[d4])  # coordinates satisfying the conditions
    xc4.extend(x[d4])
    w_fl[yc4, xc4] = 0  # initialize flow length at cells draining to NoData
    #
    i5 = fdir8[y, x + 1]  # Searching for NoData with cells draining into them, 8 directions
    d5 = numpy.where(i5 == 5)  # l
    yc5.extend(y[d5])  # coordinates satisfying the conditions
    xc5.extend(x[d5])
    w_fl[yc5, xc5] = 0  # initialize flow length at cells draining to NoData
    #
    i6 = fdir8[y - 1, x + 1]  # Searching for NoData with cells draining into them, 8 directions
    d6 = numpy.where(i6 == 6)  # lrad2
    yc6.extend(y[d6])  # coordinates satisfying the conditions
    xc6.extend(x[d6])
    w_fl[yc6, xc6] = 0  # initialize flow length at cells draining to NoData
    #
    i7 = fdir8[y - 1, x]  # Searching for NoData with cells draining into them, 8 directions
    d7 = numpy.where(i7 == 7)  # l
    yc7.extend(y[d7])  # coordinates satisfying the conditions
    xc7.extend(x[d7])
    w_fl[yc7, xc7] = 0  # initialize flow length at cells draining to NoData
    #
    i8 = fdir8[y - 1, x - 1]  # Searching for NoData with cells draining into them, 8 directions
    d8 = numpy.where(i8 == 8)  # lrad2
    yc8.extend(y[d8])  # coordinates satisfying the conditions
    xc8.extend(x[d8])
    w_fl[yc8, xc8] = 0  # initialize flow length at cells draining to NoData
    #
    # start =clock()#
    count = 1  # "0" passage already done during the previous step
    while len(yc1) or len(yc2) or len(yc3) or len(yc4) or len(yc5) or len(yc6) or len(yc7) or len(yc8) > 0:
        # Converting into array to be able to do operations
        yyc1 = numpy.asarray(yc1)
        xxc1 = numpy.asarray(xc1)
        yyc2 = numpy.asarray(yc2)
        xxc2 = numpy.asarray(xc2)
        yyc3 = numpy.asarray(yc3)
        xxc3 = numpy.asarray(xc3)
        yyc4 = numpy.asarray(yc4)
        xxc4 = numpy.asarray(xc4)
        yyc5 = numpy.asarray(yc5)
        xxc5 = numpy.asarray(xc5)
        yyc6 = numpy.asarray(yc6)
        xxc6 = numpy.asarray(xc6)
        yyc7 = numpy.asarray(yc7)
        xxc7 = numpy.asarray(xc7)
        yyc8 = numpy.asarray(yc8)
        xxc8 = numpy.asarray(xc8)
        #
        # Now I can do operations and moving towards the right cell!!!!!!!! Weighting flow length, weights are half sum of pixels weight * travelled length
        # I'm choosing the directions accordingly to Flow_dir step by step going from outlet-nodata to the ridges,
        # each time account for distance (l or l*rad2) multiplied by the half of the weights of the 2 travelled cells.
        # Then, with variables substitution I'm moving a step further, and adding the previous pixel value to the new calculated.
        #
        yyc1 = yyc1
        xxc1 = (xxc1 - 1)  # l
        yyc2 = (yyc2 + 1)
        xxc2 = (xxc2 - 1)  # lrad2
        yyc3 = (yyc3 + 1)
        xxc3 = xxc3  # l
        yyc4 = (yyc4 + 1)
        xxc4 = (xxc4 + 1)  # lrad2
        yyc5 = yyc5
        xxc5 = (xxc5 + 1)  # l
        yyc6 = (yyc6 - 1)
        xxc6 = (xxc6 + 1)  # lrad2
        yyc7 = (yyc7 - 1)
        xxc7 = xxc7  # l
        yyc8 = (yyc8 - 1)
        xxc8 = (xxc8 - 1)  # lrad2
        #
        if count == 1:  # first run zero, like TauDEM
            if len(yyc1) > 0:
                w_fl[yyc1, xxc1] = 0
            else:
                pass
            if len(yyc2) > 0:
                w_fl[yyc2, xxc2] = 0
            else:
                pass
            if len(yyc3) > 0:
                w_fl[yyc3, xxc3] = 0
            else:
                pass
            if len(yyc4) > 0:
                w_fl[yyc4, xxc4] = 0
            else:
                pass
            if len(yyc5) > 0:
                w_fl[yyc5, xxc5] = 0
            else:
                pass
            if len(yyc6) > 0:
                w_fl[yyc6, xxc6] = 0
            else:
                pass
            if len(yyc7) > 0:
                w_fl[yyc7, xxc7] = 0
            else:
                pass
            if len(yyc8) > 0:
                w_fl[yyc8, xxc8] = 0
            else:
                pass
        else:
            w_fl[yyc1, xxc1] = w_fl[yc1, xc1] + (cell_s * ((wgt[yc1, xc1] + wgt[yyc1, xxc1]) / 2))
            w_fl[yyc2, xxc2] = w_fl[yc2, xc2] + (cell_s * math.sqrt(2) * ((wgt[yc2, xc2] + wgt[yyc2, xxc2]) / 2))
            w_fl[yyc3, xxc3] = w_fl[yc3, xc3] + (cell_s * ((wgt[yc3, xc3] + wgt[yyc3, xxc3]) / 2))
            w_fl[yyc4, xxc4] = w_fl[yc4, xc4] + (cell_s * math.sqrt(2) * ((wgt[yc4, xc4] + wgt[yyc4, xxc4]) / 2))
            w_fl[yyc5, xxc5] = w_fl[yc5, xc5] + (cell_s * ((wgt[yc5, xc5] + wgt[yyc5, xxc5]) / 2))
            w_fl[yyc6, xxc6] = w_fl[yc6, xc6] + (cell_s * math.sqrt(2) * ((wgt[yc6, xc6] + wgt[yyc6, xxc6]) / 2))
            w_fl[yyc7, xxc7] = w_fl[yc7, xc7] + (cell_s * ((wgt[yc7, xc7] + wgt[yyc7, xxc7]) / 2))
            w_fl[yyc8, xxc8] = w_fl[yc8, xc8] + (cell_s * math.sqrt(2) * ((wgt[yc8, xc8] + wgt[yyc8, xxc8]) / 2))
            #
        # Reconstructing all x and y of this step and moving on upwards (Downstream if you think in GIS, right?)
        yy = []
        xx = []
        yy.extend(yyc1)
        xx.extend(xxc1)
        yy.extend(yyc2)
        xx.extend(xxc2)
        yy.extend(yyc3)
        xx.extend(xxc3)
        yy.extend(yyc4)
        xx.extend(xxc4)
        yy.extend(yyc5)
        xx.extend(xxc5)
        yy.extend(yyc6)
        xx.extend(xxc6)
        yy.extend(yyc7)
        xx.extend(xxc7)
        yy.extend(yyc8)
        xx.extend(xxc8)
        #
        yy = numpy.asarray(yy)
        xx = numpy.asarray(xx)
        #
        i1 = fdir8[yy, xx - 1]  # Searching for cells draining into them, 8 directions
        d1 = numpy.where(i1 == 1)  # l
        yc1 = yy[d1]  # coordinates satisfying the conditions, HERE i NEED TO ADD ACTUAL LENGTH VALUE + PREVIOUS ONE
        xc1 = xx[d1]
        #
        i2 = fdir8[yy + 1, xx - 1]  # Searching for cells draining into them, 8 directions
        d2 = numpy.where(i2 == 2)  # lrad2
        yc2 = yy[d2]  # coordinates satisfying the conditions
        xc2 = xx[d2]
        #
        i3 = fdir8[yy + 1, xx]  # Searching for cells draining into them, 8 directions
        d3 = numpy.where(i3 == 3)  # l
        yc3 = yy[d3]  # coordinates satisfying the conditions
        xc3 = xx[d3]
        #
        i4 = fdir8[yy + 1, xx + 1]  # Searching for cells draining into them, 8 directions
        d4 = numpy.where(i4 == 4)  # lrad2
        yc4 = yy[d4]  # coordinates satisfying the conditions
        xc4 = xx[d4]
        #
        i5 = fdir8[yy, xx + 1]  # Searching for cells draining into them, 8 directions
        d5 = numpy.where(i5 == 5)  # l
        yc5 = yy[d5]  # coordinates satisfying the conditions
        xc5 = xx[d5]
        #
        i6 = fdir8[yy - 1, xx + 1]  # Searching for cells draining into them, 8 directions
        d6 = numpy.where(i6 == 6)  # lrad2
        yc6 = yy[d6]  # coordinates satisfying the conditions
        xc6 = xx[d6]
        #
        i7 = fdir8[yy - 1, xx]  # Searching for cells draining into them, 8 directions
        d7 = numpy.where(i7 == 7)  # l
        yc7 = yy[d7]  # coordinates satisfying the conditions
        xc7 = xx[d7]
        #
        i8 = fdir8[yy - 1, xx - 1]  # Searching for cells draining into them, 8 directions
        d8 = numpy.where(i8 == 8)  # lrad2
        yc8 = yy[d8]  # coordinates satisfying the conditions
        xc8 = xx[d8]
        count = count + 1
    #
    #
    elapsed = (time.time() - start)  # computational time in secs
    print(time.strftime("%d/%m/%y %H:%M:%S    "), "Process concluded successfully \n", "%.2f" % elapsed,
          'seconds for Weighted-Flow Length calculation with ', int(count), ' iterations')  # truncating the precision
    #
    #
    w_fl = w_fl[1:w_fl.shape[0] - 1, 1:w_fl.shape[
                                           1] - 1]  # reshaping weighted flow length, we need this step to homogenize matrices dimensions!!!!!!!!!!
    #
    w_fl_ds = gdal.GetDriverByName('GTiff').Create((dir_path + "/w_flow_length.tif"), const_ar.shape[1],
                                                   const_ar.shape[0], 1, gdal.GDT_Float32)
    w_fl_ds.SetGeoTransform(geoinf)
    w_fl_ds.SetProjection(proj)
    w_fl_ds.GetRasterBand(1).WriteArray(w_fl, 0, 0)
    w_fl_ds = None
    del w_fl_ds
    del fdir8
    #
    # Downslope component, imposing 1 values where weighted flow length==0, thus avoiding division by 0
    d_down_ar = w_fl
    del w_fl
    d_down_ar[d_down_ar == 0] = 1
    d_down_ar[numpy.where(
        target_ar == 10)] = numpy.NaN  # repeating NOData setting where I have the target, otherwise it will misplace values in the output raster... Not really sure actually
    #
    """
    --------------------------------------
    WORKING ON D_UP COMPONENT
    --------------------------------------
    """
    #
    # Calculating D_inf flow accumulation to obtain dtmfilloksca
    # os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca "+filename[0:-4]+"sca.tif -nc")#windows case
    os.system(("mpiexec -n 8 AreaDinf -ang ".lower()) + filename[0:-4] + "ang.tif" + " -sca " + filename[
                                                                                                0:-4] + "sca.tif -nc")  # unix case
    #
    # dtmfilloksca file modifications--> conversion to array for manipulations
    tif_dtmsca = osgeo.gdal.Open(
        filename[0:-4] + "sca.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input tif_dtmsca dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening tif_dtmsca was successful!")
    # Array
    tif_sca_ar = tif_dtmsca.ReadAsArray()
    #
    acc_final_ar = tif_sca_ar / const_ar
    #
    # Calculating accW-->Dinf contributing area wrigthed with weigth raster
    # os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accW.tif -wg "+filewgt+" -nc")#windows case
    os.system(("mpiexec -n 8 AreaDinf -ang ".lower()) + filename[
                                                        0:-4] + "ang.tif" + " -sca " + dir_path + "/accW.tif -wg " + filewgt + " -nc")  # unix case
    #
    acc_w = osgeo.gdal.Open(dir_path + "/accW.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input acc_w dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening acc_w was successful!")
    # Array
    acc_w_ar = acc_w.ReadAsArray()
    del acc_w
    #
    # Calculating accS-->Dinf contributing area weighted with d8 slope value
    # os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accS.tif -wg "+filename[0:-4]+'s.tif'+" -nc")#windows case
    os.system(("mpiexec -n 8 AreaDinf -ang ".lower()) + filename[
                                                        0:-4] + "ang.tif" + " -sca " + dir_path + "/accS.tif -wg " + filename[
                                                                                                                     0:-4] + 's.tif' + " -nc")  # unix case
    #
    acc_s = osgeo.gdal.Open(dir_path + "/accS.tif")  # opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "couldn't open input acc_s dataset")
        sys.exit(1)
    else:
        print(time.strftime("%d/%m/%y %H:%M:%S    "), "opening acc_s was successful!")
    # Array
    acc_s_ar = acc_s.ReadAsArray()
    del acc_s
    del tif_dtmsca
    #
    # Computing C_mean as (accW+weight)/acc_final
    c_mean_ar = (acc_w_ar + tif_wgt_ar) / acc_final_ar
    #
    # Computing S mean (accS+s)/acc_final
    s_mean_ar = (acc_s_ar + tif_sd8_ar) / acc_final_ar
    #
    # Computing D_up as "%cmean.tif%" * "%smean.tif%" * SquareRoot("%ACCfinal.tif%" * "%resolution.tif%" * "%resolution.tif%")
    # i need a trick to elevate power of 2 only the positive values of the constant resolution array, otherwise all the -1 will become 1
    in_pos_cos = numpy.where(const_ar >= 0)
    const_ar[in_pos_cos] = (const_ar[in_pos_cos]) ** 2
    cell_area = const_ar  # change of variables, to be sure
    d_up_ar = c_mean_ar * s_mean_ar * numpy.sqrt(
        acc_final_ar * cell_area)  # to transform from unit values to square units
    d_up_ar[numpy.where(
        target_ar == 10)] = numpy.NaN  # repeating NoData setting where I have the target, otherwise it will misplace values in the output raster... Not really sure actually
    #
    # Computing Connectivity index
    ic_ar = numpy.log10(d_up_ar / d_down_ar)
    # ic_ar[ic_ar==0]=numpy.nan#setting NoData values, maybe incorrect for occurrences within the matrix
    #
    # saving/deleting upslope and downslope components if required, filename is fixed
    #
    del d_up_ar
    del d_down_ar
    #
    #
    o_ic_tg = out_ic_tg.replace('\\', '/')
    o_ic_tg = str(o_ic_tg)
    #
    #
    # Create the tiff raster of ic
    dst_ic_ds = gdal.GetDriverByName('GTiff').Create(o_ic_tg, const_ar.shape[1], const_ar.shape[0], 1,
                                                     gdal.GDT_Float32)
    dst_ic_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ic_ds.SetGeoTransform(geoinf)  # set coordinates and input raster projection
    dst_ic_ds.SetProjection(proj)
    dst_ic_ds.GetRasterBand(1).WriteArray(ic_ar, 0,
                                          0)  # write raster, previously only memory allocation
    dst_ic_ds.GetRasterBand(1).GetStatistics(0, 1)  # calculate statistics for visualization
    dst_ic_ds = None  # Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_ic_ds
    #
    #
    # Some cleanings removing TauDEM-created files
    os.remove(dir_path + "/accS.tif")
    os.remove(filename[0:-4] + "ang.tif")
    os.remove(filename[0:-4] + "angt.tif")
    os.remove(dir_path + "/w_flow_length.tif")
    os.remove(dir_path + "/accW.tif")
    os.remove(filename[0:-4] + "sca.tif")
    os.remove(filename[0:-4] + "ad8.tif")
    os.remove(filename[0:-4] + "s.tif")
    os.remove(filename[0:-4] + "p.tif")
    os.remove(filename[0:-4] + "sd8.tif")
    os.remove(dir_path + "/targets.tif")
    #
    elapsed_tot = time.time() - start_0
    print(time.strftime("%d/%m/%y %H:%M:%S    "), "Calculation finished in", "%.2f" % elapsed_tot,
          "seconds! \n Compliments!")
    return
    #
    #


## Execution
CavalliConnectivitytg(dtm_r_txt, cell_s_txt, tgt_file_txt, w_txt, ic_out_txt)
#
## Finish line!!!!
