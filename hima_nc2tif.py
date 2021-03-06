from osgeo import gdal
from osgeo import osr
import numpy as np


def convert(time):
    test_fpath = "E:\\RainData\\28\\NC_H08_20190828_" + \
        time + "_R21_FLDK.02401_02401.nc"
    info = gdal.Info(test_fpath)

    ds = gdal.Open(test_fpath, gdal.GA_ReadOnly)

    list_datasets = ds.GetSubDatasets()

    outfname = ''

    for i in range(17):
        if i == 3:
            continue
        band = list_datasets[i][0]
        info = gdal.Info(band)

        sub_ds = gdal.Open(band)
        data = sub_ds.ReadAsArray()

        geotransform = sub_ds.GetGeoTransform()
        if i < 3:
            outfname = test_fpath.replace('.nc', '_band' + str(i + 1) + '.tif')
        elif i > 3:
            outfname = test_fpath.replace('.nc', '_band' + str(i) + '.tif')

        driver = gdal.GetDriverByName('GTiff')

        dst_ds = gdal.GetDriverByName('GTiff').Create(
            outfname, data.shape[0], data.shape[1], 1, gdal.GDT_Float32)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(4326)

        dst_ds.SetProjection(outRasterSRS.ExportToWkt())
        dst_ds.SetGeoTransform(geotransform)

        dst_ds.GetRasterBand(1).SetNoDataValue(-99999)
        dst_ds.GetRasterBand(1).WriteArray(data)
        dst_ds.FlushCache()
        dst_ds = None

        clipped_fpath = outfname.replace('.tif', '_clipped.tif')
        gdal.Warp(destNameOrDestDS=clipped_fpath, srcDSOrSrcDSTab=outfname,  xRes=0.05, yRes=-
                  0.05, outputBounds=[104.35, 19.23, 106.1, 20.68], resampleAlg=gdal.GRA_NearestNeighbour)


for hour in range(24):
    str_hour = '0' + str(hour) if (hour >= 0 and hour <= 9) else str(hour)
    for min in range(0, 60, 10):
        str_min = '0' + str(min) if (min == 0) else str(min)
        time = str_hour + str_min
        try:
            convert(time)
        except Exception as e:
            print(e)
            print(time)
            continue
