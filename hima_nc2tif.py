# cài đặt thư viện gdal
# https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
# https://anaconda.org/conda-forge/gdal

from osgeo import gdal
from osgeo import osr
import numpy as np

# Đường dẫn đễn file ví dụ
test_fpath = r'data\NC_H08_20190831_2350_R21_FLDK.02401_02401.nc'

# Xem info của file
info = gdal.Info(test_fpath)
# print(info)

# Đọc file dữ liệu sử dụng gdal
ds = gdal.Open(test_fpath, gdal.GA_ReadOnly)

# Lấy ra danh sách các band trong file dữ liệu
list_datasets = ds.GetSubDatasets()
print(list_datasets)
# lựa chọn band1 để trích xuất, thực tế bước này cần vòng lặp để chọn từ b1 -> 16
band1 = list_datasets[0][0]

# Xem info của band 1
info = gdal.Info(band1)
# print(info)

# Đọc band 1
sub_ds = gdal.Open(band1)
data = sub_ds.ReadAsArray()

# Lấy thông tin tọa độ của band 1
geotransform = sub_ds.GetGeoTransform()

# tên file tif output
outfname = test_fpath.replace('.nc', '_band1.tif')

# tạo file tif output
driver = gdal.GetDriverByName('GTiff')
# Định nghĩa file output (đường dẫn đến file), kích thước file (shape[0], shape[1]), số band trong file (1), kiểu dữ liệu (float32)
dst_ds = gdal.GetDriverByName('GTiff').Create(
    outfname, data.shape[0], data.shape[1], 1, gdal.GDT_Float32)
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromEPSG(4326)
# ĐỊnh nghĩa hệ quy chiếu WGS84 cho file
dst_ds.SetProjection(outRasterSRS.ExportToWkt())
dst_ds.SetGeoTransform(geotransform)  # Đặt thông tin tọa độ cho file

dst_ds. GetRasterBand(1).SetNoDataValue(-99999)  # Đặt giá trị nodata
dst_ds. GetRasterBand(1).WriteArray(data)  # Ghi dữ liệu ra file trên máy
dst_ds.FlushCache()
dst_ds = None

###########################################################################################
# Cắt ảnh tif trên vùng nghiên cứu
data_path = outfname  # input là đầu ra của bước trên

# Đặt tên file output sau khi cắt
clipped_fpath = data_path.replace('.tif', '_clipped.tif')

# doc: https://gdal.org/python/ --> Modules: osgeo.gdal --> functions: Warp
# outputBounds --- output bounds as (minX, minY, maxX, maxY)
# srcDS --- a Dataset object or a filename
# destName --- Output dataset name
# resampleAlg --- resampling mode

# Cần chỉnh lại tham số outputBounds để chọn vùng nghiên cứu
# (minX, minY, maxX, maxY) là tọa độ phạm vi vùng nghiên cứu
# minX, maxX là khoảng kinh độ (longitude)
# minY, maxY là khoảng vĩ độ (latitude)
# xRes, yRes định nghĩa độ phân giải của ảnh output, đặt là (0.05, -0.05), đơn vị là degree, 0.05 degree ~ 5km
# resampleAlg: phương pháp tái tạo mẫu, chọn phương pháp láng giềng gần nhất
gdal.Warp(destNameOrDestDS=clipped_fpath, srcDSOrSrcDSTab=data_path,  xRes=0.05, yRes=-
          0.05, outputBounds=[100.1, 6.1, 111.8, 25.6], resampleAlg=gdal.GRA_NearestNeighbour)
