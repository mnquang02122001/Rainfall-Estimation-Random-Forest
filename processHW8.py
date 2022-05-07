# File name: processAERONET.py
# File type: Python script file
# Copy Right: GeoInformatic research group - University of Engineering and Technology, Vietnam National University Hanoi
# Contact: Email: thanhntn@vnu.edu.vn
# Latest update: 18/2/2022
# Description: Use this file to process and match AERONET data with VIIRS and MAIAC AOD data by the station's coordinates
#              Input: original AERINET data directory
#					  daily, monthly mean VIIRS data
#					  daily, monthly mean MAIAC data
#              Output: matched daily mean data file of AERONET, VIIRS, MAIAC AOD
#                      matched monthly mean data file of AERONET, VIIRS, MAIAC AOD
#              How to use: Open this file, then open a conda environment command prompt or regular command prompt and install the required libraries listed in the import section.
#                          After that, close the file, move to the folder containing this code file in the command prompt and run: python processAERONET.py
# Python library required: pandas, numpy, scikit-learn, scipy, gdal



import pandas as pd
import numpy as np
import time
import os

from sklearn.metrics import r2_score
from datetime import datetime as dt
from scipy.stats import pearsonr
from glob import iglob
from osgeo import gdal




def dateparse1(x): return dt.strptime(x, "%d:%m:%Y %H:%M:%S")
def dateparse2(x): return dt.strptime(x, "%Y-%m-%d %H:%M:%S")



def readAERONET(filepath):
    if filepath is None:
        return None

    def dateparse(x): return datetime.datetime.strptime(x, "%d:%m:%Y %H:%M:%S")
    aeronet = pd.read_csv(filepath, na_values=[-999],
                          parse_dates={'time': [0, 1]},
                          date_parser=dateparse1)

    # data includes 3 date and time data columns, AOT values at 500 nm wavelength and the angstrom exponent
    data = (aeronet.dropna(axis=1, how='all').dropna(axis=0, how='all'))
    data = data.loc[:, ['time', 'AOD_500nm', '440-675_Angstrom_Exponent', 'Site_Latitude(Degrees)', 'Site_Longitude(Degrees)']]
            # .loc[:, ['times', 'AOT_500', '440-675Angstrom']])

    array = data.values
    rows, cols = array.shape

    # interpolate AOD values measured at wavelength 550nm from AOD values measured at 500nm using the angstromExponent
    # convert time string to datetime
    # delete the 3rd column after calculation
    for i in range(rows):
        array[i][1] = array[i][1]/(np.e ** ((- array[i][2]) * np.log(0.5/0.55)))

    array = np.delete(array, 2, 1)
    # print(array)
    df = pd.DataFrame(array, columns = ['time', 'AOD550', 'lat', 'lon'])
    df['AOD550'] = pd.to_numeric(df['AOD550'])
    return df


# function to check if is leap year
def is_leap_year(year):
    return (((year % 4) == 0) or (((year % 100) != 0) and ((year % 400) == 0)))


# convert date index (1 to 365) in a year into a time string
def convertJDayToDate(year, jDay):
    if(is_leap_year(year)):
        day = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336]
    else:
        day = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]

    for i in reversed(range(12)):
        if (jDay - day[i] >= 0):
            return str(year), str(i + 1).zfill(2), str(jDay + 1 - day[i]).zfill(2)


# extract value from a point (pixel) on the image
def getValueAroundStation(img_data, geotransform, lat, lon, width, shape = 'CIRCLE'):
    result = []

    if img_data is None:
        return None    

    offset_y = (lon - geotransform[0])/geotransform[1]
    offset_x = (lat - geotransform[3])/geotransform[5]

    # print(offset_x, offset_y)

    if shape == 'SQUARE':
        evenWidth = ((width%2)==0)

        halfWidth = int(width/2)

        xmin = halfWidth
        ymin = halfWidth

        if evenWidth:
            # CONVENTION:
	        # If width is even, we push the center to the lower left corner 
	        # and have the "extra" in the right & above.   
            center_offset_x = int(np.ceil(offset_x))
            center_offset_y = int(np.floor(offset_y))

            xmin = xmin - 1
            ymin = ymin - 1
        elif not evenWidth:
            # CONVENTION:
	        # If width is odd, we push the center to the nearest grid point
            center_offset_x = int(np.floor(offset_x) + 0.5)
            center_offset_y = int(np.floor(offset_y) + 0.5)

        for i in range(-xmin, halfWidth + 1):
            for j in range(-ymin, halfWidth + 1):
                real_x = center_offset_x - i
                real_y = center_offset_y + j
                if img_data[real_x][real_y] != -9999:
                    result.append(img_data[real_x][real_y])

    elif shape == 'CIRCLE':
        if width == 2:
            return None

        # if the width is even, that means we are dealing with a point interpolation
        # because grid interpolation has to be odd.

        # for an ODD WIDTH the reference point is the same as the center point and is the nearest grid point
        
        # for an EVEN WIDTH, we move the "reference" point, to the lower left grid point,
        # this means offsets are stored relative to the lower left corner of the true center.
        # but we find distances based on the the true center location when determining if an
        # offset is within the circle.
        evenWidth = ((width%2)==0)
        radius = (width-1)/2.0

        # print(radius)
        # need to increase the area we look at if the width is even, because
        # some valid offset points will actually be farther from the reference point
        # than the radius, because the reference point is offset from the true
        # center of the circle.
        maxOffset = int(np.floor(radius))

        if evenWidth:
            maxOffset = maxOffset + 1
            center_offset_x = np.ceil(offset_x)
            center_offset_y = np.floor(offset_y)
        elif not evenWidth:
            center_offset_x = int(np.floor(offset_x) + 0.5)
            center_offset_y = int(np.floor(offset_y) + 0.5)

        minOffset = int(np.floor(radius) * -1)

        # print(minOffset, maxOffset)

        for y in range(minOffset, maxOffset + 1):
            for x in range(minOffset, maxOffset + 1):
                double_x = float(x)
                double_y = float(y)

                real_x = int(center_offset_x - double_x)
                real_y = int(center_offset_y + double_y)

                if evenWidth:
                    # if width is even, the reference point is actually shifted 1/2 a grid spacing down and to the left,
                    # from the true center of the circle.
                    # so when we calculate distance, we need to subtract .5 so that the distance reflects the distance from the center
                    # of the circle, instead of the distance from the reference.
                    # for example - a circle with width == 4.  The reference point is the lower left corner of the center square.
                    # the point directly below that is at (0,-1), but it's actually (-.5, -1.5) from the center of the circle.
                    # another example - same circle.  The point directly to the right of the reference point is (1,0), but it's
                    # actually (.5,-.5) from the center.
                
                    double_x = double_x - 0.5 
                    double_y = double_y - 0.5  

                    # print(double_x, double_y)
                    real_x = int(center_offset_x - (double_x + 0.5))
                    real_y = int(center_offset_y + (double_y + 0.5))

                distance = np.sqrt((double_x * double_x) + (double_y * double_y))

                

                if distance <= radius and img_data[real_x][real_y] != -9999:
                    result.append(img_data[real_x][real_y])

    # print(result)
    if result:
        return np.average(result)
    else:
        return None



AERONET_PATH = r'Data sample/AERONET/'
VIIRS_PATH = r'Data sample/VIIRS/'
MAIAC_PATH = r'Data sample/MAIAC/'

############################################################################################
# # merge AERONET data of 2018 and 2019 (not available if running on sample data)
# merge = []
# for path in iglob(AERONET_PATH + '/*.csv'):
# 	for y in ['2018', '2019']:
# 		if y in path and 'lev20' in path:
# 			print(path)
# 			file = pd.read_csv(path)
# 			file['Date'] = file['Day_of_Year'].apply(lambda x: '-'.join(convertJDayToDate(int(y), x)))
# 			file['time'] = file['Date'] + file['Time(hh:mm:ss)']
# 			# file = file.drop(columns=['Date(dd:mm:yyyy)', 'Time(hh:mm:ss)'])
# 			merge.append(file)

# merge = pd.concat(merge)
# merge.to_csv(AERONET_PATH + '/merged 2018 - 2019.csv', index=False)
############################################################################################



############################################################################################
# call on the interpolation function and save data to a new file
data = readAERONET(AERONET_PATH + 'AERONET sample.csv')
data.to_csv(AERONET_PATH + 'AERONET sample interpolated.csv', index=False)
############################################################################################



############################################################################################
# read the interpolated data file
AOD_550 = pd.read_csv(AERONET_PATH + 'AERONET sample interpolated.csv',
	parse_dates=['time'],
	date_parser=dateparse2)

# extract datetime data to new columns
AOD_550['month'] = AOD_550['time'].apply(lambda x: int(x.month))
AOD_550['year'] = AOD_550['time'].apply(lambda x: int(x.year))
AOD_550['day'] = AOD_550['time'].apply(lambda x: int(x.day))

AOD_550['hour'] = AOD_550['time'].apply(lambda x: int(x.hour))

#calculate averaged AOD values of the AERONET station by hour, day, month
mean_by_hour = AOD_550.groupby(['hour', 'day', 'month', 'year'], as_index=False).mean()

mean_by_day = mean_by_hour.groupby(['day', 'month', 'year'], as_index=False).mean()
mean_by_day = mean_by_day.drop(columns=['hour'])
############################################################################################











############################################################################################
# merge AERONET data and satellite data by day based on station coordinates
start_time = time.time()

for idx, row in mean_by_day.iterrows():
    # VIIRS AOD
    name = ''
    if row['month'] < 10:
        if row['day'] < 10:
            name = ''.join([str(int(row['year'])), str(0), str(int(row['month'])), str(0), str(int(row['day']))])
        else:
            name = ''.join([str(int(row['year'])), str(0), str(int(row['month'])), str(int(row['day']))])
    else:
        if row['day'] < 10:
            name = ''.join([str(int(row['year'])), str(int(row['month'])), str(0), str(int(row['day']))])
        else:
            name = ''.join([str(int(row['year'])), str(int(row['month'])), str(int(row['day']))])

    try:
        img = gdal.Open(VIIRS_PATH + 'daily average/' + str(int(row['year'])) + '/' + name + '.tif')
        img_data = img.GetRasterBand(1).ReadAsArray()
        img_nodata = img.GetRasterBand(1).GetNoDataValue()

        geotransform = img.GetGeoTransform()

        result = getValueAroundStation(img_data, geotransform, row['lat'], row['lon'], 1, 'SQUARE')
        mean_by_day.at[idx, 'AOD_VIIRS'] = result
    except:
        pass


    # MAIAC AOD
    name = ''
    if row['month'] < 10:
        name = str(int(row['year'])) + str(0) + str(int(row['month'])) + '_' + str(int(row['day']))
    else:
        name = str(int(row['year'])) + str(int(row['month'])) + '_' + str(int(row['day']))

    try:
        img = gdal.Open(MAIAC_PATH + 'daily average/' + str(int(row['year'])) + '/MAIAC_composite_' + name + '.tif')
        img_data = img.GetRasterBand(1).ReadAsArray()
        img_nodata = img.GetRasterBand(1).GetNoDataValue()

        geotransform = img.GetGeoTransform()

        result = getValueAroundStation(img_data, geotransform, row['lat'], row['lon'], 1, 'SQUARE') * 0.001
        mean_by_day.at[idx, 'AOD_MAIAC'] = result
    except:
        pass

end_time = time.time()
print(end_time - start_time)
print(mean_by_day)
mean_by_day.to_csv(AERONET_PATH + 'AOD_550 AERONET vs VIIRS vs MAIAC daily.csv', index=False)
############################################################################################





