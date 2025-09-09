#SCRIPT TO USE WHEN DATA ARE UPLOADED MANUALLY FROM THE USER

#Import libraries
import pandas as pd
import numpy as np
import rasterio as ra
import zipfile
from zipfile import ZipFile
import re
from pathlib import Path
import os
from PIL import Image
from osgeo import gdal_array
from matplotlib import pyplot as plt
from numpy import asarray 
import scipy
from scipy import linalg
from scipy.linalg import svd
from osgeo import gdal, osr, ogr
from xml.dom import minidom
from xml.dom.minidom import parse, parseString
from time import gmtime, strftime


method = "Strang (1988); Lundgren et al. (2004); Manzo et al. (2006); De Luca et al. (2017)"
spatial_resolution = "30 m"
folder = ""

print("***EAST_WEST AND VERTICAL EXTRACTION***\n by Occhipinti M., De Luca C., Manunta M., Monterroso M., Casu F.\n Released by: IREA-CNR")
print("\n.\n.\n.")
print("Start of the processing\n.\n.")
print("Check of input data...")
# Check process of existance of files and of their effective usage
def checklist(file):
    element = file[0]
    check = True
    for item in file:
        if element != item:
            check = False
            break
    if (check == True):
        print(NameError)
    else:
        print("Input data existing")
        return file


data_input = []
orbit_data_inu = []
cos_input_data = []
orbit_data_cos = []
map_los_vec = []
relative_orbit_number = []
print("\n.\n.\n.")

for x in os.listdir():
    if "InU" in x:
        inu_data = os.listdir(x)
        for y in inu_data:
            if ".tif" in y:
                y = x + "/" + y
                data_input.append(y)
                product_format = "GEOTIFF"

        for y in inu_data:
            if ".metadata" in y:
                metadata_path = (str(x) + "/" + y)    
                metadata = open(metadata_path, "r")
                for line in metadata:

                    if "Map_of_LOS_vector: " in line:
                        map_los_vec_line = line.strip()
                        map_los_vector = map_los_vec_line[19:] #COMMON
                        map_los_vec.append(map_los_vector)

                    if "Relative_orbit_number: " in line:
                        relative_orbit_number_line = line.strip()
                        relative_orbit_num = relative_orbit_number_line[24:]
                        relative_orbit_number.append(relative_orbit_num)


        for z in inu_data:
            if ".metadata" in z:
                metadata_path = (str(x) + "/" + z)    
                metadata = open(metadata_path, "r")  
                if "ASCENDING" in metadata.read():
                    ascending = "ascending product"
                    orbit_data_inu.append(ascending)
                else:
                    descending = "descending product"
                    orbit_data_inu.append(descending)



for y in os.listdir():
    if "Cos" in y:
        cos_data = os.listdir(y)
        for i in cos_data:
            if ".tif" in i:
                i = y + "/" + i
                cos_input_data.append(i)
                        
        for w in cos_data:
            if ".metadata" in w:
                metadata_path = (str(y) + "/" + w)    
                metadata = open(metadata_path, "r")  
                if "ASCENDING" in metadata.read():
                    ascending = "ascending product"
                    orbit_data_cos.append(ascending)
                else:
                    descending = "descending product"
                    orbit_data_cos.append(descending)

n = len(data_input)
n_cos = len(cos_input_data)
checklist(orbit_data_inu)
checklist(orbit_data_cos)
print("\n.\n.\n.")
print("Organization of data...")
if n == n_cos:
    True
    n = n_cos
else:
    print(NameError)

min_lon = np.zeros(n)
min_lat = np.zeros(n)
max_lon = np.zeros(n)
max_lat = np.zeros(n)

min_lon_cos = np.zeros(n)
min_lat_cos = np.zeros(n)
max_lon_cos = np.zeros(n)
max_lat_cos = np.zeros(n)
print("\n.\n.\n.")
print("Pre-processing of InU data: cutting of the Region Of Interest")
# Processing for Deformation Data: cutting of the common area of deformation between the n files of input
# 1) Opening of the data and extraction of shape informations
for i in range(n):
    tif_file = ra.open(data_input[i])
    x_dimension = tif_file.width
    y_dimension = tif_file.height
    box = tif_file.bounds
    min_longitude = box[0]
    max_latitude = box[3]
    spacing_lat, spacing_lon = tif_file.res

    min_lon[i] = min_longitude
    max_lat[i] = max_latitude
    max_lon[i] = min_longitude + ((x_dimension - 1)*spacing_lon)
    min_lat[i] = max_latitude - ((y_dimension - 1)*spacing_lat)

# 2) Definition of the common area
aoi_min_lon = max(min_lon)
aoi_min_lat = max(min_lat)
aoi_max_lon = min(max_lon)
aoi_max_lat = min(max_lat)


dim_range_lat = int(np.round(((aoi_max_lat-aoi_min_lat)/spacing_lat)+1))
dim_range_lon = int(np.round(((aoi_max_lon-aoi_min_lon)/spacing_lon)+1))


matrix_cut_deformation = np.empty((n, dim_range_lat, dim_range_lon))
matrix_cut_cos = np.copy(matrix_cut_deformation)
north_matrix_cos = np.copy(matrix_cut_deformation)
east_matrix_cos = np.copy(matrix_cut_deformation)
up_matrix_cos = np.copy(matrix_cut_deformation)
mask_comm = np.copy(matrix_cut_deformation)

print("InU Region Of Interest calculated correctly.")
print("\n.\n.\n.")
print("Overlapping the input data...")

# 3) Cut of the layers of input data in a common area
for i in range(n):
    tif_file = ra.open(data_input[i])
    x_dimension = tif_file.width
    y_dimension = tif_file.height
    box = tif_file.bounds
    min_longitude = box[0]
    max_latitude = box[3]
    spacing_lat, spacing_lon = tif_file.res
    min_lon[i] = min_longitude
    max_lat[i] = max_latitude
    max_lon[i] = min_longitude + ((x_dimension - 1)*spacing_lon)
    min_lat[i] = max_latitude - ((y_dimension - 1)*spacing_lat)

    x_start = int(np.round((min_lon[i] - aoi_min_lon)/spacing_lon))
    x_end = int(np.round((max_lon[i] - aoi_max_lon)/spacing_lon))
    y_start = int(np.round((min_lat[i] - aoi_min_lat)/spacing_lat))
    y_end = int(np.round((max_lat[i] - aoi_max_lat)/spacing_lat))

    tif_file = gdal.Open(data_input[i])
    array_data = np.array(tif_file.ReadAsArray())
    new_matrix_deformation = array_data[0 - y_start:y_dimension - y_end, 0 - x_start:x_dimension - x_end]

    matrix_cut_deformation[i, :, :] = new_matrix_deformation # Common area of deformation

print("Overlapping completed successfully.")
print("\n.\n.\n.")
print("Pre-processing of CosNEU data: cutting of the Region Of Interest")
# Processing for Cos data: repetition of the previous steps
for j in range(n):
    cos_file = ra.open(cos_input_data[j])
    x_dimension_cos = cos_file.width 
    y_dimension_cos = cos_file.height 
    box_cos = cos_file.bounds
    min_longitude_cos = box_cos[0]
    max_latitude_cos = box_cos[3]
    spacing_lat_cos, spacing_lon_cos = cos_file.res
    min_lon_cos[j] = min_longitude_cos
    max_lat_cos[j] = max_latitude_cos
    min_lat_cos[j] = max_latitude_cos - ((y_dimension_cos - 1)*spacing_lat_cos)
    max_lon_cos[j] = min_longitude_cos + ((x_dimension_cos - 1)*spacing_lon_cos)

    x_start_cos = int(np.round((min_lon_cos[j] - aoi_min_lon)/spacing_lon_cos))
    x_end_cos = int(np.round((max_lon_cos[j] - aoi_max_lon)/spacing_lon_cos))
    y_start_cos = int(np.round((min_lat_cos[j] - aoi_min_lat)/spacing_lat_cos))
    y_end_cos = int(np.round((max_lat_cos[j] - aoi_max_lat)/spacing_lat_cos))

    cos_file = gdal.Open(cos_input_data[j])
    array_cos = np.array(cos_file.ReadAsArray())
    print("Overlapping CosNEU data...")
    new_matrix_cos = array_cos[:, 0-y_start_cos:y_dimension_cos - y_end_cos, 0 - x_start_cos:x_dimension_cos - x_end_cos]

    
    north_matrix_cos[j, :, :] = new_matrix_cos[0, :, :]
    east_matrix_cos[j, :, :] = new_matrix_cos[1, :, :]
    up_matrix_cos[j, :, :] = new_matrix_cos[2, :, :]     # 3 Arrays (North, Up and East), but only 2 will be used (Up, East) because of the near-polar orbit of the satellite
print("Overlapping of CosNEU data completed successfully.")
print("\n.\n.\n.")
print("Hooking InU and CosNEU...")
# Construction of a common mask that contains the points in which there is effective deformation for all the deformation files  
common_mask = np.empty((n, dim_range_lat, dim_range_lon))
for i in range(n):
    mask_zeros = np.copy(matrix_cut_deformation[i])
    mask_zeros[mask_zeros != 0 ] = 1
    common_mask[i, :, :] = mask_zeros

mask_prod_zeros = np.prod(common_mask, axis = 0) # product of the arrays along the z dimensions inside common_mask -- the mask is a 2D array of 0 and 1

# Basing on the indices where the values are 1, find inside "matrix_cut_deformation" those values (ascending(x;y) and descending(x:y)) that correspond to that index with value = 1, calculate their module and their average
ind_good = np.argwhere(mask_prod_zeros == 1) # pixels where there is deformation

for i in range(n):
    points = np.abs(matrix_cut_deformation[:, [index[0] for index in ind_good], [index[1] for index in ind_good]])
    mean_values = np.mean(points, axis=0)

tmp_array = np.ma.empty(mask_prod_zeros.shape)
tmp_array.mask = True  

for i, index in enumerate(ind_good):
    tmp_array[index[0], index[1]] = mean_values[i]
    tmp_array.mask[index[0], index[1]] = False

min_value = np.amin(tmp_array)
index_points = np.argwhere(tmp_array == min_value)[0]

for i in range(n):
    hookup = matrix_cut_deformation[:, index_points[0], index_points[1]] # values that give as mean value "min_value"

hookup = np.array(hookup)

# Subtract the hookup vector to the matrix_cut_deformation to create a point of hookup
hooked_matrix_cut_deformation = np.subtract(matrix_cut_deformation, hookup[:, np.newaxis, np.newaxis])
print("InU and CosNEU have been hooked!")
print("\n.\n.\n.")


# # CASE IN WHICH THE USER INSERTS MANUALLY A HOOKUP POINT IN LAT LON
# latitude_hookup = 40.0001
# longitude_hookup = 22.1810

# posx_hookup = int(np.round(abs(longitude_hookup - aoi_min_lon)/spacing_lon))
# posy_hookup = int(np.round(abs(latitude_hookup - aoi_min_lat)/spacing_lat))


# if mask_prod_zeros[posy_hookup, posx_hookup].any() == 1:
#     hookup = matrix_cut_deformation[:, posy_hookup, posx_hookup]
#     hooked_matrix_cut_deformation = np.subtract(matrix_cut_deformation, hookup[:, np.newaxis, np.newaxis])

# else:
#     box =  10
#     box_mask = mask_prod_zeros[posy_hookup - box:posy_hookup + box, posx_hookup - box:posx_hookup + box]

#     indices = np.where(box_mask == 1)
#     if len(indices[0]) > 0:
#         index_y = indices[0][0] + (posy_hookup - box)
#         index_x = indices[1][0] + (posx_hookup - box)
#         tmp_hookup = box_mask[indices[0][0], indices[1][0]]
#         index_tmp_hookup = (index_y, index_x)

#     hookup = matrix_cut_deformation[:, index_tmp_hookup[0], index_tmp_hookup[1]]
#     hooked_matrix_cut_deformation = np.subtract(matrix_cut_deformation, hookup[:, np.newaxis, np.newaxis]) 

# Extraction of the deformation for East and Up components


matrixUp = np.zeros((dim_range_lat, dim_range_lon))
matrixEast = np.zeros((dim_range_lat, dim_range_lon))
vector_points = np.empty((1, 2))

print("Retrieval of East-West and Up information...")
for pixel in ind_good:

    matrixB = np.empty((n, 2))
    matrixB[:, 0] = east_matrix_cos[:, pixel[0], pixel[1]]
    matrixB[:, 1] = up_matrix_cos[:, pixel[0], pixel[1]]


    u, s, vh = svd(matrixB, full_matrices = False)

    diag = np.zeros((2, 2))

    if s[0] < 0.0001:
        diag[0, 0] = 0
    else:
        diag[0, 0] = s[0]**-1     

    if s[1] < 0.0001:
        diag[1, 1] = 0
    else:
        diag[1, 1] = s[1]**-1   
    

    pseudo_inverse = np.dot(vh.transpose(), np.dot(diag, u.transpose())) 
    
    vector_points = hooked_matrix_cut_deformation[:, pixel[0], pixel[1]]
    components = np.dot(pseudo_inverse, vector_points)

    matrixEast[pixel[0], pixel[1]] = components[0]
    matrixUp[pixel[0], pixel[1]] = components[1]

print("East-West and Up information ready.")

bounding_box = "%s %s %s %s" % (aoi_max_lat, aoi_max_lon, aoi_min_lat, aoi_min_lon)
bounding_box_gml = "%s %s %s %s %s %s %s %s %s %s " % (aoi_min_lon, aoi_min_lat, aoi_max_lon, aoi_min_lat, aoi_max_lon, aoi_max_lat, aoi_min_lon, aoi_max_lat, aoi_min_lon, aoi_min_lat)
reference_point = str(hookup)
product_size = dim_range_lat*dim_range_lon
bounding_box_wkt = "POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s))" 

print("\n.\n.\n.")
print("Writing of GeoTIFF for East-West data...")
# Last step: save the file
def getGeoTransform(extent, dimlat, dimlon):
    x = (extent[2] - extent[0]) / dimlon
    y = (extent[3] - extent[1]) / dimlat
    return [extent[0], x, 0, extent[3], 0, -y]

# 1) Extension of the file
extent = [aoi_min_lon, aoi_min_lat, aoi_max_lon, aoi_max_lat]

# 2) Get GDAL driver GeoTiff
driver = gdal.GetDriverByName("GTiff")
data_type = gdal.GDT_Float32

# 3) Create a temporary grid
grid_data = driver.Create("grid_data", dim_range_lon, dim_range_lat, 1, data_type) 

# 4) Write data bands
grid_data.GetRasterBand(1).WriteArray(matrixEast)

# 5) Reference System
rs = osr.SpatialReference()
rs.ImportFromEPSG(4326)
grid_data.SetProjection(rs.ExportToWkt())
grid_data.SetGeoTransform(getGeoTransform(extent, dim_range_lat, dim_range_lon))

# 6) Set the name of the file and save the file
file_name_matrixEast = "EW_.tif"
driver.CreateCopy(file_name_matrixEast, grid_data, 0)

# 7) Close the file
driver = None
grid_data = None
 
# 8) Delete the temporary grid
import os                
os.remove("grid_data")


# Metadata writer
print("Metadata creation -- WORK IN PROGRESS PART, IF ANY ERROR OCCURS, COMMENT FROM LINE 384 TO LINE 1121")
ddss_id = "EW"
product_id = "EW_XXX"
preview_url = "xxx"
legend_url = "xxx"
product_url = "xxx"
kmz_url = "xxx"
Lookup_table_from_radar_coordinates_to_ground_coordinates = "xxx"
code_space = "xxx"
code_value = "xxx"



for x in os.listdir():
    if "InU" in x:
        inu_data = os.listdir()
        for y in inu_data:
            if ".metadata" in y:
                metadata_path = (str(x) + "/" + y)    
                metadata = open(metadata_path, "r")
                for line in metadata:
                    if "Sensor: " in line:
                        sensor_line = line.strip()
                        sensor = sensor_line[8:]
                    if "License: " in line:
                        license_line = line.strip()
                        license = license_line[9:]
                    if "User_ID: " in line:
                        user_id_line = line.strip()
                        user_id = user_id_line[9:]
                    if "Software_version: " in line:
                        software_version_line = line.strip()
                        software_version = software_version_line[18:]
                    if "Geographic_CS_type_code: " in line:
                        reference_system_line = line.strip()
                        reference_system = reference_system_line[25:]
                    if "Applied_algorithm_description: " in line:
                        applied_algorithm_description_line = line.strip()
                        applied_algorithm_description = applied_algorithm_description_line[31:]
                    if "Main_reference: " in line:
                        main_reference_line = line.strip()
                        main_reference = main_reference_line[17:]
                    if "Service_used_for_generation: " in line:
                        service_used_for_generation_line = line.strip()
                        service_generation = service_used_for_generation_line[39:]
                    if "Used_DEM: " in line:
                        used_dem_line = line.strip()
                        used_dem = used_dem_line[10:]
                    if "Applied_unwrapping_algorithm: " in line:
                        applied_unwrapping_algorithm_line = line.strip()
                        applied_unwrapping_algorithm = applied_unwrapping_algorithm_line[30:]
                    if "Mode: " in line:
                        mode_line = line.strip()
                        mode = mode_line[6:]
                    if "Antenna_side: " in line:
                        antenna_side_line = line.strip()
                        antenna_side = antenna_side_line[14:]
                    if "Wavelength: " in line:
                        wavelength_line = line.strip()
                        wavelength = wavelength_line[11:]
                    if "Value_unit: " in line:
                        value_unit_line = line.strip()
                        value_unit = value_unit_line[12:]
                    if "AffiliationIdentifier: " in line:
                        affiliationidentifier_line = line.strip()
                        affiliationidentifier = affiliationidentifier_line[23:]

                    
ddss_id_string = "DDSS_ID: " + ddss_id
product_id_string = "Product_ID: " + product_id
product_format_string = "Product_format: " + product_format
product_size_string = "Product_size: " + str(product_size)
preview_url_string = "Preview_url: " + preview_url
legend_url_string = "Legend_url: " + legend_url 
product_url_string = "Product_url: " + product_url 
bounding_box_string = "Bounding_box: " + bounding_box
bounding_box_wkt_string = "Bounding_box_wkt: POLYGON((" + bounding_box_wkt
license_string = "License: " + license
user_id_string = "User_ID: " + user_id
software_version_string = "Software_version: " + software_version
applied_algorithm_description_string = "Applied_algorithm_description: " + applied_algorithm_description 
main_reference_string = "Main_reference: " + main_reference 
# date_of_measurement_start_string = "Date_of_measurement_start: " + date_of_measurement_start
# date_of_measurement_end_string = "Date_of_measurement_end: " + date_of_measurement_end
date_of_production_string = "Date_of_production: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
service_used_for_generation_string = "Service_used_for_generation: " + service_generation
reference_system_string = "Geographic_CS_type_code: " + reference_system
used_dem_string = "Used_DEM: " + used_dem
perpendicular_baseline_string = "Perpendicular_baseline: N/a"
parallel_baseline_string = "Parallel_baseline: N/a"
along_track_baseline_string = "Along_track_baseline: N/a"
map_los_vec_string = "Map_of_LOS_vector: " + str(map_los_vec)
applied_unwrapping_algorithm_string = "Applied_unwrapping_algorithm: " + applied_unwrapping_algorithm
reference_point_string = "Reference_point: " + reference_point
spatial_resolution_string = "Spatial_resolution: "  + spatial_resolution #30 m
sensor_string = "Sensor: " + sensor
mode_string = "Mode: " + mode
antenna_side_string = "Antenna_side: " + antenna_side
relative_orbit_number_string = "Relative_orbit_number: " + str(relative_orbit_number) # serie di orbite combinate se ci sono nei metadati delle linee di vista (parametrizzato su n input data)
wavelength_string = "Wavelength: " + wavelength
value_unit_string = "Value_unit: " + value_unit 
number_of_looks_azimuth_string = "Number_of_looks_azimuth: N/a"
number_of_looks_range_string = "Number_of_looks_range: N/a" 
applied_filter_string = "Applied_filter: No_Filter"
kmz_url_string = "kmz_url: " + kmz_url
# affiliationidentifier_string = "AffiliationIdentifier: " + affiliationidentifier
method_string = "MethodDescription: " + method

metadata = [ddss_id_string, product_id_string, product_format_string, product_size_string, preview_url_string, legend_url_string, product_url_string, bounding_box_string, bounding_box_wkt_string, license_string, user_id_string, software_version_string, applied_algorithm_description_string,  main_reference_string, date_of_production_string, service_used_for_generation_string, reference_system_string, used_dem_string, perpendicular_baseline_string, parallel_baseline_string, along_track_baseline_string, map_los_vec_string, applied_unwrapping_algorithm_string, reference_point_string, spatial_resolution_string, sensor_string, mode_string, antenna_side_string, relative_orbit_number_string, wavelength_string, value_unit_string, number_of_looks_azimuth_string, number_of_looks_range_string, applied_filter_string, kmz_url_string, method_string]

metadata_file_path = "EW.metadata"

with open(metadata_file_path, "w") as metadata_file:
    for x in metadata:
        metadata_file.write(x + "\n")



root = minidom.Document()
xml = root.createElement("feed")
xml.setAttribute("xmlns", "http://www.w3.org/2005/Atom")
root.appendChild(xml)

title_child = root.createElement("title")
title_child.setAttribute("type= 'text'", "Template for ingestion of an EPOS product")
xml.appendChild(title_child)

id = root.createElement("id")
id.setAttribute("", "http://catalog.terradue.com/gep-epos/")
xml.appendChild(id)

entry = root.createElement("entry")
xml.appendChild(entry)

identifier = root.createElement("identifier")
identifier.setAttribute("xmlsn = http://purl.org/dc/elements/1.1/", product_id)
entry.appendChild(identifier)

title_2 = root.createElement("title")
title_2.setAttribute("type= 'text'", product_id)
entry.appendChild(title_2)

summary = root.createElement("summary type= 'html'")
entry.appendChild(summary)

data = root.createElement("![CDATA[")
summary.appendChild(data)

table = root.createElement("table")
summary.appendChild(table)

body = root.createElement("tbody")
table.appendChild(body)

tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
tr.appendChild(td)

a = root.createElement("a")
a.setAttribute("href=" + preview_url + " target= '_blank' title = 'View preview image'", None)
td.appendChild(a)

img = root.createElement("img")
img.setAttribute("align= 'left' border= '0' width= '100px' title= 'View legend' src=" + preview_url, None)
a.appendChild(img)

td = root.createElement("td")
tr.appendChild(td)

a = root.createElement("a")
a.setAttribute("href=" + legend_url + " target= '_blank' title = 'View preview image'", None)
td.appendChild(a)

img = root.createElement("img")
img.setAttribute("align= 'left' border= '0' width= '100px' title= 'View legend' src=" + legend_url, None)
a.appendChild(img)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Product Identifier")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", product_id)
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Platform")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", sensor_string[7:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Sensor Type")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", "RADAR")
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Product Type")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", ddss_id_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Processing Center")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", user_id_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Date of measurement start")
tr.appendChild(td)

# strong = root.createElement("strong")
# strong.setAttribute("", date_of_measurement_start)
# td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Date of measurement end")
tr.appendChild(td)

# strong = root.createElement("strong")
# strong.setAttribute("", date_of_measurement_end)
# td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Processor Name")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", service_used_for_generation_string[28:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Unit of measure")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", value_unit_string[11:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Referemce System")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", reference_system_string[24:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "License")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", license_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Reference point [lon lat]")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", reference_point_string[16:])
td.appendChild(strong)

tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "DEM")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", used_dem_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Method")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", applied_algorithm_description[20:])
td.appendChild(strong)


link1 = root.createElement("link")
link1.setAttribute("href =" + product_url + "rel= 'enclosure' title= 'Product URL' type='text/html'", None)
entry.appendChild(link1)
link2 = root.createElement("link")
link2.setAttribute("href =" + preview_url + "rel ='icon' title = 'Preview URL' type ='image/x-icon'", None)
entry.appendChild(link2)
link3 = root.createElement("link")
link3.setAttribute("rel = 'enclosure' type = 'application/octet-stream' title = 'LOS vector (tif)' href = %s" % map_los_vec, None)
entry.appendChild(link3)
link4 = root.createElement("link")
link4.setAttribute("rel = 'enclosure' type = 'application/octet-stream' title = 'Browse image (kmz)' href = " + kmz_url, None)
entry.appendChild(link4)


georss = root.createElement("georss:polygon")
georss.setAttribute("xmlns:georss = 'http://www.georss.org/georss'", bounding_box_gml)
entry.appendChild(georss)

box_xmlns = root.createElement("box")
box_xmlns.setAttribute("xmlns = 'http://www.georss.org/georss'", bounding_box_string[14:])
entry.appendChild(box_xmlns)

date = root.createElement("date")
date.setAttribute("xmlns = 'http://purl.org/dc/elements/1.1/'", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
entry.appendChild(date)

spatial = root.createElement("spatial")
spatial.setAttribute("xmlns = 'http://purl.org/dc/terms/'", bounding_box_wkt_string)
entry.appendChild(spatial)

license_child = root.createElement("license")
license_child.setAttribute("", license_string[8:])
entry.appendChild(license_child)

product_type = root.createElement("Product Type")
product_type.setAttribute("xmlns = 'http://www.opengis.net/eop/2.1'", ddss_id_string[9:])
entry.appendChild(product_type)

sar = root.createElement("sar:EarthObservation")
sar.setAttribute("xmlns:eop= 'http://www.opengis.net/eop/2.1'\n \t \t \t \t xmlns:gml = 'http://www.opengis.net/gml/3.2'\n \t \t \t \t xmlns:om = 'http://www.opengis.net/om'\n \t \t \t \t xmlns:ows = 'http://www.opengis.net/ows/2.0'\n \t \t \t \t xmlns:swe = 'http://www.opengis.net/swe/1.0'\n \t \t \t \t xmlns:sar = 'http://www.opengis.net/sar/2.1'\n \t \t \t \t xmlns:xlink = 'http://www.w3.org/1999/xlink'\n \t \t \t \t xmlns:xsi = 'http://www.w3.org/2001/XMLSchema-instance'\n \t \t \t \t xsi:schemaLocation = 'http://www.opengis.net/sar/2.1 xsd/sar.xsd http://www.opengis.net/om  http://schemas.opengis.net/om/2.0/observation.xsd'", None)
entry.appendChild(sar)

phenomenon_time = root.createElement("phenomenonTime")
phenomenon_time.setAttribute("xmlns='http://www.opengis.net/om/2.0'", None)
sar.appendChild(phenomenon_time)

gml_period = root.createElement("gml:TimePeriod")
phenomenon_time.appendChild(gml_period)

# gml_begin = root.createElement("gml:beginPosition")
# gml_begin.setAttribute("", date_of_measurement_start)
# gml_period.appendChild(gml_begin)

# gml_end = root.createElement("gml:endPosition")
# gml_end.setAttribute("", date_of_measurement_end)
# gml_period.appendChild(gml_end)

procedure = root.createElement("procedure")
procedure.setAttribute("xmlns= 'http://www.opengis.net/om/2.0'", None)
sar.appendChild(procedure)

eop_EOE = root.createElement("eop:EarthObservationEquipment")
procedure.appendChild(eop_EOE)

eop_plat = root.createElement("eop:platform")
eop_EOE.appendChild(eop_plat)

eop_PLAT = root.createElement("eop:Platform")
eop_plat.appendChild(eop_PLAT)

eop_name = root.createElement("eop:shortName")
eop_name.setAttribute("", sensor_string[7:])
eop_PLAT.appendChild(eop_name)

eop_inst = root.createElement("eop:instrument")
procedure.appendChild(eop_inst)

eop_INST = root.createElement("eop:Instrument")
eop_inst.appendChild(eop_INST)

eop_shortname = root.createElement("eop:shortName")
eop_shortname.setAttribute("", "SAR")
eop_INST.appendChild(eop_shortname)

eop_desc = root.createElement("eop:description")
eop_desc.setAttribute("", "Synthetic Aperture Radar")
eop_INST.appendChild(eop_desc)

eop_sen = root.createElement("eop:sensor")
procedure.appendChild(eop_sen)

eop_SEN = root.createElement("eop:Sensor")
eop_sen.appendChild(eop_SEN)

eop_sentype = root.createElement("eop:sensorType")
eop_sentype.setAttribute("", "RADAR")
eop_SEN.appendChild(eop_sentype)

eop_OM = root.createElement("eop:operationalMode")
eop_OM.setAttribute("", mode_string[5:])
eop_SEN.appendChild(eop_OM)

eop_wave = root.createElement("eop:wavelengthInformation")
eop_SEN.appendChild(eop_wave)

eop_WAVE = root.createElement("eop:WavelengthInformation")
eop_wave.appendChild(eop_WAVE)

eop_disc = root.createElement("eop:discreteWavelength")
eop_disc.setAttribute("uom = m", wavelength_string[11:])

eop_acqpar = root.createElement("eop:acquisitionParameters")
procedure.appendChild(eop_acqpar)

sar_acq = root.createElement("sar:Acquisition")
eop_acqpar.appendChild(sar_acq)

sar_look = root.createElement("sar:antennaLookDirection")
sar_look.setAttribute("", antenna_side_string[13:])
sar_acq.appendChild(sar_look)

foi = root.createElement("featureOfInterest")
foi.setAttribute("xmlns='http://www.opengis.net/om/2.0'", None)
sar.appendChild(foi)

eop_fp = root.createElement("eop_Footprint")
foi.appendChild(eop_fp)

eop_multiexof = root.createElement("multiExtentOf")
eop_fp.appendChild(eop_multiexof)

gml_multisur = root.createElement("gml:MultiSurface")
eop_multiexof.appendChild(gml_multisur)

gml_surmem = root.createElement("gml:surfaceMembers")
gml_multisur.appendChild(gml_surmem)

gml_pol = root.createElement("gml:Polygon")
gml_surmem.appendChild(gml_pol)

gml_ext = root.createElement("gml:exterior")
gml_pol.appendChild(gml_ext)

gml_lr = root.createElement("gml:LinearRing")
gml_ext.appendChild(gml_lr)

gml_poslist = root.createElement("gml:posList")
gml_poslist.setAttribute("count = 5", bounding_box_gml)
gml_lr.appendChild(gml_poslist)

eop_metprop = root.createElement("eop:metaDataProperty")
sar.appendChild(eop_metprop)

eop_EOM = root.createElement("eop:EarthObservationMetaData")
eop_metprop.appendChild(eop_EOM)

eop_id = root.createElement("eop:identifier")
eop_id.setAttribute("", product_id)
eop_EOM.appendChild(eop_id)

eop_pt = root.createElement("eop:productType")
eop_pt.setAttribute("", ddss_id_string[9:])
eop_EOM.appendChild(eop_pt)

eop_stat = root.createElement("eop:status")
eop_stat.setAttribute("", "ARCHIVED")
eop_EOM.appendChild(eop_stat)

eop_statsubt = root.createElement("eop:statusSubType")
eop_statsubt.setAttribute("", "ONLINE")
eop_EOM.appendChild(eop_statsubt)

eop_doi = root.createElement("eop:doi")
eop_doi.setAttribute("", main_reference_string[16:])
eop_EOM.appendChild(eop_doi)

eop_prodqs = root.createElement("eop:productQualityStatus")
eop_prodqs.setAttribute("", "NOMINAL")
eop_EOM.appendChild(eop_prodqs)

eop_processing = root.createElement("eop:processing")
eop_EOM.appendChild(eop_processing)

eop_pi = root.createElement("eop:ProcessingInformation")
eop_processing.appendChild(eop_pi)

eop_pc = root.createElement("eop:processingCenter")
eop_pc.setAttribute("", user_id_string[9:])
eop_pi.appendChild(eop_pc)

eop_met = root.createElement("eop:method")
eop_met.setAttribute("", applied_algorithm_description_string[20:])
eop_met.setAttribute("Used DEM", used_dem_string[9])
eop_pi.appendChild(eop_met)

eop_soft = root.createElement("eop:processorVersion")
eop_soft.setAttribute("", software_version_string[17:])
eop_pi.appendChild(eop_soft)

# eop_procdate = root.createElement("eop:processingDate")
# eop_procdate.setAttribute("", date_of_measurement)
# eop_pi.appendChild(eop_procdate)

eop_procname = root.createElement("eop:processorName")
eop_procname.setAttribute("", service_used_for_generation_string[28:])
eop_pi.appendChild(eop_procname)

eop_proclev = root.createElement("eop:processingLevel")
eop_proclev.setAttribute("", "EPOS_L1")
eop_pi.appendChild(eop_proclev)

eop_natfor = root.createElement("eop:nativeProductFormat")
eop_natfor.setAttribute("", product_format_string[16:])
eop_pi.appendChild(eop_natfor)

eop_aux = root.createElement("eop:auxiliaryDataSetFileName")
eop_aux.setAttribute("", map_los_vec_string[18:])
eop_pi.appendChild(eop_aux)

# eop_aux1 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux1.setAttribute("", Lookup_table_from_radar_coordinates_to_ground_coordinates)
# eop_pi.appendChild(eop_aux1)

# eop_aux2 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux2.setAttribute("", DEM_radar_geometry)
# eop_pi.appendChild(eop_aux2)

# eop_aux3 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux3.setAttribute("", APS_from_global_model)
# eop_pi.appendChild(eop_aux3)

eop_ven = root.createElement("eop:vendorSpecific")
eop_EOM.appendChild(eop_ven)

eop_spec = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec)

eop_locat = root.createElement("eop:localAttribute")
eop_locat.setAttribute("", "xGroundSPatialRwsolution")
eop_spec.appendChild(eop_locat)

eop_loval = root.createElement("eop:localValue")
eop_loval.setAttribute("", "N/a")
eop_spec.appendChild(eop_loval)

eop_spec1 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec1)

eop_locat1 = root.createElement("eop:localAttribute")
eop_locat1.setAttribute("", "yGroundSPatialRwsolution")
eop_spec1.appendChild(eop_locat1)

eop_loval1 = root.createElement("eop:localValue")
eop_loval1.setAttribute("", "N/a")
eop_spec1.appendChild(eop_loval1)

eop_spec2 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec2)

# eop_locat2 = root.createElement("eop:localAttribute")
# eop_locat2.setAttribute("", "PerpendicularBaseline_[m]")
# eop_spec2.appendChild(eop_locat2)

# eop_loval2 = root.createElement("eop:localValue")
# eop_loval2.setAttribute("", Perpendicular_baseline)
# eop_spec2.appendChild(eop_loval2)

# eop_spec3 = root.createElement("eop:SpecificInformation")
# eop_ven.appendChild(eop_spec3)

# eop_locat3 = root.createElement("eop:localAttribute")
# eop_locat3.setAttribute("", "ParallelBaseline_[m]")
# eop_spec3.appendChild(eop_locat3)

# eop_loval3 = root.createElement("eop:localValue")
# eop_loval3.setAttribute("", Parallel_baseline)
# eop_spec3.appendChild(eop_loval3)

# eop_spec4 = root.createElement("eop:SpecificInformation")
# eop_ven.appendChild(eop_spec4)

# eop_locat4 = root.createElement("eop:localAttribute")
# eop_locat4.setAttribute("", "Along_trak_baseline_[m]")
# eop_spec4.appendChild(eop_locat4)

# eop_loval4 = root.createElement("eop:localValue")
# eop_loval4.setAttribute("", Along_track_baseline)
# eop_spec4.appendChild(eop_loval4)

eop_spec5 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec5)

eop_locat5 = root.createElement("eop:localAttribute")
eop_locat5.setAttribute("", "Reference Point")
eop_spec5.appendChild(eop_locat5)

eop_loval5 = root.createElement("eop:localValue")
eop_loval5.setAttribute("", reference_point_string[16:])
eop_spec5.appendChild(eop_loval5)

eop_spec6 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec6)

eop_locat6 = root.createElement("eop:localAttribute")
eop_locat6.setAttribute("", "Legend_url")
eop_spec6.appendChild(eop_locat6)

eop_loval6 = root.createElement("eop:localValue")
eop_loval6.setAttribute("", legend_url)
eop_spec6.appendChild(eop_loval6)

result = root.createElement("result")
result.setAttribute("xmlns= 'http://www.opengis.net/om/2.0'", None)
sar.appendChild(result)

eop_eor = root.createElement("eop:EarthObservationResult")
result.appendChild(eop_eor)

eop_product = root.createElement("eop:product")
eop_eor.appendChild(eop_product)

eop_prodinf = root.createElement("eop:ProductInformation")
eop_product.appendChild(eop_prodinf)

eop_filename = root.createElement("eop:fileName")
eop_prodinf.appendChild(eop_filename)

ows_sr = root.createElement("ows:ServiceReference")
ows_sr.setAttribute("xlink:href=" + product_url, None)
eop_filename.appendChild(ows_sr)

ows_req = root.createElement("ows:RequestMessage")
ows_sr.appendChild(ows_req)

eop_size = root.createElement("eop:size")
eop_size.setAttribute("uom = 'Byte", product_size_string[14:])
eop_prodinf.appendChild(eop_size)

# eop_rsi = root.createElement("eop:referenceSystemIdentifier")
# eop_rsi.setAttribute("codeSpace = " + code_space, code_value)

eop_browse = root.createElement("eop:browse")
eop_eor.appendChild(eop_browse)

eop_browseinf = root.createElement("eop:BrowseInformation")
eop_browse.appendChild(eop_browseinf)

eop_filename2 = root.createElement("eop:fileName")
eop_browseinf.appendChild(eop_filename2)

ows_seref = root.createElement("ows:ServiceReference")
ows_seref.setAttribute("xlink:href =" + preview_url, None)
eop_filename2.appendChild(ows_seref)

ow_reqmes = root.createElement("ows:RequestMessage")
ows_seref.appendChild(ow_reqmes)

eop_refsysid = root.createElement("eop:ReferenceSystemIdentifier")
# eop_refsysid.setAttribute("codeSpace = " + code_space, code_value)
# eop_browseinf.appendChild(eop_refsysid)

eop_type = root.createElement("eop:type")
eop_type.setAttribute("", "QUICKLOOK")
eop_browseinf.appendChild(eop_type)

eop_parameter = root.createElement("eop:parameter")
eop_eor.appendChild(eop_parameter)

eop_parinf = root.createElement("eop:ParameterInformation")
eop_parameter.appendChild(eop_parinf)

eop_uom = root.createElement("eop:unitOfMeasure")
eop_uom.setAttribute("uom =" + value_unit_string[11:], None)
eop_parinf.appendChild(eop_uom)

eop_phen = root.createElement("eop:phenomenon")
eop_parinf.appendChild(eop_phen)

swe_phen = root.createElement("swe:Phenomenon")
swe_phen.setAttribute("xmlns = 'http://www.opengis.net/gml' ns1:id='phenom1'", None)
eop_phen.appendChild(swe_phen)

name = root.createElement("name")
name.setAttribute("xmlns = 'http://www.opengis.net/gml'", "Displacement")
swe_phen.appendChild(name)

xml_str = root.toprettyxml(indent="\t")

with open("ew.xml", "w") as xml_file:
    root.writexml(xml_file, indent = "\n \t", addindent= "\t")

print("East-West GeoTIFF ready!")
print("\n.\n.\n.")
print("Writing of Up GeoTIFF...")
# Repeat for UP
# 1) Extension of the file
extent = [aoi_min_lon, aoi_min_lat, aoi_max_lon, aoi_max_lat]

# 2) Get GDAL driver GeoTiff
driver = gdal.GetDriverByName("GTiff")
data_type = gdal.GDT_Float32

# 3) Create a temporary grid
grid_data = driver.Create("grid_data", dim_range_lon, dim_range_lat, 1, data_type) 

# 4) Write data bands
grid_data.GetRasterBand(1).WriteArray(matrixUp)

# 5) Reference System
rs = osr.SpatialReference()
rs.ImportFromEPSG(4326)
grid_data.SetProjection(rs.ExportToWkt())
grid_data.SetGeoTransform(getGeoTransform(extent, dim_range_lat, dim_range_lon))

# 6) Set the name of the file and save the file
file_name_matrixUp = "UP_.tif"
driver.CreateCopy(file_name_matrixUp, grid_data, 0)

# 7) Close the file
driver = None
grid_data = None
 
# 8) Delete the temporary grid
import os                
os.remove("grid_data")

# Metadata writing 
print("Metadata creation -- WORK IN PROGRESS PART, IF ANY ERROR OCCURS, COMMENT FROM LINE 1159 TO LINE 1898")
ddss_id = "UP"
product_id = "UP_XXX"
preview_url = "xxx"
legend_url = "xxx"
product_url = "xxx"
kmz_url = "xxx"
Lookup_table_from_radar_coordinates_to_ground_coordinates = "xxx"
code_space = "xxx"
code_value = "xxx"



for x in os.listdir():
    if "InU" in x:
        inu_data = os.listdir(x)
        for y in inu_data:
            if ".metadata" in y:
                metadata_path = (str(x) + "/" + y)    
                metadata = open(metadata_path, "r")
                for line in metadata:
                    if "Sensor: " in line:
                        sensor_line = line.strip()
                        sensor = sensor_line[8:]
                    if "License: " in line:
                        license_line = line.strip()
                        license = license_line[9:]
                    if "User_ID: " in line:
                        user_id_line = line.strip()
                        user_id = user_id_line[9:]
                    if "Software_version: " in line:
                        software_version_line = line.strip()
                        software_version = software_version_line[18:]
                    if "Geographic_CS_type_code: " in line:
                        reference_system_line = line.strip()
                        reference_system = reference_system_line[25:]
                    if "Applied_algorithm_description: " in line:
                        applied_algorithm_description_line = line.strip()
                        applied_algorithm_description = applied_algorithm_description_line[31:]
                    if "Main_reference: " in line:
                        main_reference_line = line.strip()
                        main_reference = main_reference_line[17:]
                    if "Service_used_for_generation: " in line:
                        service_used_for_generation_line = line.strip()
                        service_generation = service_used_for_generation_line[39:]
                    if "Used_DEM: " in line:
                        used_dem_line = line.strip()
                        used_dem = used_dem_line[10:]
                    if "Applied_unwrapping_algorithm: " in line:
                        applied_unwrapping_algorithm_line = line.strip()
                        applied_unwrapping_algorithm = applied_unwrapping_algorithm_line[30:]
                    if "Mode: " in line:
                        mode_line = line.strip()
                        mode = mode_line[6:]
                    if "Antenna_side: " in line:
                        antenna_side_line = line.strip()
                        antenna_side = antenna_side_line[14:]
                    if "Wavelength: " in line:
                        wavelength_line = line.strip()
                        wavelength = wavelength_line[11:]
                    if "Value_unit: " in line:
                        value_unit_line = line.strip()
                        value_unit = value_unit_line[12:]
                    if "AffiliationIdentifier: " in line:
                        affiliationidentifier_line = line.strip()
                        affiliationidentifier = affiliationidentifier_line[23:]

                    
ddss_id_string = "DDSS_ID: " + ddss_id
product_id_string = "Product_ID: " + product_id
product_format_string = "Product_format: " + product_format
product_size_string = "Product_size: " + str(product_size)
preview_url_string = "Preview_url: " + preview_url
legend_url_string = "Legend_url: " + legend_url 
product_url_string = "Product_url: " + product_url 
bounding_box_string = "Bounding_box: " + bounding_box
bounding_box_wkt_string = "Bounding_box_wkt: POLYGON((" + bounding_box_wkt
license_string = "License: " + license
user_id_string = "User_ID: " + user_id
software_version_string = "Software_version: " + software_version
applied_algorithm_description_string = "Applied_algorithm_description: " + applied_algorithm_description 
main_reference_string = "Main_reference: " + main_reference 
# date_of_measurement_start_string = "Date_of_measurement_start: " + date_of_measurement_start
# date_of_measurement_end_string = "Date_of_measurement_end: " + date_of_measurement_end
date_of_production_string = "Date_of_production: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
service_used_for_generation_string = "Service_used_for_generation: " + service_generation
reference_system_string = "Geographic_CS_type_code: " + reference_system
used_dem_string = "Used_DEM: " + used_dem
perpendicular_baseline_string = "Perpendicular_baseline: N/a"
parallel_baseline_string = "Parallel_baseline: N/a"
along_track_baseline_string = "Along_track_baseline: N/a"
map_los_vec_string = "Map_of_LOS_vector: " + str(map_los_vec)
applied_unwrapping_algorithm_string = "Applied_unwrapping_algorithm: " + applied_unwrapping_algorithm
reference_point_string = "Reference_point: " + reference_point
spatial_resolution_string = "Spatial_resolution: "  + spatial_resolution #30 m
sensor_string = "Sensor: " + sensor
mode_string = "Mode: " + mode
antenna_side_string = "Antenna_side: " + antenna_side
relative_orbit_number_string = "Relative_orbit_number: " + str(relative_orbit_number) # serie di orbite combinate se ci sono nei metadati delle linee di vista (parametrizzato su n input data)
wavelength_string = "Wavelength: " + wavelength
value_unit_string = "Value_unit: " + value_unit 
number_of_looks_azimuth_string = "Number_of_looks_azimuth: N/a"
number_of_looks_range_string = "Number_of_looks_range: N/a" 
applied_filter_string = "Applied_filter: No_Filter"
kmz_url_string = "kmz_url: " + kmz_url
# affiliationidentifier_string = "AffiliationIdentifier: " + affiliationidentifier
method_string = "MethodDescription: " + method

metadata = [ddss_id_string, product_id_string, product_format_string, product_size_string, preview_url_string, legend_url_string, product_url_string, bounding_box_string, bounding_box_wkt_string, license_string, user_id_string, software_version_string, applied_algorithm_description_string,  main_reference_string, date_of_production_string, service_used_for_generation_string, reference_system_string, used_dem_string, perpendicular_baseline_string, parallel_baseline_string, along_track_baseline_string, map_los_vec_string, applied_unwrapping_algorithm_string, reference_point_string, spatial_resolution_string, sensor_string, mode_string, antenna_side_string, relative_orbit_number_string, wavelength_string, value_unit_string, number_of_looks_azimuth_string, number_of_looks_range_string, applied_filter_string, kmz_url_string, method_string]

metadata_file_path = folder + "UP.metadata"

with open(metadata_file_path, "w") as metadata_file:
    for x in metadata:
        metadata_file.write(x + "\n")



root = minidom.Document()
xml = root.createElement("feed")
xml.setAttribute("xmlns", "http://www.w3.org/2005/Atom")
root.appendChild(xml)

title_child = root.createElement("title")
title_child.setAttribute("type= 'text'", "Template for ingestion of an EPOS product")
xml.appendChild(title_child)

id = root.createElement("id")
id.setAttribute("", "http://catalog.terradue.com/gep-epos/")
xml.appendChild(id)

entry = root.createElement("entry")
xml.appendChild(entry)

identifier = root.createElement("identifier")
identifier.setAttribute("xmlsn = http://purl.org/dc/elements/1.1/", product_id)
entry.appendChild(identifier)

title_2 = root.createElement("title")
title_2.setAttribute("type= 'text'", product_id)
entry.appendChild(title_2)

summary = root.createElement("summary type= 'html'")
entry.appendChild(summary)

data = root.createElement("![CDATA[")
summary.appendChild(data)

table = root.createElement("table")
summary.appendChild(table)

body = root.createElement("tbody")
table.appendChild(body)

tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
tr.appendChild(td)

a = root.createElement("a")
a.setAttribute("href=" + preview_url + " target= '_blank' title = 'View preview image'", None)
td.appendChild(a)

img = root.createElement("img")
img.setAttribute("align= 'left' border= '0' width= '100px' title= 'View legend' src=" + preview_url, None)
a.appendChild(img)

td = root.createElement("td")
tr.appendChild(td)

a = root.createElement("a")
a.setAttribute("href=" + legend_url + " target= '_blank' title = 'View preview image'", None)
td.appendChild(a)

img = root.createElement("img")
img.setAttribute("align= 'left' border= '0' width= '100px' title= 'View legend' src=" + legend_url, None)
a.appendChild(img)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Product Identifier")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", product_id)
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Platform")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", sensor_string[7:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Sensor Type")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", "RADAR")
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Product Type")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", ddss_id_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Processing Center")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", user_id_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Date of measurement start")
tr.appendChild(td)

# strong = root.createElement("strong")
# strong.setAttribute("", date_of_measurement_start)
# td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Date of measurement end")
tr.appendChild(td)

# strong = root.createElement("strong")
# strong.setAttribute("", date_of_measurement_end)
# td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Processor Name")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", service_used_for_generation_string[28:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Unit of measure")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", value_unit_string[11:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Referemce System")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", reference_system_string[24:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "License")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", license_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Reference point [lon lat]")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", reference_point_string[16:])
td.appendChild(strong)

tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "DEM")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", used_dem_string[9:])
td.appendChild(strong)


tr = root.createElement("tr")
body.appendChild(tr)

td = root.createElement("td")
td.setAttribute("", "Method")
tr.appendChild(td)

strong = root.createElement("strong")
strong.setAttribute("", applied_algorithm_description[20:])
td.appendChild(strong)


link1 = root.createElement("link")
link1.setAttribute("href =" + product_url + "rel= 'enclosure' title= 'Product URL' type='text/html'", None)
entry.appendChild(link1)
link2 = root.createElement("link")
link2.setAttribute("href =" + preview_url + "rel ='icon' title = 'Preview URL' type ='image/x-icon'", None)
entry.appendChild(link2)
link3 = root.createElement("link")
link3.setAttribute("rel = 'enclosure' type = 'application/octet-stream' title = 'LOS vector (tif)' href = %s" % map_los_vec, None)
entry.appendChild(link3)
link4 = root.createElement("link")
link4.setAttribute("rel = 'enclosure' type = 'application/octet-stream' title = 'Browse image (kmz)' href = " + kmz_url, None)
entry.appendChild(link4)


georss = root.createElement("georss:polygon")
georss.setAttribute("xmlns:georss = 'http://www.georss.org/georss'", bounding_box_gml)
entry.appendChild(georss)

box_xmlns = root.createElement("box")
box_xmlns.setAttribute("xmlns = 'http://www.georss.org/georss'", bounding_box_string[14:])
entry.appendChild(box_xmlns)

date = root.createElement("date")
date.setAttribute("xmlns = 'http://purl.org/dc/elements/1.1/'", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
entry.appendChild(date)

spatial = root.createElement("spatial")
spatial.setAttribute("xmlns = 'http://purl.org/dc/terms/'", bounding_box_wkt_string)
entry.appendChild(spatial)

license_child = root.createElement("license")
license_child.setAttribute("", license_string[8:])
entry.appendChild(license_child)

product_type = root.createElement("Product Type")
product_type.setAttribute("xmlns = 'http://www.opengis.net/eop/2.1'", ddss_id_string[9:])
entry.appendChild(product_type)

sar = root.createElement("sar:EarthObservation")
sar.setAttribute("xmlns:eop= 'http://www.opengis.net/eop/2.1'\n \t \t \t \t xmlns:gml = 'http://www.opengis.net/gml/3.2'\n \t \t \t \t xmlns:om = 'http://www.opengis.net/om'\n \t \t \t \t xmlns:ows = 'http://www.opengis.net/ows/2.0'\n \t \t \t \t xmlns:swe = 'http://www.opengis.net/swe/1.0'\n \t \t \t \t xmlns:sar = 'http://www.opengis.net/sar/2.1'\n \t \t \t \t xmlns:xlink = 'http://www.w3.org/1999/xlink'\n \t \t \t \t xmlns:xsi = 'http://www.w3.org/2001/XMLSchema-instance'\n \t \t \t \t xsi:schemaLocation = 'http://www.opengis.net/sar/2.1 xsd/sar.xsd http://www.opengis.net/om  http://schemas.opengis.net/om/2.0/observation.xsd'", None)
entry.appendChild(sar)

phenomenon_time = root.createElement("phenomenonTime")
phenomenon_time.setAttribute("xmlns='http://www.opengis.net/om/2.0'", None)
sar.appendChild(phenomenon_time)

gml_period = root.createElement("gml:TimePeriod")
phenomenon_time.appendChild(gml_period)

# gml_begin = root.createElement("gml:beginPosition")
# gml_begin.setAttribute("", date_of_measurement_start)
# gml_period.appendChild(gml_begin)

# gml_end = root.createElement("gml:endPosition")
# gml_end.setAttribute("", date_of_measurement_end)
# gml_period.appendChild(gml_end)

procedure = root.createElement("procedure")
procedure.setAttribute("xmlns= 'http://www.opengis.net/om/2.0'", None)
sar.appendChild(procedure)

eop_EOE = root.createElement("eop:EarthObservationEquipment")
procedure.appendChild(eop_EOE)

eop_plat = root.createElement("eop:platform")
eop_EOE.appendChild(eop_plat)

eop_PLAT = root.createElement("eop:Platform")
eop_plat.appendChild(eop_PLAT)

eop_name = root.createElement("eop:shortName")
eop_name.setAttribute("", sensor_string[7:])
eop_PLAT.appendChild(eop_name)

eop_inst = root.createElement("eop:instrument")
procedure.appendChild(eop_inst)

eop_INST = root.createElement("eop:Instrument")
eop_inst.appendChild(eop_INST)

eop_shortname = root.createElement("eop:shortName")
eop_shortname.setAttribute("", "SAR")
eop_INST.appendChild(eop_shortname)

eop_desc = root.createElement("eop:description")
eop_desc.setAttribute("", "Synthetic Aperture Radar")
eop_INST.appendChild(eop_desc)

eop_sen = root.createElement("eop:sensor")
procedure.appendChild(eop_sen)

eop_SEN = root.createElement("eop:Sensor")
eop_sen.appendChild(eop_SEN)

eop_sentype = root.createElement("eop:sensorType")
eop_sentype.setAttribute("", "RADAR")
eop_SEN.appendChild(eop_sentype)

eop_OM = root.createElement("eop:operationalMode")
eop_OM.setAttribute("", mode_string[5:])
eop_SEN.appendChild(eop_OM)

eop_wave = root.createElement("eop:wavelengthInformation")
eop_SEN.appendChild(eop_wave)

eop_WAVE = root.createElement("eop:WavelengthInformation")
eop_wave.appendChild(eop_WAVE)

eop_disc = root.createElement("eop:discreteWavelength")
eop_disc.setAttribute("uom = m", wavelength_string[11:])

eop_acqpar = root.createElement("eop:acquisitionParameters")
procedure.appendChild(eop_acqpar)

sar_acq = root.createElement("sar:Acquisition")
eop_acqpar.appendChild(sar_acq)

sar_look = root.createElement("sar:antennaLookDirection")
sar_look.setAttribute("", antenna_side_string[13:])
sar_acq.appendChild(sar_look)

foi = root.createElement("featureOfInterest")
foi.setAttribute("xmlns='http://www.opengis.net/om/2.0'", None)
sar.appendChild(foi)

eop_fp = root.createElement("eop_Footprint")
foi.appendChild(eop_fp)

eop_multiexof = root.createElement("multiExtentOf")
eop_fp.appendChild(eop_multiexof)

gml_multisur = root.createElement("gml:MultiSurface")
eop_multiexof.appendChild(gml_multisur)

gml_surmem = root.createElement("gml:surfaceMembers")
gml_multisur.appendChild(gml_surmem)

gml_pol = root.createElement("gml:Polygon")
gml_surmem.appendChild(gml_pol)

gml_ext = root.createElement("gml:exterior")
gml_pol.appendChild(gml_ext)

gml_lr = root.createElement("gml:LinearRing")
gml_ext.appendChild(gml_lr)

gml_poslist = root.createElement("gml:posList")
gml_poslist.setAttribute("count = 5", bounding_box_gml)
gml_lr.appendChild(gml_poslist)

eop_metprop = root.createElement("eop:metaDataProperty")
sar.appendChild(eop_metprop)

eop_EOM = root.createElement("eop:EarthObservationMetaData")
eop_metprop.appendChild(eop_EOM)

eop_id = root.createElement("eop:identifier")
eop_id.setAttribute("", product_id)
eop_EOM.appendChild(eop_id)

eop_pt = root.createElement("eop:productType")
eop_pt.setAttribute("", ddss_id_string[9:])
eop_EOM.appendChild(eop_pt)

eop_stat = root.createElement("eop:status")
eop_stat.setAttribute("", "ARCHIVED")
eop_EOM.appendChild(eop_stat)

eop_statsubt = root.createElement("eop:statusSubType")
eop_statsubt.setAttribute("", "ONLINE")
eop_EOM.appendChild(eop_statsubt)

eop_doi = root.createElement("eop:doi")
eop_doi.setAttribute("", main_reference_string[16:])
eop_EOM.appendChild(eop_doi)

eop_prodqs = root.createElement("eop:productQualityStatus")
eop_prodqs.setAttribute("", "NOMINAL")
eop_EOM.appendChild(eop_prodqs)

eop_processing = root.createElement("eop:processing")
eop_EOM.appendChild(eop_processing)

eop_pi = root.createElement("eop:ProcessingInformation")
eop_processing.appendChild(eop_pi)

eop_pc = root.createElement("eop:processingCenter")
eop_pc.setAttribute("", user_id_string[9:])
eop_pi.appendChild(eop_pc)

eop_met = root.createElement("eop:method")
eop_met.setAttribute("", applied_algorithm_description_string[20:])
eop_met.setAttribute("Used DEM", used_dem_string[9])
eop_pi.appendChild(eop_met)

eop_soft = root.createElement("eop:processorVersion")
eop_soft.setAttribute("", software_version_string[17:])
eop_pi.appendChild(eop_soft)

# eop_procdate = root.createElement("eop:processingDate")
# eop_procdate.setAttribute("", date_of_measurement)
# eop_pi.appendChild(eop_procdate)

eop_procname = root.createElement("eop:processorName")
eop_procname.setAttribute("", service_used_for_generation_string[28:])
eop_pi.appendChild(eop_procname)

eop_proclev = root.createElement("eop:processingLevel")
eop_proclev.setAttribute("", "EPOS_L1")
eop_pi.appendChild(eop_proclev)

eop_natfor = root.createElement("eop:nativeProductFormat")
eop_natfor.setAttribute("", product_format_string[16:])
eop_pi.appendChild(eop_natfor)

eop_aux = root.createElement("eop:auxiliaryDataSetFileName")
eop_aux.setAttribute("", map_los_vec_string[18:])
eop_pi.appendChild(eop_aux)

# eop_aux1 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux1.setAttribute("", Lookup_table_from_radar_coordinates_to_ground_coordinates)
# eop_pi.appendChild(eop_aux1)

# eop_aux2 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux2.setAttribute("", DEM_radar_geometry)
# eop_pi.appendChild(eop_aux2)

# eop_aux3 = root.createElement("eop:auxiliaryDataSetFileName")
# eop_aux3.setAttribute("", APS_from_global_model)
# eop_pi.appendChild(eop_aux3)

eop_ven = root.createElement("eop:vendorSpecific")
eop_EOM.appendChild(eop_ven)

eop_spec = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec)

eop_locat = root.createElement("eop:localAttribute")
eop_locat.setAttribute("", "xGroundSPatialRwsolution")
eop_spec.appendChild(eop_locat)

eop_loval = root.createElement("eop:localValue")
eop_loval.setAttribute("", "N/a")
eop_spec.appendChild(eop_loval)

eop_spec1 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec1)

eop_locat1 = root.createElement("eop:localAttribute")
eop_locat1.setAttribute("", "yGroundSPatialRwsolution")
eop_spec1.appendChild(eop_locat1)

eop_loval1 = root.createElement("eop:localValue")
eop_loval1.setAttribute("", "N/a")
eop_spec1.appendChild(eop_loval1)

eop_spec2 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec2)

# eop_locat2 = root.createElement("eop:localAttribute")
# eop_locat2.setAttribute("", "PerpendicularBaseline_[m]")
# eop_spec2.appendChild(eop_locat2)

# eop_loval2 = root.createElement("eop:localValue")
# eop_loval2.setAttribute("", Perpendicular_baseline)
# eop_spec2.appendChild(eop_loval2)

# eop_spec3 = root.createElement("eop:SpecificInformation")
# eop_ven.appendChild(eop_spec3)

# eop_locat3 = root.createElement("eop:localAttribute")
# eop_locat3.setAttribute("", "ParallelBaseline_[m]")
# eop_spec3.appendChild(eop_locat3)

# eop_loval3 = root.createElement("eop:localValue")
# eop_loval3.setAttribute("", Parallel_baseline)
# eop_spec3.appendChild(eop_loval3)

# eop_spec4 = root.createElement("eop:SpecificInformation")
# eop_ven.appendChild(eop_spec4)

# eop_locat4 = root.createElement("eop:localAttribute")
# eop_locat4.setAttribute("", "Along_trak_baseline_[m]")
# eop_spec4.appendChild(eop_locat4)

# eop_loval4 = root.createElement("eop:localValue")
# eop_loval4.setAttribute("", Along_track_baseline)
# eop_spec4.appendChild(eop_loval4)

eop_spec5 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec5)

eop_locat5 = root.createElement("eop:localAttribute")
eop_locat5.setAttribute("", "Reference Point")
eop_spec5.appendChild(eop_locat5)

eop_loval5 = root.createElement("eop:localValue")
eop_loval5.setAttribute("", reference_point_string[16:])
eop_spec5.appendChild(eop_loval5)

eop_spec6 = root.createElement("eop:SpecificInformation")
eop_ven.appendChild(eop_spec6)

eop_locat6 = root.createElement("eop:localAttribute")
eop_locat6.setAttribute("", "Legend_url")
eop_spec6.appendChild(eop_locat6)

eop_loval6 = root.createElement("eop:localValue")
eop_loval6.setAttribute("", legend_url)
eop_spec6.appendChild(eop_loval6)

result = root.createElement("result")
result.setAttribute("xmlns= 'http://www.opengis.net/om/2.0'", None)
sar.appendChild(result)

eop_eor = root.createElement("eop:EarthObservationResult")
result.appendChild(eop_eor)

eop_product = root.createElement("eop:product")
eop_eor.appendChild(eop_product)

eop_prodinf = root.createElement("eop:ProductInformation")
eop_product.appendChild(eop_prodinf)

eop_filename = root.createElement("eop:fileName")
eop_prodinf.appendChild(eop_filename)

ows_sr = root.createElement("ows:ServiceReference")
ows_sr.setAttribute("xlink:href=" + product_url, None)
eop_filename.appendChild(ows_sr)

ows_req = root.createElement("ows:RequestMessage")
ows_sr.appendChild(ows_req)

eop_size = root.createElement("eop:size")
eop_size.setAttribute("uom = 'Byte", product_size_string[14:])
eop_prodinf.appendChild(eop_size)

# eop_rsi = root.createElement("eop:referenceSystemIdentifier")
# eop_rsi.setAttribute("codeSpace = " + code_space, code_value)

eop_browse = root.createElement("eop:browse")
eop_eor.appendChild(eop_browse)

eop_browseinf = root.createElement("eop:BrowseInformation")
eop_browse.appendChild(eop_browseinf)

eop_filename2 = root.createElement("eop:fileName")
eop_browseinf.appendChild(eop_filename2)

ows_seref = root.createElement("ows:ServiceReference")
ows_seref.setAttribute("xlink:href =" + preview_url, None)
eop_filename2.appendChild(ows_seref)

ow_reqmes = root.createElement("ows:RequestMessage")
ows_seref.appendChild(ow_reqmes)

eop_refsysid = root.createElement("eop:ReferenceSystemIdentifier")
# eop_refsysid.setAttribute("codeSpace = " + code_space, code_value)
# eop_browseinf.appendChild(eop_refsysid)

eop_type = root.createElement("eop:type")
eop_type.setAttribute("", "QUICKLOOK")
eop_browseinf.appendChild(eop_type)

eop_parameter = root.createElement("eop:parameter")
eop_eor.appendChild(eop_parameter)

eop_parinf = root.createElement("eop:ParameterInformation")
eop_parameter.appendChild(eop_parinf)

eop_uom = root.createElement("eop:unitOfMeasure")
eop_uom.setAttribute("uom =" + value_unit_string[11:], None)
eop_parinf.appendChild(eop_uom)

eop_phen = root.createElement("eop:phenomenon")
eop_parinf.appendChild(eop_phen)

swe_phen = root.createElement("swe:Phenomenon")
swe_phen.setAttribute("xmlns = 'http://www.opengis.net/gml' ns1:id='phenom1'", None)
eop_phen.appendChild(swe_phen)

name = root.createElement("name")
name.setAttribute("xmlns = 'http://www.opengis.net/gml'", "Displacement")
swe_phen.appendChild(name)

xml_str = root.toprettyxml(indent="\t")

with open("UP.xml", "w") as xml_file:
    root.writexml(xml_file, indent = "\n \t", addindent= "\t")


print("Up GeoTIFF ready!")
print("\n.\n.\n.")
print("***End of the Processing!***")


