**EAST-WEST AND VERTICAL DISPLACEMENT FIELD EXTRAPOLATION USING EPOSAR PRODUCTS**
**by Occhipinti M., De Luca C., Manunta M., Monterroso F., Casu F.**
**d**

This README.md will help you on handling the Python scripts for the extraction of the East-West and Vertical components of the displacement field using EPOSAR database in EPOS Analysis platform. The folder _components_combination_ contains different files: Python files (.py) with the "pure" Python scripts for the computation of EW and Vertical displacement outside EPOS Analysis platform; Jupyter Notebook files (.ipynb) containing the ready-to-use kernel for the performance of the script in EPOS Analysis platform. Each script is performable also outside EPOS Analysis portal by using any type of interpreter; the single condition to respect is the usage of EPOSAR products and their organization in the proper folders. This aspect will be further discussed at point 1.

1) Input data organization
2) "ew_vert_epos_complete.py"
3) "ew_vert_epos_kernel.ipynb"
4) "ew_vert_epos_edit.py"
5) "requirements.txt"

1. **INPUT DATA ORGANIZATION**

The most simple example for input data consists of one ascending and one descending acquisition. Each acquisition will consists of: _unwrapped interferograms_ (InU files); _map of LOS vectors_ (CosNEU). The script will automatically recognize the mutual InU-CosNEU association for the correct performance of the processing, and for this reason, it is very important to respect the name of the folders when downloading them without changing any name (download available here: https://www.ics-c.epos-eu.org/). This means that each InU file folder MUST have the name of the respective InU file (e.g., InU_CNRIREA_123G folder will contain InU_CNRIREA_123G.tif, InU_CNRIREA_123G.xml, InU_CNRIREA_123G.metadata ...), as for the CosNEU folder (e.g., CosNEU_CNRIREA_123G folder will contain CosNEU_CNRIREA_123G.tif, CosNEU_CNRIREA_123G.metadata ...).

IMPORTANT: **if you are planning to perform the Python script outside EPOS analysis portal**, all the folders (InU and Cos) must be stored in one single macro-folder, which path must be inserted in the variable "folder" at line 27 of the Python script (e.g., C:/navre/personal/epos_test/).

**if you're performing the Python script outside EPOS analysis portal**, be sure to insert the proper path in "folder" variable (e.g., folder  "C:/navre/personal/epos_test").

**If you're performing the Python script in EPOS analysis portal and you have manually uploaded these testing data in your environment**, please upload each InU and Cos folder in the main branch of your environment, since "folder" variable must be empty. More information on how to organize your input data will be given in 2. section
   
2. **EW_VERT_EPOS_COMPLETE.PY**

This script must be used when _the script is run outside EPOS Analysis tool_.
When using this .py file and when uploading manually the input data, you must manually create different folders with the names of the products to process, and therefore, upload manually each data in the respective folder. 

This is how the folders should appear.

<img width="503" height="192" alt="image" src="https://github.com/user-attachments/assets/33cf6395-3a09-4b94-9074-1aa8d9ab99fe" />

and this is the list of the products inside each folder:

<img width="502" height="115" alt="image" src="https://github.com/user-attachments/assets/678dffdf-5b56-49c5-85e0-76f86b1c5e59" />

_this is a LOS vector maps folder content_

<img width="501" height="245" alt="image" src="https://github.com/user-attachments/assets/454580b3-44b6-43e2-ac07-da012a07db8e" />

_this is a unwrapped interferogram folder content_


Once you organized your data in such way, it's time to process the data with the .py file. Please store all the folders in a single macro folder.
In the variable folder, the "" is referred to the actual folder in which you are located. If you organize your data in a macro folder (e.g., "eposar_test"), therefore the variable folder = "eposar_test".

A small recommendation regards the part of Metadata writing, which is work in progress since each metadata element is extrapolated from the input data. Some of the elements, in fact, are left empty for this reason, and they are planned to be defined as soon as possible. Some other elements, since are derived from the input data, are created using strings manipulation. If some error occurs during the writing of these elements, please just comment the metadata creation lines. In the script, they are located at lines 384 and 1159. We apologize for this inconvenient.

3. **EW_VERT_EPOS_KERNEL.IPYNB**

This file consists of the transposition of the file _ew_vert_epos_complete.py_ in a Jupyter Notebook kernel ready to use in EPOS Analysis platform. it's advisable to use this file for the performance of the software due to the completeness of the information furnished. The script is organized in cells of debugging depending on the function to be performed (e.g., data organization, unwrapped interferogram processing, CosNEU file processing...). Also in this case, beware to the issues that may occurr when importing the libraries. Please refer to the tips reported in 2. section.

4. **REQUIREMENTS.TXT**
List of the necessary modules for the running of the script.
