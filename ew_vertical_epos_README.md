This README.md will help you on handling the Python scripts "ew_vert" for the extraction of the East-West and Vertical components of the displacement field using EPOSAR database. These scripts must be downloaded, copied and pasted in the Jupyter Notebook environmnets of EPOS analysis section.
As you can notice, two different Python scripts are present in the branch:
1) "ew_vert_epos_complete.py"
2) "ew_vert_epos_kernel.ipynb"
3) "ew_vert_epos_edit.py"

**EW_VERT_EPOS_COMPLETE.PY**

This script must be used when _the input data are inserted manually_, and not by ingesting them directly from EPOS.
In this case, you must create different folders with the names of the products to process, and upload inside them each object that would be contained in a standard EPOSAR folder
This is how the folders should appear

<img width="503" height="192" alt="image" src="https://github.com/user-attachments/assets/33cf6395-3a09-4b94-9074-1aa8d9ab99fe" />

and this is the list of the products inside each folder:

<img width="502" height="115" alt="image" src="https://github.com/user-attachments/assets/678dffdf-5b56-49c5-85e0-76f86b1c5e59" />

_this is a LOS vector maps folder content_

<img width="501" height="245" alt="image" src="https://github.com/user-attachments/assets/454580b3-44b6-43e2-ac07-da012a07db8e" />

_this is a unwrapped interferogram folder content_



Once you organized your data in such way, it's time to process the data with the .py file. But first, some problems with library import have to be solved. For this reason, copy-paste only the lines of the script with the import of the library on a Jupyter notebook kernel, and run to check which library must be imported
To import successfully a library, type the following command:

_pip install library_

and restart the kernel once the process has been concluded. Usually, the import of the library osgeo is the most tricky. Use the command:

_conda install -c conda-forge gdal_

and restart the kernel. After that, the import stage should be finalized and you can copy-paste the entire script. In the variable folder, the "" is referred to the actual folder in which you are located. In you organize the data as reported in the images above, you must not change this variable. In contrary case, if you organize your data in another sub-folder (e.g., "eposar_test"), therefore the variable folder = "eposar_test".
At this stage, the script should be ready for the use. A small recommendation regards the part of Metadata writing, which is work in progress since each metadata element is extrapolated from the input data. Some of the elements, in fact, are left empty for this reason, and they are planned to be defined as soon as possible. Some other elements, since are derived from the input data, are created using strings manipulation. If some error occurs during the writing of these elements, please just comment the metadata creation lines. In the script, they are located at lines 384 and 1159. We apologize for this inconvenient.

**EW_VERT_EPOS_KERNEL.IPYNB**

This file consists of the transposition of the file _ew_vert_epos_complete.py_ in a Jupyter Notebook kernel ready to use. The script is organized in cells of debugging depending on the function to be performed (e.g., data organization, unwrapped interferogram processing, CosNEU file processing...). 

**EW_VERT_EPOS_EDIT.PY**

This script is compatible with the import of the data from EPOS portal, without uploading them manually. The data are stored in the folder "data/staginghistory/", and then you have to check manually in which folder the data is stores (e.g., stage.0001)
(_script actually work in progress for permission issues in extracting zip files_)
