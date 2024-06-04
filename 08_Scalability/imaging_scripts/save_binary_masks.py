# %%
# import glob
import os
import sys
import numpy as np

import czifile as czi
from skimage.filters import threshold_otsu
from skimage import feature
import numpy as np

from skimage import morphology, measure, segmentation, feature
from scipy import ndimage as ndi

from skimage import graph
from skimage import img_as_float

from skimage.morphology import (erosion, dilation, opening, closing)
from skimage.morphology import black_tophat

from skimage.morphology import remove_small_holes
from skimage.morphology import disk  

# %%
def save_object(obj, filename):

    import pickle

    try:
        with open(filename + ".pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)


# %%

footprint = morphology.disk(30)

imgs_dict = {}
folder = sys.argv[1]
print(folder)

if 'Imaging_to_confluency' in folder:
    folder_date_time = folder.split('/')[-4] + '_t' + folder.split('/')[-3] + '+confluency'
elif 'Generation' in folder:
    folder_date_time = folder.split('/')[-4] + '_t' + folder.split('/')[-3] + '+generation'
else:
    folder_date_time = folder.split('/')[-3] + '_t' + folder.split('/')[-2]
print(folder_date_time)

# %%
files = [i for i in os.listdir(folder) if i.endswith('.czi')]
imgs_dict['file_names'] = files

nimgs = len(files)

img = czi.CziFile(folder + '/' + files[0]).asarray()[0]
img = np.squeeze(img)
height, width = img.shape

#img_array = np.zeros((nimgs, height, width))
labels_array = np.zeros((nimgs, height, width))


# %%
i = -1

for file in files:    
    if file.endswith('.czi'):
        i += 1
        print(file)
        file_name = file.strip('.czi')

        img = czi.CziFile(folder + '/' + file).asarray()[0]
        img = np.squeeze(img)
        
        #img_array[i,:,:] = img
        
        img = img_as_float(img)
        gimage = segmentation.inverse_gaussian_gradient(img)
        
        edges = morphology.black_tophat(gimage, footprint=np.ones((3, 3)))

        print("Edges detected")
        
        threshold_value = threshold_otsu(edges)
        binary_image = edges > threshold_value
        
        labels = measure.label(binary_image)
        labels = morphology.closing(segmentation.clear_border(labels), footprint)
        
        labels = measure.label(morphology.remove_small_holes(labels, 1200))
        
        # remove labels smaller than 5000 (off target labels)
        labels_areas = np.unique(labels, return_counts=True)[1]
        area_label_remove = labels_areas < 5000
        
        labels_to_remove = np.unique(labels, return_counts=True)[0][area_label_remove]

        mask = np.in1d( labels, labels_to_remove ).reshape( labels.shape )
        labels[mask] = 0 
        
        labels_array[i,:,:] = labels > 0


# imgs_dict['img_array'] = img_array
imgs_dict['labels_array'] = labels_array

save_object(imgs_dict, f'../../iPSC_imaging/masks/{folder_date_time}_img_labels')
