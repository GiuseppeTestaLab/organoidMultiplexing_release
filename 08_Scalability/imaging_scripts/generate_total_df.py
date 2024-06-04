# %%
from skimage import measure
import pandas as pd
import numpy as np
import os
import datetime

def load_object(filename):

    import pickle

    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported): ", ex)

def save_object(obj, filename):

    import pickle

    try:
        with open(filename + ".pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)

# %%
total_df = pd.DataFrame()

folder = '../../iPSC_imaging/'

for file in os.listdir(folder + 'masks/'):
    

    if file.endswith('.pickle'):
        file_name = file.strip('.pickle')
        
        print(file)

        imgs_dict = load_object(f'{folder}masks/{file}')
        
        if 'img_array' in imgs_dict:
            del imgs_dict['img_array']
            save_object(imgs_dict, f'{folder}masks/{file_name}')

        transposed = {i: {} for i in imgs_dict['file_names']}

        for label, i in zip(imgs_dict["labels_array"], transposed):
            
            #transposed[i]['img_array'] = img
            #transposed[i]['labels_array'] = label
            transposed[i]['total_area'] = np.sum(label)
            transposed[i]['perc_area'] = transposed[i]['total_area'] / 9519251 * 100
            labeled_area = measure.label(label)
            label_and_area = np.unique(labeled_area, return_counts=True)
            if len(label_and_area[0]) > 1:
                colonies_areas = label_and_area[1][1:]
                transposed[i]['mean_area_per_colony'] = colonies_areas.mean()
                transposed[i]['n_colonies'] = len(colonies_areas)
                del labeled_area

            else:
                colonies_areas = [0]
                transposed[i]['mean_area_per_colony'] = 0
                transposed[i]['n_colonies'] = 0
                del labeled_area

        df = pd.DataFrame.from_dict(transposed).T
        df['time_point'] = file.strip('.pickle')
        total_df = pd.concat([total_df, df])
      
        del imgs_dict
        del transposed
        del colonies_areas
        del df

total_df.to_csv(folder + 'quantification_raw.csv')

# %%
total_df = pd.read_csv(folder + 'quantification_raw.csv', index_col=0)
total_df

# %%
total_df['time_point'] = total_df['time_point'].str.strip('_img_labels')
total_df['confluency/generation'] = total_df['time_point'].apply(lambda x: 'confluency' if 'confluency' in x else ('generation'))
total_df['time_point'] = total_df['time_point'].apply(lambda x: x.split('+')[0] if 'confluency' in x or 'generation' in x else (x))

total_df['hour'] = total_df['time_point'].apply(lambda x: int(x.split('_')[-1].strip('t')))
total_df['month'] = total_df['time_point'].apply(lambda x: x.split('_')[1]) 
total_df['day'] = total_df['time_point'].apply(lambda x: x.split('_')[0])
total_df['line'] = total_df.reset_index()['index'].apply(lambda x: x.split('_')[0]).values
total_df.line = total_df.line.replace('CTL04', 'CTL04E')
total_df['datetime'] = [datetime.datetime(year = 2023, month = m, day = d, hour=h) for m, d, h in zip(total_df['month'].astype(int), total_df['day'].astype(int), total_df['hour'].astype(int))]

# %%
mean_df_time_point = total_df.groupby(['time_point']).mean('perc_area')
mean_df_time_point_dict = {i:j for i, j in zip(mean_df_time_point.index, mean_df_time_point.perc_area)}
mean_df_time_point_dict

# %%
total_df['norm_factor'] = total_df.time_point.map(mean_df_time_point_dict)
total_df['perc_area_norm'] = total_df['perc_area'] / total_df['norm_factor']
total_df['err_bar_mean'] = total_df['perc_area'] - total_df['perc_area_norm']
total_df['n_split'] = pd.Series(total_df.index).apply(lambda x: x.split('_')[1]).values

total_df['split_time'] = total_df.groupby(['line', 'n_split'])['datetime'].transform(lambda x: (x - x.min()).dt.total_seconds())
total_df['split_time'] = total_df['split_time'] / 3600 # convert seconds to hours

# %%
total_df.n_split = total_df.n_split.replace('day', '1')


total_df.to_csv(folder + 'quantification.csv')
