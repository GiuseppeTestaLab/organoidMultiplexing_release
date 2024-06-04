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

def preprocess(summary_df, original_v = 'perc_area', final_output = 'mean', scale = False, **kwargs):

    import numpy as np
    import pandas as pd
    import datetime
    
    stds = summary_df[['time_point', original_v]].groupby('time_point').std()

    if final_output == 'mean':
        mean_df = summary_df[['time_point', original_v]].groupby('time_point').mean()
    
    if final_output == 'area_sum':
        mean_df = summary_df[['time_point', original_v]].groupby('time_point').sum()

    stds_dict = {i:j for i, j in zip(stds.index, stds[original_v])}
    mean_dict = {i:j for i, j in zip(mean_df.index, mean_df[original_v])}

    summary_df = pd.DataFrame.from_dict(mean_dict, columns=[final_output], orient='index')
    summary_df['stds'] = summary_df.index.map(stds_dict)
    summary_df = summary_df.reset_index()
    summary_df.columns = ['time_point', final_output, 'stds']

    summary_df['hour'] = summary_df['time_point'].apply(lambda x: int(x.split('_')[-1].strip('t')))
    summary_df['month'] = summary_df['time_point'].apply(lambda x: x.split('_')[1]) 
    summary_df['day'] = summary_df['time_point'].apply(lambda x: x.split('_')[0])
    summary_df['datetime'] = [datetime.datetime(year = 2023, month = m, day = d, hour=h) for m, d, h in zip(summary_df['month'].astype(int), summary_df['day'].astype(int), summary_df['hour'].astype(int))]
    summary_df = summary_df.sort_values(by = 'datetime')

    summary_df['split_time'] = summary_df['datetime'].transform(lambda x: (x - x.min()).dt.total_seconds())
    summary_df['split_time'] = summary_df['split_time'] / 3600

    if scale:

        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler(**kwargs)

        ydata = np.array(summary_df[final_output].values).reshape(-1, 1)
        scaler.fit(ydata)
        ydata = scaler.transform(ydata)
        summary_df[final_output] = ydata

    return summary_df

def highlight_growth_curves(all_lines, 
                            xlabel = 'Hours from split', 
                            ylabel = 'Total area normalized on t0', 
                            lines = None,
                           fontsize = 30):

    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    #plt.style.use('seaborn-v0_8-whitegrid')

    if lines is None:
        lines = all_lines.keys()
    
    nrows = int(len(lines) / 4)

    if len(lines) % 4 > 0:
        nrows += 1\

    ncols = 4

    split_palette = {'1': '#264653', '2': '#2a9d8f', '3': '#8ab17d', '4': '#e9c46a', '5': '#f4a261', '6': '#e76f51'}

    fig, ax = plt.subplots(nrows, ncols, figsize = (ncols * 8, nrows * 7), 
                           gridspec_kw= {'hspace': 0.6, 'wspace': .5})
    ax = ax.flatten().T

    for l, x in zip(lines, ax):
        
        all_keys_to_color = [i for i in all_lines.keys() if l in i]
        #print(all_keys_to_color)
        all_keys_not_to_color = [i for i in all_lines.keys() if l not in i]
        #print(all_keys_not_to_color)

        for k in all_keys_to_color:
            split_n = k.split('_')[-1]
            line = all_lines[k]
            new_line = mlines.Line2D(line.get_xdata(), line.get_ydata(), color=split_palette[split_n])
            new_line.set_label(k)
            x.add_line(new_line)

        for k in all_keys_not_to_color:
            line = all_lines[k]
        # Create a new Line2D object with the same properties
            new_line = mlines.Line2D(line.get_xdata(), line.get_ydata(), color='lightgray', alpha = .3)
            x.add_line(new_line)
            
        x.autoscale()
        x.set_xlabel(xlabel, fontdict={'size': fontsize})
        x.set_ylabel(ylabel, fontdict={'size': fontsize})
        x.set_xticklabels(x.get_xticklabels(), fontdict={'size': fontsize}, rotation = 90)
        x.set_yticklabels(x.get_yticklabels(), fontdict={'size': fontsize})
        x.set_title(l, fontdict={'size': fontsize + 10})
        x.legend(fontsize= fontsize - 10)

    fig.show()