import numpy as np
import statsmodels.api as sm
import mooda as md
import pandas as pd
import re
import matplotlib.pyplot as plt
from scipy import stats, interpolate
from io import StringIO
from datetime import datetime

class Kdunio_config:

    def __init__(self, instrument, start_date, start_time, end_date, end_time):
        self.instrument = instrument
        self.start_date = start_date
        self.start_time = start_time
        self.end_date = end_date
        self.end_time = end_time

    def duration(self):
        return datetime(self.end_time) - datetime(self.start_time)

def get_depths(df):
    """
    Find all depths that are included in column names
    and return a list of this depths
    Returns
    -------
        depth_list: list
            list of depths extracted from column names
    """
    # get list of columns values
    pre_depth_list = str(df.columns)
    # extract numbers of this values
    depth_list = re.findall(r'\d+\.\d+', pre_depth_list)
    # convert numbers to float
    depth_list = [float(item) for item in depth_list]
    # convert list to a Series
    # depth_list = pd.Series(depth_list)

    return depth_list

def slice_df(start_date_time, end_date_time, df_clear):
    """
    Function to slice time of dataframe in the specific range using mask
    Parameters
    ----------
        start_date_time: datetime
            Start datetime of our range
        end_date_time: datetime
            End datetime of our range
        df_clear: Pandas Dataframe
            Dataframe filtered by full spectral columns
    Returns
    -------
        df_clear: Pandas Dataframe
            Dataframe filtered by full spectral columns and time
    """
    mask = (df_clear.index >= start_date_time) & (df_clear.index <= end_date_time)
    df_clear = df_clear.loc[mask]

    return df_clear

# calculate mean standard deviation and log
def calc_df(df_clear):
    """
    Function to calculate mean, standard deviation and log of counts
    Parameters
    ----------
        df_clear: Pandas Dataframe
            Dataframe filtered by full spectral columns and time
    Returns
    -------
        df_clear: Pandas Dataframe
            Dataframe with new columns of "mean", "stdev" (standard deviation), and "log"
    """
    df_clear['mean'] = df_clear.mean(axis=1)
    df_clear['stdev'] = df_clear.std(axis=1)
    log_row = np.log(df_clear['mean'])
    df_clear['log_row'] = log_row

    return df_clear

def Kd_anomaly_checker(depths, log_irradiance, threshold=0.1, subset_size=7):
    """
    Function to identify anomalous shifts in the data being used for Kd estimate
    Parameters
    ----------
        depths: An array of depths
        log_irradiance: an array of logged irradiance values
        threshold: Std_err threshold at which to raise flag in first test
        subset_size: number of points to use for subsample in second test
    Returns
    -------
        anomaly_flags: A pair of True/False flags for the 2 anomaly tests
    """
    # using all the data test against the threshold
    x = depths
    y = log_irradiance

    # get coeffs of linear fit
    results = sm.OLS(y1, x1).fit()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    if std_err>threshold:
        print('anomaly_flag1 raised')
        anomaly_flag_1=True
    else:
        anomaly_flag_1=False
    if (subset_size % 2) == 0 or subset_size<3:
        print("subset size must be an odd integer >= 3")
        raise ValueError
    last_point=len(x)-1
    #Now see if the relationship is consistent with different subsets
    x_1, y_1=x[0:subset_size], y[0:subset_size]
    x_1=sm.add_constant(x_1)
    res1=sm.OLS(y_1, x_1).fit()
    #print(res1.summary())
    x_2, y_2=x[int(last_point/2-np.floor(subset_size/2)):int(last_point/2+np.ceil(subset_size/2))], y[int(last_point/2-np.floor(subset_size/2)):int(last_point/2+np.ceil(subset_size/2))]
    x_2=sm.add_constant(x_2)
    res2=sm.OLS(y_2, x_2).fit()
    #print(res2.summary())
    x_3, y_3=x[last_point-subset_size:], y[last_point-subset_size:]
    x_3=sm.add_constant(x_3)
    res3=sm.OLS(y_3, x_3).fit()
    #print(res3.summary())
    all_res=[res1, res2, res3]
    smin=np.argmin([res1.params['depths'],res2.params['depths'],res3.params['depths']])
    smax=np.argmax([res1.params['depths'],res2.params['depths'],res3.params['depths']])

    #For z stat use see Paternoster, Raymond, et al. "Using the correct statistical test for the equality of regression coefficients." Criminology 36.4 (1998): 859-866.
    z=(all_res[smax].params['depths']-all_res[smin].params['depths'])/(np.sqrt(all_res[smax].bse['depths']**2+all_res[smin].bse['depths']**2))
    #print(z)
    # Significant at 95% on a 1-tail test?
    if z> 1.6448536269514722:
        print('anomaly_flag2 raised')
        anomaly_flag_2=True
    else:
        anomaly_flag_2=False
    print('')
    anomaly_flags=[anomaly_flag_1, anomaly_flag_2]
    return anomaly_flags

#sub functions required for the analysis_kduino function
def merge_metadata(dict1, dict2):
    # Merge dictionaries
    dict3 = {**dict1, **dict2}
    # Iterate over items in new dictionary
    for key, value in dict3.items():
        # If keys are in both dictionaries
        if key in dict1 and key in dict2:
            # If dictionary contains list of elements
            if isinstance(value, list):
                # If values of new dict and values from parameter dict are different,
                # and not included in the new dict
                if (dict1[key] not in value) and (set(dict1[key]) != set(value)):
                    dict3[key].append(dict1[key])
                elif (dict2[key] not in value) and (set(dict2[key]) != set(value)):
                    dict3[key].append(dict2[key])
            # If dictionary not contains list of elements
            else:
                if value != dict1[key]:
                    dict3[key] = [value, dict1[key]]
                elif value != dict2[key]:
                    dict3[key] = [value, dict2[key]]
    return dict3

def create_wf(metadata_list, data_list, index):
    df = data_list[index]
    metadata = metadata_list[index]
    df.columns = range(df.shape[1])
    # Delete unused columns
    if len(df.columns) > 4:
        df_copy = df.copy()
        ncol = len(df.columns)
        x = range(4, ncol)
        df = df.drop(x, axis=1)
    df.columns = range(df.shape[1])
    # Creation of WaterFrame
    wf = md.WaterFrame()
    # Copy metadata to waterframe
    wf.metadata = metadata
    depth = wf.metadata["depth"]
    # Set name of parameters
    param_red = f'RED_{depth}'
    param_green = f'GREEN_{depth}'
    param_blue = f'BLUE_{depth}'
    param_clear = f'CLEAR_{depth}'
    # Set name of QC parameters
    param_red_qc = f'RED_{depth}_QC'
    param_green_qc = f'GREEN_{depth}_QC'
    param_blue_qc = f'BLUE_{depth}_QC'
    param_clear_qc = f'CLEAR_{depth}_QC'
    # Init data of waterframe
    wf.data[param_red] = df[0]
    wf.data[param_green] = df[1]
    wf.data[param_blue] = df[2]
    wf.data[param_clear] = df[3]
    # Create vocabulary
    for param in [param_red, param_green, param_blue, param_clear]:
        wf.vocabulary[param] = {'units': "counts"}
    # Resample to seconds
    wf.resample('S')

    # Delete last index because it is a minute that we are not going to use
    wf.data.drop(wf.data.tail(1).index, inplace=True)

    # Extract data of the dataframe df. Put all counts in the proper column
    red_list = []
    green_list = []
    blue_list = []
    clear_list = []
    for j in range(len(df_copy.index) - 1):
        for i in range(len(df_copy.columns)):
            if i % 4 == 0:
                red_list.append(df_copy[i][j])
                green_list.append(df_copy[i + 1].iloc[j])
                blue_list.append(df_copy[i + 2].iloc[j])
                clear_list.append(df_copy[i + 3].iloc[j])

    red_array = np.array(red_list)
    green_array = np.array(green_list)
    blue_array = np.array(blue_list)
    clear_array = np.array(clear_list)

    wf.data[param_red] = red_array
    wf.data[param_green] = green_array
    wf.data[param_blue] = blue_array
    wf.data[param_clear] = clear_array

    # Init waterframe QC data
    wf.data[param_red_qc] = 0
    wf.data[param_green_qc] = 0
    wf.data[param_blue_qc] = 0
    wf.data[param_clear_qc] = 0

    return wf

def make_depth_plots(depths_row_clear, row_clear, slope, _intercept, r_value, index, save=False):
    # plot depths - values clear
    plt.plot(depths_row_clear, row_clear, marker='o', linestyle="")
    plt.plot(depths_row_clear, slope * depths_row_clear + _intercept)

    plt.xlabel("Depth (m)")
    plt.ylabel("ln(Light)")
    plt.legend([f"Kd = {np.around(slope * (-1), 3)}",
                f'r\N{SUPERSCRIPT TWO} = {np.around(r_value ** 2, 3)}'],
               loc='best', markerscale=0)
    plt.title(f"Kd at time: {index}")
    if save==True:
        plt.savefig(f'./Output_depth_plot_{index}.png')
    else:
        plt.show()

def plot_Kd_timeseries(wf_clear, save=False):
    # Plot Kd
    print("\n\nPlot time series Kd PAR")
    fig, ax = plt.subplots()
    ax.set_ylabel('$K_d$ ($m^{-1}$)')
    ax.xaxis.set_tick_params(rotation=45)
    ax.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.plot(wf_clear.data['Kd'], marker='o', linestyle="-")
    plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=10))
    plt.xlabel("Time (minutes)")
    # plt.ylabel("Kd (m-1)")
    ax.legend(["Kd"], bbox_to_anchor=(1.25, 0.95), borderaxespad=1., ncol=1, fontsize=8, title_fontsize=10,
              prop={'size': 10})
    plt.title("Timeseries $K_d$")
    # plt.rcParams['figure.figsize'] = [20/2.54, 16/2.54]
    if save == True:
        plt.savefig(f'./Output_Kd_ts_plot_{wf_clear.data.index}.png')
    else:
        plt.show()

def plot_Kd_timeseries_withr2(wf_clear, save=False):
    # Plot Kd
    print("\n\nPlot time series Kd PAR with coefficient of determination r2")
    plt.close('all')
    fig, ax = plt.subplots()
    twin1 = ax.twinx()
    twin1.set_ylabel('$r^2$')
    ax.set_ylabel('$K_d$ ($m^{-1}$)')
    ax.xaxis.set_tick_params(rotation=45)
    ax.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    p1, = ax.plot(wf_clear.data['Kd'], color="black", linestyle='-',
                  label="$K_d$")
    p2, = twin1.plot(wf_clear.data['r2'], color="peru", linestyle='-', label="$r^2$")

    plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=10))
    ax.legend(handles=[p1, p2], bbox_to_anchor=(1.35, 0.95), borderaxespad=1.,
              ncol=1,fontsize=10, title_fontsize=10, prop={'size': 10})
    ax.set_xlabel("Time (minutes)")
    plt.title("Timeseries $K_d$")
    # plt.rcParams['figure.figsize'] = [20/2.54, 16/2.54]
    if save == True:
        plt.savefig(f'./Output_Kd_ts_r2_plot_{wf_clear.data.index}.png')
    else:
        plt.show()

def analysis_kduino(instrument, datafile_list):
    # Function to analyse data from KdUINO to obtain Kd
    if instrument != 'KduPRO':
        print("KdUINO instrument not configured well or unknown")
    else:
        print("KduPRO")
        # print(datetime_start)
        # print(datetime_stop)
        # Definitions for regular expression paterns
        start_string_metadata = r"METADATA"
        stop_string_metadata = r"DATA"
        start_string_data = r"\bDATA"
        stop_string_data = r"METADATA"
        last_start_string_data = r"\bDATA"
        end_string_data = r'$(?![\r\n])'
        metadata_patron = re.compile(r'{}(?P<length>)\s*(?P<table>[\s\S]*?){}'.format(
            start_string_metadata, stop_string_metadata))
        data_patron = re.compile(r'{}(?P<length>)\s*(?P<table>[\s\S]*?){}'.format(
            start_string_data, stop_string_data))
        last_data_patron = re.compile(r'{}(?P<length>)\s*(?P<table>[\s\S]*?){}'.format(
            last_start_string_data, end_string_data))
        selected_info = ""
        metadata_list = []
        metadata = {}
        data_list = []
        kdupro_list = []
        for k, v in datafile_list():
            with open(k) as reader:
                content = reader.read()

                # Regular expression to find the metadata patron
                for m in re.finditer(metadata_patron, content):
                    selected_info = m.group('table')
                    metadata = {}
                    lines = selected_info.splitlines()

                    for line in lines:
                        key = line.split(":")[0]
                        if line.count(":") > 1:
                            date_splitted = (line.rsplit(":")[-3:])
                            date_splitted = " ".join(date_splitted)
                            value = date_splitted
                            metadata[key] = value
                        else:
                            value = line.split(":")[1]
                            metadata[key] = value.strip()
                        metadata_list.append(metadata)

                # Regular expression to find the data patron
                for d in re.finditer(data_patron, content):
                    selected_info_data = d.group('table')
                    data = StringIO(selected_info_data)
                    df = pd.read_csv(data, skipinitialspace=True, skiprows=1, header=None,
                                     parse_dates={'TIME': [0]}, delimiter=' ',
                                     engine='python').set_index('TIME')
                    data_list.append(df)

                for m in re.finditer(last_data_patron, content):
                    selected_info_data = m.group('table')
                    index_last_data = selected_info_data.rfind('DATA')
                    selected_info_last_data = selected_info_data[index_last_data:]
                    data = StringIO(selected_info_last_data)
                    df = pd.read_csv(data, skipinitialspace=True, skiprows=2, header=None,
                                     parse_dates={'TIME': [0, 1]}, delimiter=' ',
                                     engine='python').set_index('TIME')
                    data_list.append(df)

                for index, df in enumerate(data_list):
                    if datetime_start in df.index or datetime_stop in df.index:
                        kdupro_list.append(index)

        # process the files we have selected
        wf_list = []

        for index in kdupro_list:
            wf = create_wf()
            wf_list.append(wf)

            # Declare lists
            names = []
            depths = []

            # Create unique waterframe
            wf_all = md.WaterFrame()

            # Concat all waterframes
            for index, wf in enumerate(wf_list):
                if index == 0:
                    wf_all = wf.copy()
                else:
                    # Concat data
                    wf_all.data = pd.concat([wf_all.data, wf.data], axis=1)
                    # Add metadata
                    wf_all.metadata = merge_metadata(wf.metadata, wf_all.metadata)
                    # Add vocabulary
                    for param in wf.parameters:
                        wf_all.vocabulary[param] = wf.vocabulary[param]

            # Append names and depths to each list
            for wf in wf_list:
                if wf is None:
                    pass
                else:
                    name = wf.metadata["name"]
                    names.append(name)
                    depth = wf.metadata["depth"]
                    depths.append(depth)

                # Slice time
                mask = (wf_all.data.index >= datetime_start) & (
                        wf_all.data.index <= datetime_stop)
                wf_all.data = wf_all.data.loc[mask]

                # Resample time
                wf_all.resample("T")

                # Convert depths list elements to float, and save it as a numpy array
                depths = np.array(list(map(float, depths)))

                # Save CLEAR data in one new waterframe
                wf_clear = md.WaterFrame()
                wf_clear.metadata = wf_all.metadata
                wf_clear.vocabulary = wf_all.vocabulary

                match_CLEAR = [s for s in wf_all.data if ("CLEAR" in s) and ("QC" not in s)]
                wf_clear.data = wf_all.data.filter(match_CLEAR)

                # Calculate Kd in new column
                kd_list = []
                r2_list = []
                wf_clear.data['Kd'] = 0
                wf_clear.data['r2'] = 0

                for index, _row in wf_clear.data.iterrows():
                    # CLEAR
                    row_clear = wf_clear.data.loc[index, match_CLEAR].tolist()
                    # calculate Kd from linear regression
                    slope, _intercept, r_value, _p_value, _std_err = stats.linregress(
                        depths, np.log(row_clear))
                    kd_list.append(slope * (-1))
                    r2_list.append(r_value ** 2)

                    depths_row_clear = np.array(depths)
                    row_clear = np.array(np.log(row_clear))
                    make_depth_plots(depths_row_clear, row_clear, slope, _intercept, r_value, index)
                    plt.close('all')

                    wf_clear.data['Kd'] = kd_list
                    wf_clear.data['r2'] = r2_list
                    print(wf_clear.data)
                    print("\n")
                    avgKd = wf_clear.data.loc[wf_clear.data['r2'] > 0.9, 'Kd'].mean()
                    print("Average Kd for samples with r2 > 0.9: {}".format(np.around(avgKd, 3)))

                plot_Kd_timeseries(wf_clear)
                plot_Kd_timeseries_withr2(wf_clear)



