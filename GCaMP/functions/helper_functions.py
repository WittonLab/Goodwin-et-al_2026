import numpy as np
import pandas as pd

def bin_data_trace(coordinates, traces, bin_edges):
    """
    Bins data traces based on corresponding coordinate values and specified bin edges.
    
    Parameters:
    coordinates (np.ndarray): 1D array of coordinates used to bin the data.
    traces (np.ndarray): 1D or 2D array where each column is a data trace to be binned.
                         Shape is (num_data_points, num_traces).
    bin_edges (np.ndarray): 1D array specifying the edges of the bins.
    
    Returns:
    np.ndarray: 2D array where each row represents the binned data for a corresponding trace.
                Shape is (num_bins, num_traces).
    """

    # Check if traces is a 1D array and reshape to 2D if needed
    if traces.ndim == 1:
        traces = traces.reshape(-1, 1)
    
    binned_traces = []
    # Iterate over each trace to bin the data
    for i in range(traces.shape[1]):
        trace = traces[:, i]
        
        # Bin the coordinates based on the provided bin edges
        coord_bins = np.digitize(coordinates, bins=bin_edges, right=True)
        
        # Initialize an array to store the binned data
        binned_data = np.zeros(len(bin_edges) - 1)
        
        # Calculate the mean value of the trace data for each bin
        for bin_num in range(1, len(bin_edges)):
            bin_indices = np.where(coord_bins == bin_num)[0]  # Indices of data points in the current bin
            
            if len(bin_indices) > 0:
                bin_mean = np.mean(trace[bin_indices])  # Calculate the mean for the bin
            else:
                bin_mean = 0  # Set to zero if no data falls into the bin
                
            binned_data[bin_num - 1] = bin_mean  # Store the mean in the binned data array
        
        binned_traces.append(binned_data)  # Append the binned data for the current trace

    # Convert the list of binned traces to a 2D array and transpose it
    return np.array(binned_traces).T


def boxcar_filter(data, windowsize):
    """
    Applies a boxcar (moving average) filter to a list of 1D data arrays.

    Parameters:
    - data (list of np.ndarray): A list where each element is a 1D numpy array representing a signal.
    - windowsize (int): The size of the moving average window.

    Returns:
    - np.ndarray: A numpy array where each element is the filtered version of the corresponding input data array.
    """
    
    # Create a boxcar (moving average) window
    window = np.ones(windowsize)
    boxcar = window / windowsize
    
    output = []
    # Apply the boxcar filter to each 1D array in the data list
    for t in data:
        # Convolve the current data array with the boxcar window
        filtered = np.convolve(t, boxcar, mode='same')
        output.append(filtered)
    
    # Return the filtered data as a numpy array
    return np.array(output)


def calculate_stdev_snr(transient_traces, raw_traces):
    """
    Calculates the Signal-to-Noise Ratio (SNR) for each transient trace using the 
    standard deviation of the baseline noise.

    Parameters:
    transient_traces (np.ndarray): 2D array of shape (num_traces, num_data_points) 
                                   where each row is a significant transients trace.
    raw_traces (np.ndarray): 2D array of shape (num_traces, num_data_points) 
                             where each row is the corresponding raw trace used for 
                             noise calculation.

    Returns:
    list: A list of SNR values, one for each transient trace.
    """

    snr = []
    # Iterate over each transient trace
    for i in range(transient_traces.shape[0]):
        
        # Extract the current transient trace and raw trace
        transient = transient_traces[i,:]
        trace = raw_traces[i,:]
        # Determine baseline noise by selecting portions of the trace where the transient_trace is zero
        baseline = trace[transient == 0]
        # Calculate the noise as the standard deviation of the baseline
        noise = np.std(baseline)
        # Determine the signal as the maximum value of the transient
        signal = np.max(transient)
        # Calculate the SNR and append to the list
        snr.append(signal/noise)

    return snr


def calculate_mad_snr(transient_traces, raw_traces):
    """
    Calculates the Signal-to-Noise Ratio (SNR) for each transient trace using the 
    Median Absolute Deviation (MAD) as the noise estimate.

    Parameters:
    transient_traces (np.ndarray): 2D array of shape (num_traces, num_data_points) 
                                   where each row is a significant transients trace.
    raw_traces (np.ndarray): 2D array of shape (num_traces, num_data_points) 
                             where each row is the corresponding raw trace used for 
                             noise calculation.

    Returns:
    list: A list of SNR values, one for each transient trace.
    """
    
    # calculated the signal to noise ratio based on the MAD median absolute deviation
    snr = []
    for i in range(transient_traces.shape[0]): 
        
        transient = transient_traces[i,:]
        trace = raw_traces[i,:]
        noise = np.median(np.abs(trace- np.median(trace)))
        #not sure which is better
        signal = np.max(transient)
        #signal = np.median(transient[transient>0])
        snr.append(signal/noise)

    return snr


def baseline_correction(trace, sl_window, percentile_value):
    """
    Perform baseline correction on the input trace.

    Args:
    - trace (numpy.ndarray): Input trace data.
    - sl_window (int): sliding window length.
    - percentile_value (float): lower percentile to define baseline.

    Returns:
    - corrected_trace (numpy.ndarray): Trace data after baseline correction.
    """
    # Pre-calculate indices for sliding window
    num_frames = trace.shape[0]
    start_indices = np.maximum(0, np.arange(num_frames) - sl_window)
    end_indices = np.minimum(num_frames, np.arange(num_frames) + sl_window)

    # Initialize an array to store window values
    window_values = np.zeros(num_frames)

    # Iterate over each frame
    for j in range(num_frames):
        # Determine the start and end indices of the current frame
        start_idx = start_indices[j]
        end_idx = end_indices[j]

        # Calculate the window value as the percentile of trace values within the frame
        window_values[j] = np.percentile(
            trace[start_idx:end_idx], 100 * percentile_value
        )

    # Perform baseline correction by subtracting window values and normalizing
    #corrected_trace = (trace - window_values) / window_values
    corrected_trace = trace  / window_values

    median = np.median(corrected_trace)
    corrected_trace = (corrected_trace - median) / median # This is now deltaF/F0 

    return corrected_trace

# --- Function for ΔF/F ---
def deltaF_over_F(trace):
    """
    Compute ΔF/F using the median of the whole trace as F0.
    
    Args:
    - trace (numpy.ndarray): Input fluorescence trace.
    
    Returns:
    - dFF (numpy.ndarray): ΔF/F trace
    """
    F0 = np.median(trace)
    if F0 == 0:
        F0 = np.finfo(float).eps  # prevent divide by zero
    dFF = (trace - F0) / F0
    return dFF

def get_significant_transients(Fc):
    """
    Identify and extract significant transients from calcium fluorescence traces.

    This function processes calcium fluorescence data to detect significant transients
    by identifying positive and negative transients, and filtering them based on their 
    significance level. The output is a modified trace where only significant transients
    are preserved and baseline set to 0.

    Args:
        traces (numpy.ndarray): 2D array of calcium fluorescence data, where each row
                                represents a trace over time. (number_frames, number_traces).

    Returns:
        numpy.ndarray: A 2D array with the same shape as `traces`, containing only the
                       significant transients, with the values expressed as ΔF/F0.

    Key Processing Steps:
        1. **Transient Detection:** Positive and negative transients are identified using
           double thresholding based on a predefined start and end threshold.
        2. **Transient Filtering:** The detected transients are filtered based on their
           significance level by comparing their duration and amplitude against noise thresholds.
        3. **Reconstruction:** The function reconstructs the traces, preserving only the
           significant transients, and converts the results back into ΔF/F0 format.
    """

    
    # combines lists of different lens into an np.array padding with zeros
        # combines lists of different lens into an np.array padding with zeros
    def lists2array(list_of_lists):
        """
        Converts a list of lists into a 2D numpy array.

        Args:
        - list_of_lists (list): List containing lists of variable lengths.

        Returns:
        - result_array (numpy.ndarray): 2D numpy array with padded rows to ensure equal length.
        """
        # Find the maximum length of sublists
        max_length = max(len(sublist) for sublist in list_of_lists)
        
        # Initialize the result array with zeros
        result_array = np.zeros((len(list_of_lists), max_length))

        # Fill the result array with values from sublists
        for i, sublist in enumerate(list_of_lists):
            result_array[i, :len(sublist)] = sublist
            
        return result_array


    def double_thresh(data, thresh1, thresh2):
        """
        Perform double thresholding on the input data array.

        Args:
        - data (numpy.ndarray): Input data array.
        - thresh1 (float): Upper threshold.
        - thresh2 (float): Lower threshold.

        Returns:
        - state (numpy.ndarray): Array indicating the state based on thresholds.
        """
        # Initialize the state array with zeros
        state = np.zeros_like(data)
        
        # Apply thresholding conditions
        state[data > thresh1] = 1
        state[data <= thresh2] = 0
        
        # Apply double thresholding logic
        for i in range(1, len(data)):
            if state[i - 1] == 0 and data[i] > thresh1:
                state[i] = 1
            elif state[i - 1] == 1 and data[i] < thresh2:
                state[i] = 0
        
        return state


    def get_transients(Fc2_std, str_threshold = 2.0, end_threshold = 0.5):
        """
        Extracts transients from the standard deviation calcium fluorescence data.

        Args:
        - Fc2_std (numpy.ndarray): Standard deviation calcium fluorescence data.

        Returns:
        - transient_max (numpy.ndarray): Array containing the maximum value (amplitude) of each transient in std.
        - transient_duration (numpy.ndarray): Array containing the duration of each transient.
        - transient_start (numpy.ndarray): Array containing the start index of each transient.
        """
        # Initialize transient state
        transient_state = np.zeros_like(Fc2_std)

        # Find transient state using double thresholding
        transient_state = double_thresh(Fc2_std, str_threshold, end_threshold)
        
        # Find start and end indices of transients
        transient_start = np.where(np.diff(transient_state) == 1)[0]
        transient_end = np.where(np.diff(transient_state) == -1)[0]

        # Process transients if any found
        if transient_start.size > 0:
            # Handle edge cases
            if transient_start[0] > transient_end[0]:
                transient_start = np.insert(transient_start, 0, 0)
            if transient_start[-1] > transient_end[-1]:
                transient_end = np.append(transient_end, len(Fc2_std) - 1)
            
            # Calculate maximum value and duration of each transient
            transient_max = [np.max(Fc2_std[transient_start[j]:transient_end[j] + 1]) for j in range(len(transient_start))]
            transient_duration = [transient_end[j] - transient_start[j] for j in range(len(transient_start))]
            
            return transient_max, transient_duration, transient_start
        else:
            # Return zeros if no transients found
            return np.zeros(0), np.zeros(0), np.zeros(0)
    


    # Upper limit for noise
    noise_stdev_mult = 2

    # Transient start and end thresholds
    trans_str_thresh = 2.0
    trans_end_thresh = 0.5

    # Transpose Fc
    Fc_trans = Fc.T
    Fc2 = Fc_trans.copy()

    # Number of traces
    num_traces = Fc2.shape[1]

    num_frames = Fc2.shape[0]

    # Arrays to store standard deviations and transient states
    Fc2_std = np.zeros_like(Fc2)
    noise_thresh = np.zeros(num_traces)

    # Lists to store transient information
    pos_trans_max_list = []
    pos_trans_dur_list = []
    neg_trans_max_list = []
    neg_trans_dur_list = []
    pos_trans_start_list = []
    neg_trans_start_list = []    


    # Define baseline and extract transient information
    for i in range(num_traces):
        # Calculate noise threshold for the baseline
        noise_thresh[i] = noise_stdev_mult * np.std(Fc2[:, i])
        base_idx = np.where(Fc2[:, i] < noise_thresh[i])[0]
        base = Fc2[base_idx, i]
        base_med = np.median(base)
        base_std = np.std(base)
        
        # Normalize to baseline standard deviation
        Fc2_std[:, i] = (Fc2[:, i] - base_med) / base_std

        # Get positive transients
        pos_trans_max, pos_trans_dur, pos_trans_start = get_transients(Fc2_std[:, i], 
                                                                    trans_str_thresh, 
                                                                    trans_end_thresh)
        pos_trans_max_list.append(pos_trans_max)
        pos_trans_dur_list.append(pos_trans_dur)
        pos_trans_start_list.append(pos_trans_start)

        # Get negative transients
        neg_trans_max, neg_trans_dur, neg_trans_start = get_transients(-Fc2_std[:, i], 
                                                                    trans_str_thresh, 
                                                                    trans_end_thresh)
        neg_trans_max_list.append(neg_trans_max)
        neg_trans_dur_list.append(neg_trans_dur)
        neg_trans_start_list.append(neg_trans_start)


    # Convert lists to arrays
    neg_trans_max_ar = lists2array(neg_trans_max_list)
    neg_trans_dur_ar = lists2array(neg_trans_dur_list)
    neg_trans_start_ar = lists2array(neg_trans_start_list).astype(int)

    pos_trans_max_ar = lists2array(pos_trans_max_list)
    pos_trans_dur_ar = lists2array(pos_trans_dur_list)
    pos_trans_start_ar = lists2array(pos_trans_start_list).astype(int)

    ## ========================================================================================================

    ## ===================================
    ## Negative/Positive transient ratio
    ## ===================================

    # Different transient durations 
    X = np.arange(1, 26)

    # Initialize histogram arrays
    num_bins = len(X) - 1
    hist_pos_trans = np.zeros((6, num_bins))
    hist_neg_trans = np.zeros((6, num_bins))

    # Iterate over different standard deviations
    for std_idx in range(1, 6):
        # Calculate the current standard deviation
        std = std_idx + 1
        
        if std_idx == 5:
            # For the last standard deviation, consider all transients above the threshold
            pos_trans_mask = (pos_trans_max_ar > std)
            neg_trans_mask = (neg_trans_max_ar > std)
        else:
            # For other standard deviations, consider transients within the threshold range
            pos_trans_mask = np.logical_and(pos_trans_max_ar < (std + 1), pos_trans_max_ar > std)
            neg_trans_mask = np.logical_and(neg_trans_max_ar < (std + 1), neg_trans_max_ar > std)
        
        # Calculate transient durations for positive and negative transients
        pos_trans_dur = pos_trans_dur_ar * pos_trans_mask
        neg_trans_dur = neg_trans_dur_ar * neg_trans_mask
        
        # Calculate histograms of transient durations
        hist_pos_trans[std_idx, :], _ = np.histogram(pos_trans_dur, bins=X)
        hist_neg_trans[std_idx, :], _ = np.histogram(neg_trans_dur, bins=X)

    # Calculate the ratio of negative to positive transient durations (error rate)
    error_rate = hist_neg_trans / hist_pos_trans

    ## ========================================================================================================

    ## ===================================
    ## Find significant transient duration threshold
    ## ===================================

    # Threshold for significance ratio
    significance_thresh = 0.01

    # Initialize array to store minimum durations below threshold for each standard deviation
    min_duration_thresh = np.zeros(6)

    # Iterate over each standard deviation
    for std_idx in range(1, 6):
        # Find indices where the ratio is below the significance threshold
        below_thresh_indices = np.where(error_rate[std_idx, :] < significance_thresh)[0]
        
        if below_thresh_indices.size > 0:
            # If there are durations below threshold, find the minimum duration
            min_below_thresh_value = np.min(X[below_thresh_indices])
            min_duration_thresh[std_idx] = min_below_thresh_value
        else:
            # If no durations are below threshold, set the minimum duration to infinity
            min_duration_thresh[std_idx] = np.inf

    ## ========================================================================================================

    ## ===================================
    ## Filter the traces to get just significant transients
    ## ===================================

    # Initialize array to store results
    Fc3 = np.zeros_like(Fc2)

    # Loop through each trace
    for trace_idx in range(num_traces):
        # Initialize array to store significant transient flags for each standard deviation
        std_sig = np.zeros((num_frames, 6))
        
        # Loop through each standard deviation (only consider 3 std deviations and above)
        for std_idx in range(2, 6):
            std = std_idx + 1  # Calculate the current standard deviation
            
            # Determine significant transients exceeding the minimum duration threshold
            if std_idx == 5:  # For the last standard deviation
                significant_transients = (pos_trans_max_ar[trace_idx] > std) * pos_trans_dur_ar[trace_idx]
                significant_transients = significant_transients * (significant_transients > min_duration_thresh[std_idx])
            else:
                significant_transients = (np.logical_and(pos_trans_max_ar[trace_idx] < (std + 1), pos_trans_max_ar[trace_idx] > std)) * pos_trans_dur_ar[trace_idx]
                significant_transients = significant_transients * (significant_transients > min_duration_thresh[std_idx])
            
            # Loop through each transient start index
            for trans_idx in range(len(pos_trans_start_ar[trace_idx])):
                
                # if there is a transient and this transient is a significant transient
                if pos_trans_start_ar[trace_idx, trans_idx] > 0 and significant_transients[trans_idx] > 0:
                    # Mark the duration of significant transients for each standard deviation
                    start_idx = pos_trans_start_ar[trace_idx, trans_idx] + 1
                    end_idx = pos_trans_start_ar[trace_idx, trans_idx] + int(significant_transients[trans_idx])
                    std_sig[start_idx:end_idx, std_idx] = 1
        
        # Sum up the significant transients for each frame and multiply by the standardized trace
        Fc3[:, trace_idx] = np.sum(std_sig, axis=1) * Fc2_std[:, trace_idx]

    ## ========================================================================================================

    ## ===================================
    ## Convert traces to back into DeltaF / F0
    ## ===================================

    # Convert the traces back to DF/F
    Fc3_DF = np.zeros_like(Fc2)

    for i in range(Fc2.shape[1]):
        baseline_index = np.where(Fc2[:, i] < noise_thresh[i])[0]
        baseline = Fc2[baseline_index, i]
        base_stdev = np.std(baseline)
        Fc3_DF[:, i] = Fc3[:, i] * base_stdev


    return Fc3_DF.T


def get_spatial_information(time_per_bin, activity_per_bin):
    """
    Calculate the spatial information of fluorescence activity across spatial bins.

    This function computes the spatial information, a measure often used in neuroscience
    to quantify the relationship between a neuron's activity and an external variable (e.g., location).
    The formula used is based on a modification of the Shannon information theory concept.

    Parameters
    ----------
    time_per_bin : array-like
        An array representing the amount of time spent in each spatial bin.
        Each entry corresponds to the time spent in that particular bin.
    
    activity_per_bin : array-like
        An array representing the fluorescence activity in each spatial bin.
        Each entry corresponds to the activity recorded in that bin.
    
    Returns
    -------
    spatial_info : float
        The calculated spatial information value. This value reflects how much 
        information about the fluorescence activity is conveyed by knowing the position.
        Higher values indicate a stronger relationship between position and activity.
    
    f_mean : float
        The calculated mean activity rate of the cell across all bins 
        can sometimes be refered to as lamda.
    
    Notes
    -----
    The spatial information is calculated as:
    
        spatial_info = (1 / mean_activity) * Σ(p_i * f_i * log2(f_i / mean_activity))
    
    where:
    - p_i is the probability of being in bin i (time spent in bin i divided by total time).
    - f_i is the activity in bin i.
    - mean_activity is the mean activity across all bins.

    This metric is typically used in spatial encoding studies to assess how well 
    a signal, such as neural activity, encodes information about the position.
    """

    number_bins = len(activity_per_bin)
    T_bins = time_per_bin #time per bin
    T_total = np.sum(time_per_bin) #time total
    f_bins = activity_per_bin #fluor per bin
    f_mean = np.mean(activity_per_bin) #fluor mean
    #print(f_mean)
    spaceinfo_bin = np.zeros_like(activity_per_bin)

    for i in range(number_bins):
        
        fi = f_bins[i]
        p = T_bins[i]/T_total

        spaceinfo_bin[i] = p * fi * np.log2(fi/f_mean)

    spatial_info = np.sum(spaceinfo_bin)/f_mean
    
    return spatial_info, f_mean


def shuffle_transient_trace(trace):
    """
    Randomly shuffles non-zero transient segments within a 1D array without overlap.

    This function takes a 1D numpy array `trace` containing transient segments 
    (continuous non-zero values) and shuffles these segments into new random positions
    within the array. The shuffled segments do not overlap with each other in the 
    resulting array. The rest of the array remains filled with zeros.

    Parameters
    ----------
    trace : numpy.ndarray
        A 1D numpy array containing transient segments to be shuffled. 
        Transients are defined as consecutive non-zero values.

    Returns
    -------
    numpy.ndarray
        A 1D numpy array of the same length as `trace`, containing the original 
        transient segments shuffled into new positions without overlap. Non-transient 
        positions remain as zeros.

    Notes
    -----
    - Transient segments are identified as continuous sequences of non-zero values 
      in the input array.
    - The function ensures that no two transient segments overlap in the output array, 
      by randomly selecting positions and checking for overlaps before placing each segment.
    - The length and content of each transient segment remain unchanged.
    - The random placement of transients may lead to different outputs each time the function is called.
    """

    import time  # At the top of your script if not already imported

    
    # Finds cases where transient starts or ends at first or last datapoint
    if trace[0] > 0:
        trace[0] = 0
    
    if trace[-1] > 0:
        trace[-1] = 0

    
    # Find the indices of the start and end of the transients
    trace_logic = trace.copy()
    trace_logic[trace_logic>0] = 1

    start_index = np.where(np.diff(trace_logic) == 1)[0]
    end_index = np.where(np.diff(trace_logic) == -1)[0]


    # Create an empty array of zeros
    trace_length = trace.shape[0]
    shuff_trace = np.zeros(trace_length)

    # Function to check for overlap
    def is_overlap(start_pos, block_length, occupied_positions):
        # Check if any of the positions in the block are already occupied
        return any(pos in occupied_positions for pos in range(start_pos, start_pos + block_length))

    # To track the positions already occupied by blocks
    occupied_positions = set()  

    # Insert the blocks into the array without overlap
    # for str, end in zip(start_index, end_index):
        
    #     transient_len = end - str
    #     inserted = False

    #     while not inserted:
    #         # Randomly choose a start position
    #         start_pos = np.random.randint(0, trace_length - transient_len)
            
    #         # Check if the block can fit without overlapping
    #         if not is_overlap(start_pos, transient_len, occupied_positions):
    #             # Insert the block
    #             shuff_trace[start_pos:start_pos + transient_len] = trace[str:end]
    #             # Mark these positions as occupied
    #             occupied_positions.update(range(start_pos, start_pos + transient_len))
    #             inserted = True

    max_attempts = 1000  # Prevent infinite loops

    for str, end in zip(start_index, end_index):
        transient_len = end - str
        inserted = False
        attempts = 0
    
        while not inserted and attempts < max_attempts:
            start_pos = np.random.randint(0, trace_length - transient_len)
            if not is_overlap(start_pos, transient_len, occupied_positions):
                shuff_trace[start_pos:start_pos + transient_len] = trace[str:end]
                occupied_positions.update(range(start_pos, start_pos + transient_len))
                inserted = True
            attempts += 1
    
        if not inserted:
            print(f"Warning: Could not insert transient of length {transient_len} after {max_attempts} attempts.")


    return shuff_trace

import numpy as np


    
def sort_ratemaps(rates):
    """
    Sorts rate maps by the bin with the highest rate in each row. Returns the sorted 
    indices and the corresponding sorted rate maps.

    Parameters:
    rates (pd.DataFrame): A DataFrame where each row represents a rate map, 
                          and columns represent different bins.

    Returns:
    tuple:
        - sort_indx (pd.Index): Indices of rate maps sorted by the highest rate bin.
        - sorted_ratemaps (pd.DataFrame): Rate maps sorted by the highest rate bin.
    """

    # get number of ratemaps
    num_ratemaps = rates.shape[0]
    
    # iterate through each ratemap
    maxrate_bin = np.zeros(num_ratemaps, dtype=int)
    for i in range(num_ratemaps):
        # find the bin with the maximum rate
        max_bin = np.argmax(rates.iloc[i, :])
        #store bin index in maxrate_bin
        maxrate_bin[i] = max_bin

    # sort maxrate_bins in ascending order
    id = np.argsort(maxrate_bin)[::-1]
    
    # retrieve the sorted row indices from the original DataFrame
    sort_indx = rates.index[id]

    # sort the original ratemaps by max bin
    sorted_ratemaps = rates.reindex(sort_indx).fillna(0)

    return np.array(sort_indx), np.array(sorted_ratemaps)

def calculate_placecell_pcnt(p_values, alpha = 0.05):
    
    p_val_df = pd.DataFrame(p_values)
    place_cell_count = int(p_val_df[p_val_df<=alpha].count()[0])

    total_cells = p_val_df.shape[0]
    pc_pcnt = (place_cell_count/total_cells) * 100

    return print(f"{place_cell_count} out of {total_cells} are place cells: {pc_pcnt:.2f}%")