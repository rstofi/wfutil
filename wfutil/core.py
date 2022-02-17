"""The core functions and definition of the WFdata class which holds the 
time-frequency data and the associated axes

The idea is to have a structured format to work with data in the
time-frequency-polarisation space. The data can be e.g. (complex) gains,
or flag-related data.

As such I try to make a data model using both a pythonic and more 'physical'
mindset.

Similar solutions might exist, but I have not found a working data model that
could handle masks and complex numbers at the same time... maybe I am wrong here,
but after a quick look at existing structures like xarray, I decided it might be
better to make m own minimalist solution for my specific needs.

TO DO: add logging

TO DO: add a function that creates the mapping index given a merged 1D array and
        a 1D array to be mapped => check maybe numpy alread has this...

"""

__all__ = ['WFdata', 'load_WFdata', 'map_and_merge_1D', 'map_data_arrays']


import os
import copy
import numpy as np
import numpy.ma as ma

#=== Globals ===
global _WF_DESC_COMPULSORY

_WF_DESC_COMPULSORY = ['Masked', 'DataUnit', 'TimeUnit', 'ChanUnit', 'PolUnit',
                        'DataFrame', 'TimeFrame', 'ChanFrame', 'PolFrame']

global _SUPPORTED_DATA_FRAMES
global _SUPPORTED_TIME_FRAMES
global _SUPPORTED_CAN_FRAMES
global _SUPPORTED_POL_FRAMES


_SUPPORTED_DATA_FRAMES = ['CGain','Flag']
_SUPPORTED_TIME_FRAMES = ['MJD', 'UNIX']
_SUPPORTED_CAN_FRAMES = ['freq', 'vel']
_SUPPORTED_POL_FRAMES = ['XY', 'Stokes', 'Intensity']

#*******************************************************************************
#=== Functions I ===
def map_and_merge_1D(arrays_to_map):
    """Core function to map a list of 1D arrays onto a common (sorted) array.

    This is the main function to join data sets and should work on both the
    time and channel axes/arrays

    NOTE: it should not used on the polarisation axis for obvious reasons...

    The routine used is maybe slow, but is general and can handle:

        - arrays that are not sorted
        - arrays that have non monotonly increasing values (e.g. gaps)
        - arrays with overlapping values (currently these are mapped on the
            same index on the merged_array)

    Parameters
    ==========
    arrays_to_map: list of <numpy.NdArray>
        This is a list of 1D numpy arrays to merge. Each array should be real-valued.

    Returns
    =======
    merged_array: <numpy.NdArray>
        The merged and sorted 1D array, i.e. the unique union of the input arrays

    mapped_index_list: list of <numpy.NdArray>
        A list containing the mapping iundices for each input array. i.e this is
        the same size (containing same sized arrays) as `arrays_to_map` withe ach
        array containing the indices of which the values are mapped on the
        `merged_array`

    """
    mapped_index_list = []

    #The np.unique() returns a sorted (!) unique list
    #It also returns the mapping
    merged_array, mapped_indices = np.unique(np.concatenate(arrays_to_map),
                                            return_inverse=True)

    #Non optimal way to create the mapping to the input format
    i = j = 0
    for arr in arrays_to_map:
        arr_size = np.size(arr)
        j += arr_size

        index_map = mapped_indices[i:j]
        mapped_index_list.append(index_map)

        i = j

    return merged_array, mapped_index_list

def map_data_arrays(merged_data_array,
                    merged_time_array,
                    merged_chan_array,
                    data_arrays_to_merge_list,
                    time_mapped_index_list,
                    chan_mapped_index_list):

    """The core merging function for the data array polarisation slices.

    Note that merging along the polarisation axis (i.e. merge data with different
    polarisaton frames) is meaningless. And converting between polarisation frames
    is out of the scope of this package
    
    This function maps input 2D time-frequency arrays of [Ni x Mi]
    size to the [N x M] sized input array, which can already have some values.

    The N sized time and M sized channel arrays are the merged arrays and the
    mapped index list (e.g. created by map_and_merge_1D()) are the indices
    used for mapping

    merged_data_array: np 2d array [N x M] size
    merged_time_array: np 1d array N size
    merged_chan_array: np 1d array M size

    the rest are lists of the same length with each element:

    data_arrays_to_merge_list[i]: np 2d array [n x m] size
    time_mapped_index_list[i]: 1d np array of indices with n size
    chan_mapped_index_list[i]: 1d np array of indices with m size


    NOTE: if the input data arrays are mapped to the same indice taht indiex value
        is overwritten, and so the final value will be the one given last
    
    NOTE: If the input `merged_data_array` is masked, then the output will be
        masked as well, however if the input array is not masked, but any of the
        arrays in the `data_arrays_to_merge_list` are masked, the mask *will be 
        ignored*!

    #TO DO: vectorise as much as possible from the marging loops

    Parameters
    ==========
    merged_data_array: <numpy.NdArray> or <numpy.ma.core.MaskedArray>
        The 2D array with [N x M] size to wich the other data arrays will be
        mapped. It can be masked.

    merged_time_array: <numpy.NdArray>
        The time axis array with [N x 1] size

    merged_chan_array: <numpy.NdArray>
        The channel axis array with [M x 1] size        

    data_arrays_to_merge_list: list of <numpy.NdArray> or <numpy.ma.core.MaskedArray>
        The list of the data arrays to be mapped and merged onto the `merged_data_array`

        The ith element of the array is a 2D array with [n x m] size which can be
        masked.

        NOTE: if the `merged_data_array` is not masked these masks will be ignored

    time_mapped_index_list: list of <numpy.NdArray>
        The time axis map indices (e.g. generated by `map_and_merge_1D`) for 
        each data to be merged stored in a list.

        The ith element of the array is a 1D array with n length.

    chan_mapped_index_list: list of <numpy.NdArray>
        The channel axis map indices (e.g. generated by `map_and_merge_1D`) for 
        each data to be merged stored in a list.

        The ith element of the array is a 1D array with m length.


    Returns
    =======
    merged_data_array: <numpy.NdArray> or <numpy.ma.core.MaskedArray>
        The `merged_data_array` that is appended with the data mapped.
    """
    #Check for masking
    if not isinstance(merged_data_array, np.ma.MaskedArray):
            masked = False

    else:
        masked = True
        #Create the master mask for merging and mapping the masks
        master_mask = copy.deepcopy(ma.getmaskarray(merged_data_array))

    #Some sanity checks
    if np.shape(merged_data_array)[0] != np.size(merged_time_array):
        raise ValueError('The merged time array and input data array time \
dimension sizes are different!')

    if np.shape(merged_data_array)[1] != np.size(merged_chan_array):
        raise ValueError('The merged chan array and input data array chan \
dimension sizes are different!')

    if len(data_arrays_to_merge_list) != len(time_mapped_index_list):
        raise ValueError('Input data array list and time index array list are not equal!')

    if len(data_arrays_to_merge_list) != len(chan_mapped_index_list):
        raise ValueError('Input data array list and chan index array list are not equal!')

    # The mapping loop
    for i in range(0,len(data_arrays_to_merge_list)):

        #Some size checks 
        if np.shape(data_arrays_to_merge_list[i])[0] != np.size(time_mapped_index_list[i]):
            raise ValueError('The merged time index array and input data array time \
dimension sizes are different for the {0:d} element of the input list'.format(i))

        if np.shape(data_arrays_to_merge_list[i])[1] != np.size(chan_mapped_index_list[i]):
            raise ValueError('The merged chan index array and input data array chan \
dimension sizes are different for the {0:d} element of the input list'.format(i))

        #Separate the mask and value matrices if the input is masked
        if masked:
            if isinstance(data_arrays_to_merge_list[i], np.ma.MaskedArray):
                #Initialise the mapped mask to map the mask onto
                mapped_mask = np.full((np.size(merged_time_array),
                                        np.size(merged_chan_array)),
                                        False, dtype=bool)
            
                #Get the mask
                data_array_mask = ma.getmaskarray(data_arrays_to_merge_list[i])

                #This is just a namespace for convinience
                data_array_values = data_arrays_to_merge_list[i]

            else:
                #TO DO: raise warning and fix exception
                #raise ValueError('Input data array {0:d} is not masked!'.format(i))
                raise Warning('Input data array {0:d} is not masked!'.format(i))

                #Use an empty mask
                data_array_mask = np.full((np.size(time_mapped_index_list[i]),
                                            np.size(chan_mapped_index_list[i])),
                                        False, dtype=bool)

                #This is needed so mask joining is working
                data_array_values = ma.masked_array(data_arrays_to_merge_list[i],
                                                    data_array_mask)
        else:
            data_array_values = data_arrays_to_merge_list[i]

        #TO DO: vectorise as much as possible from the double loop below

        #Loop through time
        for t in range(0,np.size(time_mapped_index_list[i])):
            #Loop trough channels
            for c in range(0,np.size(chan_mapped_index_list[i])):
                #Map and merge the data arrays
                merged_data_array[time_mapped_index_list[i][t],
                                chan_mapped_index_list [i][c]] = \
                                ma.getdata(data_array_values[t,c])

                #Map and merge the masks
                if masked:
                    mapped_mask[time_mapped_index_list[i][t],
                                chan_mapped_index_list [i][c]] = \
                                data_array_mask[t,c]

        #Create the joint master mask and current mask
        if masked:
            master_mask = ma.mask_or(master_mask, mapped_mask)
            #merged_data_array = ma.masked_array(merged_data_array,master_mask)

    #Apply the master (merged) mask to the merged data
    if masked:
        merged_data_array = ma.masked_array(merged_data_array,master_mask)

    return merged_data_array

#*******************************************************************************
#=== Classes ===
class WFdata(object):
    """The main data model: data in the time-frequency(-polarisation) space.

    The data model consist of:

        - [P x N x M] sized data array
        - [P x 1] sized polarisation axis array (this is an array of strings)
        - [N x 1] sized time axis array
        - [M x 1] sized channel (frequency) axis array
        - dictionary with the units and frames

    This is basically a structured wrapper around these data arrays and associated
    info.

    
    NOTE: only the data array can be masked, and if masked axis arrays are provided,
        then the mask is removed.

    TO DO: add optional but automatic masking for the axis arrays

    Keyword Arguments
    =================

    data_array: <numpy.NdArray> or <numpy.ma.core.MaskedArray>
        [P x N x M] sized data array can be real, complex and it also can be a
        masked array. The shape has to be (pol, time, chan)

    pol_array: <numpy.NdArray>
        [P x 1] sized polariastion axis array.

    time_array: <numpy.NdArray>
        [N x 1] sized time axis array.

    pol_array: <numpy.NdArray>
        [N x 1] sized channel axis array.

    wf_desc: dict, optional
        A dictionary containing the data units and the axis frames. It can be used
        add additional information to the WFdata object.
        It has to have the following keys (these should be self-explanatory):

        ['Masked', DataUnit', 'TimeUnit', 'ChanUnit', 'PolUnit',
        'DataFrame', 'TimeFrame', 'ChanFrame', 'PolFrame']

        The walid key walues are only defined for the array frames as the units
        *should* be derived from them (although the units shoul be stored properly
        as well)

        _SUPPORTED_DATA_FRAMES = ['CGain','Flag']
        _SUPPORTED_TIME_FRAMES = ['MJD', 'UNIX']
        _SUPPORTED_CAN_FRAMES = ['freq', 'vel']
        _SUPPORTED_POL_FRAMES = ['XY', 'Stokes', 'Intensity']
    
        NOTE: this is an optional parameter, and if not specified, the default 
        keys and values will be used and so the metadata can be corrupted this way

        TO DO: add an automatic check for units and conversions when one makes a
                frame conversion.
    """
    def __init__(self,
                data_array,
                pol_array,
                time_array,
                chan_array,
                wf_desc = {'Masked' : 'False',
                    'DataUnit' : 'amp',
                    'TimeUnit' : 's',
                    'ChanUnit' : 'Hz',
                    'PolUnit' : 'Cm^-2',
                    'DataFrame' : 'CGain',
                    'TimeFrame' : 'UNIX',
                    'ChanFrame' :'freq',
                    'PolFrame': 'XY'}):

        #Set up instances
        self.data_array = data_array
        self.pol_array = pol_array
        self.time_array = time_array
        self.chan_array = chan_array
        self.wf_desc = wf_desc

        #Derived instances
        self.axes = [self.pol_array, self.time_array, self.chan_array]


        #=== Check the wf desc dictionary values
        for WF_DESC in _WF_DESC_COMPULSORY:
            if WF_DESC not in set(self.wf_desc.keys()):
                raise ValueError('The key: {0:s} is not defined in the \
WF description (wf_descr)!'.format(WF_DESC))

        if self.wf_desc['DataFrame'] not in _SUPPORTED_DATA_FRAMES:
            raise TypeError('The data frame: {0:s} is not supported!'.format(
                            self.wf_desc['DataFrame']))

        if self.wf_desc['TimeFrame'] not in _SUPPORTED_TIME_FRAMES:
            raise TypeError('The time frame: {0:s} is not supported!'.format(
                            self.wf_desc['TimeFrame']))

        if self.wf_desc['ChanFrame'] not in _SUPPORTED_CAN_FRAMES:
            raise TypeError('The channel frame: {0:s} is not supported!'.format(
                            self.wf_desc['ChanFrame']))

        if self.wf_desc['PolFrame'] not in _SUPPORTED_POL_FRAMES:
            raise TypeError('The polarisation frame: {0:s} is not supported!'.format(
                            self.wf_desc['PolFrame']))

        #=== Check the data dimensions
        self.data_shape = np.shape(self.data_array)

        #TO do add special case for polarisation handling for Intensity (i.e. no pol)
        if self.data_shape[0] != np.size(self.pol_array):
            raise ValueError('Data array pol dimension and pol array dimension are not equal!')

        if self.data_shape[1] != np.size(self.time_array):
            raise ValueError('Data array time dimension and time array dimension are not equal!')

        if self.data_shape[2] != np.size(self.chan_array):
            raise ValueError('Data array chan dimension and chan array dimension are not equal!')

        #=== Check for masks
        
        #NOTE: this part of the code is not tested

        if isinstance(self.data_array, np.ma.MaskedArray):
            self.wf_desc['Masked'] = 'True'
        else:
            self.wf_desc['Masked'] = 'False'

        #Convert possible masked axes to non-masked values
        for WFax in self.axes:
            if isinstance(WFax, np.ma.MaskedArray):
                WFax = ma.getdata(Wfax) #This converts the ma array to numpy array in theory
            else:
                pass

    #=== Functions ===
    def get_pol_slice(self,pol_val):
        """Routine to sub-select a polariastion from the `data_array` 
        """
        print(np.where(self.pol_array == pol_val))

        #HERE IS THE BUG I NEED TO WORK ON...

        p = np.where(self.pol_array == pol_val)[0][0]
        
        print(p)
        pol_slice_data_array = self.data_array[p,...]

        return pol_slice_data_array

    def save_WFdata(self, output_path, WFdata_name, overwrite=True, compressed=True):
        """Write the WFdata to disc. The basic format is a .npz file, i.e. a
        binary format for labeled numpy arrays.

        Compression is set by default, but optional.

        NOTE: numpy cannot handle masked arrays with savez, e.g. see:
                
                https://github.com/numpy/numpy/issues/18134 

        this was one of the reasons for implementing the WFutil package...

        As such this code supports this option. If the data array is masked the
        mask is saved as an additional array, and it has to be applied when reading
        the data from disc
        
        Parameters
        ==========
        output_path: str
            Folder path where the `WFdata` opbject is saved

        WFdata_name: str
            Name of the file to be saved *without* the `.npz` extension. It will
            be added to the file automatically

        overwrite: bool, optional
            If True an existing file will be overwritten, otherwise the function
            halts if the file already exists

        compressed: bool, optional
            If True, the data will be saved ina  compressed format.
            For more info see: `numpy.savez()`
        
        Returns
        =======
        Save the `WFdata` to a `.npz` binary file

        """
        #The save routine will automatically azz an .npz extension to the name
        outfile = os.path.join(output_path, WFdata_name)

        if os.path.isdir(outfile) and overwrite == False: 
            raise TypeError('WFdata file already exist, and the \
overwrite parameters is set to False!')

        #Convert dict to numpy array:
        # This will produce a  [[key1, value1], ... [keyN, valueN]] array
        wf_desc_array = np.array(list(self.wf_desc.items()))

        #Save the data
        if self.wf_desc['Masked']:
            #Create the mask and retrieve the original array
            if compressed:
                np.savez_compressed(outfile,
                        data_array = ma.getdata(self.data_array),
                        dmask_array = ma.getmaskarray(self.data_array), #see documentation
                        time_array = self.time_array,
                        chan_array = self.chan_array,
                        pol_array = self.pol_array,
                        wf_desc = wf_desc_array)
            else:
                np.savez(outfile,
                        data_array = ma.getdata(self.data_array),
                        dmask_array = ma.getmaskarray(self.data_array), #see documentation
                        time_array = self.time_array,
                        chan_array = self.chan_array,
                        pol_array = self.pol_array,
                        wf_desc = wf_desc_array)

        else: 
            #Save to a numpy binary file
            if compressed:
                np.savez_compressed(outfile,
                        data_array = self.data_array,
                        time_array = self.time_array,
                        chan_array = self.chan_array,
                        pol_array = self.pol_array,
                        wf_desc = wf_desc_array)
            else:            
                np.savez(outfile,
                        data_array = self.data_array,
                        time_array = self.time_array,
                        chan_array = self.chan_array,
                        pol_array = self.pol_array,
                        wf_desc = wf_desc_array)


        #TO DO: add slicing routine

#*******************************************************************************
#=== Functions II ===
def load_WFdata(input_path, input_fname):
    """Routine to read and create a `WFdata` object from an `.npz` filr, which 
    was created by the `save_WFdata()` function.

    Parameters
    ==========
    input_path: str
        Folder path where the `WFdata` opbject is saved

    input_fname: str
        Name of the file to be read *without* the `.npz` extension. It will
        be added to the file automatically

    Returns
    =======
    Create an `WFdata` object

    """
    #Append WFdata_name
    input_fname += '.npz'

    input_WFname = os.path.join(input_path, input_fname)
    raw_WFdata = np.load(input_WFname, mmap_mode='r')

    #Convert numpy array to dict
    wf_desc_dict = {}
    raw_wf_desc_array = raw_WFdata['wf_desc']
    for k, v in zip(raw_wf_desc_array[:,0], raw_wf_desc_array[:,1]):
        wf_desc_dict[k] = v

    #Check if mask is saved and apply it if yes
    input_WF_keys = [k for k in raw_WFdata.files]

    if 'dmask_array' in input_WF_keys:
        data_array = ma.masked_array(raw_WFdata['data_array'],
                                    raw_WFdata['dmask_array'])
    else:
        data_array = raw_WFdata['data_array']

    wf_obj = WFdata(data_array = data_array,
                    time_array = raw_WFdata['time_array'],
                    chan_array = raw_WFdata['chan_array'],
                    pol_array = raw_WFdata["pol_array"],
                    wf_desc = wf_desc_dict)

    #Close the file to prevent memory leak
    raw_WFdata.close()

    return wf_obj

def merge_WFdata(WFdata_to_merge_list):
    """The top level function to merge several `WFdata` object and create a merged
    object.

    TO DO: add some more documentation to this docstring

    """

    #=== Create the merged time and channel and polarisation axis
    #NOTE: this step can be a bottleneck if the arrays are the same (e.g. only 
    #   time axis merging is needed but the channels are identical)


    #Get the polarisation axis
    #Here I use the first `WFdata` object but all should have the same pol dimensions!
    merged_pol_array = copy.deepcopy(WFdata_to_merge_list[0].pol_array)

    #Get the merged axes
    time_arrays_to_map = []
    chan_arrays_to_map = []

    for WFD in WFdata_to_merge_list:
        #Check for polarisation sizes
        if np.size(WFD.pol_array) != np.size(merged_pol_array):
            raise ValueError('The {0:d} th WFdata has incompatible polarisation \
with size {1:d} and pol frame of {2:s}!'.format(i, np.size(WFD.pol_array),
                                                WFD.wf_desc['PolFrame']))

        #Otherwise append the time and channel arrays to map
        else:
            time_arrays_to_map.append(WFD.time_array)
            chan_arrays_to_map.append(WFD.chan_array)

    merged_time_array, mapped_time_indices = map_and_merge_1D(time_arrays_to_map)
    merged_chan_array, mapped_chan_indices = map_and_merge_1D(chan_arrays_to_map)


    #TO DO: add a check for the data frames and units for the different `WFdata`

    #=== Create the empty merged data array and copy the description dictionary

    #Copy the first `WFdata` desc dictionary
    merged_wf_desc = copy.deepcopy(WFdata_to_merge_list[0].wf_desc)

    #Check if the first `WFdata` data_array is masked or not
    #The type is determined by numpy ma => see the documentation

    #TO DO: perform a data type check

    if WFdata_to_merge_list[0].wf_desc['Masked'] == 'True':
        merged_data_array = \
        np.ma.array(np.zeros((np.size(merged_pol_array),
                            np.size(merged_time_array),
                            np.size(merged_chan_array))),
        mask=False)

    else:
        merged_data_array = \
        np.zeros(np.shape((np.size(merged_pol_array),
                            np.size(merged_time_array),
                            np.size(merged_chan_array))))

    #=== Create an empty WFdata object
    merged_WFdata = WFdata(data_array = merged_data_array,
                time_array = merged_time_array,
                chan_array = merged_chan_array,
                pol_array = merged_pol_array,
                wf_desc = merged_wf_desc)

    #=== Map each object onto the merged object


    #TO DO: put this into a distictive function and optimize it more...

    #Loop through the polarisation
    for pol_val in merged_WFdata.pol_array:
        #Get the pol data slice
        
        merged_data_array_pol_slice = merged_WFdata.get_pol_slice(pol_val)

        print(merged_data_array_pol_slice)

        #Get the list of the arrays
        data_arrays_to_map = []

        for WFD in WFdata_to_merge_list:
            WFD_data_array_pol_slice = WFD.get_pol_slice(pol_val)
            data_arrays_to_map.append(WFD_data_array_pol_slice)


        merged_data_array_pol_slice = map_data_arrays(
                    merged_data_array = merged_data_array_pol_slice,
                    merged_time_array = merged_WFdata.time_array,
                    merged_chan_array = merged_WFdata.chan_array,
                    data_arrays_to_merge_list = data_arrays_to_map,
                    time_mapped_index_list=mapped_time_indices,
                    chan_mapped_index_list=mapped_chan_indices)

        #Now owerwrite the merged data array
        p = np.where(merged_WFdata.pol_array == pol_val)[0]
        merged_WFdata.data_array[p,...] = merged_data_array_pol_slice

    return merged_WFdata

#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass

#*******************************************************************************
#=== Quick and dirty testing ===

#TO DO write real testing module

#Quick and dirty testing
test_saving = False
test_map_and_merge = False
test_map_data_arrays = False
test_merge_WFdata = True

if test_merge_WFdata:

    #First WFData
    d1 = np.array([[[1,0],[0,1]],[[0,1],[1,0]]])
    m1 = np.array([[[True,False],[False,False]],
                    [[False,True],[False,False]]])

    d1 = ma.masked_array(d1,m1)

    p1 = ['a','b']
    t1 = np.array([1,2])
    c1 = np.array([1,3])

    WFd1 = WFdata(data_array = d1,
                time_array = t1,
                chan_array = c1,
                pol_array = p1)

    #Second WFData
    d2 = np.array([[[1,0],[0,1]],[[0,1],[1,0]]])
    m2 = np.array([[[True,False],[False,False]],
                    [[False,True],[False,False]]])

    d2 = ma.masked_array(d2,m2)

    p2 = ['a','b']
    t2 = np.array([3,4])
    c2 = np.array([2,3])

    WFd2 = WFdata(data_array = d2,
                time_array = t2,
                chan_array = c2,
                pol_array = p2)

    #Merge them
    WFd_to_merge_list = [WFd1, WFd2]

    merged_WFd = merge_WFdata(WFd_to_merge_list)

    print(merged_WFd.pol_array)
    print(merged_WFd.time_array)
    print(merged_WFd.chan_array)
    print(merged_WFd.data_array)


#Tests
if test_map_data_arrays:
    
    d1 = np.array([[1,3,1],[1,1,1]])
    t1 = np.array([2,1])
    c1 = np.array([10,20,30])

    d1_mask = np.array([[True, False, False],[False, False, False]])

    d1 = ma.masked_array(d1,d1_mask)

    d2 = np.array([[2,2],[2,0]])
    t2 = np.array([0,2])
    c2 = np.array([10,20])

    d2_mask = np.array([[False, False],[False, True]])

    d2 = ma.masked_array(d2,d2_mask)

    d3 = np.array([[4],[5],[4]])
    t3 = np.array([0,2,4])
    c3 = np.array([10])

    d3_mask = np.array([[False],[False],[False]]) 

    d3 = ma.masked_array(d3,d3_mask)

    data_arrays_to_map = [d1,d2,d3]
    time_arrays_to_map = [t1,t2,t3]
    chan_arrays_to_map = [c1,c2,c3]

    #Create merged axes
    merged_time, mapped_time_indices = map_and_merge_1D(time_arrays_to_map)
    merged_chan, mapped_chan_indices = map_and_merge_1D(chan_arrays_to_map)

    #print(mapped_time_indices)
    #print(mapped_chan_indices)

    #Create merged empty data array
    merged_darr = ma.masked_array(np.zeros((np.size(merged_time),np.size(merged_chan))))

    #print(merged_time)
    #print(merged_chan)
    #print(merged_darr)

    merged_darr = map_data_arrays(merged_data_array=merged_darr,
                    merged_time_array=merged_time,
                    merged_chan_array=merged_chan,
                    data_arrays_to_merge_list=data_arrays_to_map,
                    time_mapped_index_list=mapped_time_indices,
                    chan_mapped_index_list=mapped_chan_indices)

    print(merged_darr)


if test_map_and_merge:
    a = np.array([1,2,5])
    b = np.array([2,3])
    c = np.array([4,1])

    arrays_to_map = [a,b,c]

    merged_array, mapped_indices = map_and_merge_1D(arrays_to_map)

    print(type(merged_array), merged_array)
    print(type(mapped_indices), mapped_indices)

#Test saving features
if test_saving:
    working_dir = './'
    test_file_name = 'test'

    datarr = np.array([[[1,1],[1,0],[0,0]],
                        [[1,0],[0,1],[0,0]]],dtype=complex)

    dmask = np.array([[[True,False],[False,False],[False,False]],
            [[False,False],[False,False],[False,True]]], dtype=bool)

    masked_datarr = ma.masked_array(datarr,dmask)
    #masked_datarr = datarr

    tarr = np.array([0,1,2])
    farr = np.array([0,1])
    polarr = np.array(['XX','YY'])

    print(type(datarr[0,0]))

    a = WFdata(data_array = masked_datarr,
                time_array = tarr,
                chan_array = farr,
                pol_array = polarr)

    print(type(a.data_array))

    a.save_WFdata(working_dir,test_file_name,overwrite=True)

    b = load_WFdata(working_dir, test_file_name)

    print(type(b.data_array))
    print(type(b.data_array[0,0,0]))
    print(type(b.data_array[0,0,1]))

    print(np.shape(b.data_array))

    print(b.data_array[0,...])

    print(ma.getdata(b.data_array[0,...]))
