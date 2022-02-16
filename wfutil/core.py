"""The core functions and definition of the WFdata class which holds the 
time-frequency data and the associated axes

NOTE: currently only 2D data arrays are supported
"""

__all__ = ['WFdata', 'load_WFdata', 'map_and_merge_1D', 'map_data_arrays']


import os
import copy
import numpy as np
import numpy.ma as ma

#=== Globals ===
global _WF_DESC_COMPULSORY

_WF_DESC_COMPULSORY = ['DataUnit', 'TimeUnit', 'ChanUnit', 'PolUnit',
                        'DataFrame', 'TimeFrame', 'ChanFrame', 'PolFrame']

global _SUPPORTED_DATA_FRAMES
global _SUPPORTED_TIME_FRAMES
global _SUPPORTED_CAN_FRAMES
global _SUPPORTED_POL_FRAMES


_SUPPORTED_DATA_FRAMES = ['CGain','Flag']
_SUPPORTED_TIME_FRAMES = ['MJD', 'UNIX']
_SUPPORTED_CAN_FRAMES = ['freq', 'vel']
_SUPPORTED_POL_FRAMES = ['XY', 'Stokes', 'Intensity']


#=== Functions ===
def decompose_complex_matrix(cmatrix):
    """

    """
    cmatrix_shape = np.shape(cmatrix)

    decomposed_cmatrix = np.zeros((cmatrix_shape[0],
                                cmatrix_shape[1],
                                2))

    decomposed_cmatrix[...,0] = np.real(cmatrix)
    decomposed_cmatrix[...,1] = np.imag(cmatrix)

    return decomposed_cmatrix

def compose_complex_matrix(decomposed_cmatrix):
    """

    """
    decomposed_cmatrix_shape = np.shape(decomposed_cmatrix)

    cmatrix = np.zeros((decomposed_cmatrix_shape[0],
                        decomposed_cmatrix_shape[1]),
                        dtype=complex)

    cmatrix = np.add(decomposed_cmatrix[...,0]) #Add real part
    cmatrix = np.add(np.multiply(1j, decomposed_cmatrix[...,1])) #Add imaginary part

    return cmatrix

#=== Classes ===
class WFdata(object):
    """

    """
    def __init__(self,
                data_array,
                time_array,
                chan_array,
                pol_array,
                wf_desc = {'DataUnit' : 'amp',
                    'TimeUnit' : 's',
                    'ChanUnit' : 'Hz',
                    'PolUnit' : 'Cm^-2',
                    'DataFrame' : 'CGain',
                    'TimeFrame' : 'UNIX',
                    'ChanFrame' :'freq',
                    'PolFrame': 'XY'}):
        """

        Data array shape:

        (pol, time, chan)

        this is how the input data must look like!

        #Note that the data array is the only one should be masked!

        """

        self.data_array = data_array
        self.time_array = time_array
        self.chan_array = chan_array
        self.pol_array = pol_array
        self.wf_desc = wf_desc


        #Check the wf desc dictionary values
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

        #Check the data dimensions
        self.data_shape = np.shape(self.data_array)

        #TO do add special case for polarisation handling for Intensity (i.e. no pol)
        if self.data_shape[0] != np.size(self.pol_array):
            raise ValueError('Data array pol dimension and pol array dimension are not equal!')

        if self.data_shape[1] != np.size(self.time_array):
            raise ValueError('Data array time dimension and time array dimension are not equal!')

        if self.data_shape[2] != np.size(self.chan_array):
            raise ValueError('Data array chan dimension and chan array dimension are not equal!')

    def save_WFdata(self, output_path, WFdata_name, overwrite=True, compressed=True):
        """

        """
        #The save routine will automatically azz an .npz extension to the name
        outfile = os.path.join(output_path, WFdata_name)

        if os.path.isdir(outfile) and overwrite == False: 
            raise TypeError('WFdata file already exist, and the \
overwrite parameters is set to False!')

        #Convert dict to numpy array:
        # This will produce a  [[key1, value1], ... [keyN, valueN]] array
        wf_desc_array = np.array(list(self.wf_desc.items()))

        #Numpy cannot handle masked arrays with savez:
        #See: https://github.com/numpy/numpy/issues/18134

        #Handle this case:

        #If the data array is masked
        if isinstance(self.data_array, np.ma.MaskedArray):
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

#=== Functions II ===
def load_WFdata(input_path, input_fname):
    """

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

    #Check for saved mask
    #TO DO add better documentation

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

    raw_WFdata.close()

    return wf_obj


def map_and_merge_1D(arrays_to_map):
    """Core function to map a list of 1D arrays onto a common (sorted) array.

    This is a core function to join data sets and should work on both the
    time_array and the chan_array

    This is maybe slow, but is general and can handle:

        - arrays that are not sorted
        - arrays that have non monotonly increasing values (e.g. gaps)
        - arrays with overlapping values (currently these are mapped on the
            same index on the merged_array)

    input should be a list and each element should be an array with real values

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
                    chan_mapped_index_list,
                    masked = False):

    """
    The core merging function: it maps input 2D time-frequency arrays of [Ni x Mi]
    size to the [N x M] sized input array (that could already have some values!)

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


    NOTE: if the input data arrays are mapped to the same indice it is overwritten,
            and so the final value will be the one given last

            + if the input data array has a value there it is overwritten as well

    Masked is a quick and dirty option to joint he masks instead of using the last
    mask provided. This is a bug if masked arrays are provided

    """
    #Some sanity checks
    if masked:
        if not isinstance(merged_data_array, np.ma.MaskedArray):
            raise TypeError('The input merged data array is not a masked array,\
but mask merging is specified!')
        else:
            master_mask = copy.deepcopy(ma.getmaskarray(merged_data_array))

            #print(master_mask)

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

    #Now the big mapping loop
    for i in range(0,len(data_arrays_to_merge_list)):

        #Some size checks 
        if np.shape(data_arrays_to_merge_list[i])[0] != np.size(time_mapped_index_list[i]):
            raise ValueError('The merged time index array and input data array time \
dimension sizes are different for the {0:d} element of the input list'.format(i))

        if np.shape(data_arrays_to_merge_list[i])[1] != np.size(chan_mapped_index_list[i]):
            raise ValueError('The merged chan index array and input data array chan \
dimension sizes are different for the {0:d} element of the input list'.format(i))

        if masked:
            #HERE the masks get tangled up

            #To do: build the mask and the value array separately and apply the mask in the end

            #i.e. when the masked arrays are merged the actual values not written and
            # when the actual values are written, then the masks are not written

            #Soultion do the basic stuff on a non-masked array and build a mask if needed in parallel and apply it in the end


            if isinstance(data_arrays_to_merge_list[i], np.ma.MaskedArray):
                mapped_mask = np.full((np.size(merged_time_array),
                                        np.size(merged_chan_array)),
                                        False, dtype=bool)
            
                data_array_mask = ma.getmaskarray(data_arrays_to_merge_list[i])

                #Thois is just passing the name but using actually a masked object
                data_array_values = data_arrays_to_merge_list[i]

            else:
                #TO DO: raise warning and fix exception
                #raise ValueError('Input data array {0:d} is not masked!'.format(i))
                raise Warning('Input data array {0:d} is not masked!'.format(i))

                data_array_mask = np.full((np.size(time_mapped_index_list[i]),
                                            np.size(chan_mapped_index_list[i])),
                                        False, dtype=bool)

                #This is needed so mask joining is working
                data_array_values = ma.masked_array(data_arrays_to_merge_list[i],
                                                    data_array_mask)
        else:
            data_array_values = data_arrays_to_merge_list[i]

        #TO DO vectorise as much as possible from the double loop below

        #Loop through time
        for t in range(0,np.size(time_mapped_index_list[i])):
            #Loop trough channels
            for c in range(0,np.size(chan_mapped_index_list[i])):
                merged_data_array[time_mapped_index_list[i][t],
                                chan_mapped_index_list [i][c]] = \
                                data_array_values[t,c]

                if masked:
                    mapped_mask[time_mapped_index_list[i][t],
                                chan_mapped_index_list [i][c]] = \
                                data_array_mask[t,c]

        #Apply the joint master mask and current mask
        if masked:
            master_mask = ma.mask_or(master_mask, mapped_mask)
            merged_data_array = ma.masked_array(merged_data_array,master_mask)

        #For quick and dirty testing for masked arrays
        print(merged_data_array)
        #print(type(merged_data_array))
        #For masked arrays case
        print(ma.getdata(merged_data_array))

#=== MAIN ===

#TO DO write real testing module

#Quick and dirty testing
test_saving = False
test_map_and_merge = False
test_map_data_arrays = True

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

    d3 = np.array([[4],[0],[4]])
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

    map_data_arrays(merged_data_array=merged_darr,
                    merged_time_array=merged_time,
                    merged_chan_array=merged_chan,
                    data_arrays_to_merge_list=data_arrays_to_map,
                    time_mapped_index_list=mapped_time_indices,
                    chan_mapped_index_list=mapped_chan_indices,
                    masked=True)


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
