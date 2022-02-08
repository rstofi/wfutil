"""The core functions and definition of the WFdata class which holds the 
time-frequency data and the associated axes

NOTE: currently only 2D data arrays are supported
"""

__all__ = ['WFdata']


import os
import numpy as np


#=== Globals ===
global _WF_DESC_COMPULSORY

_WF_DESC_COMPULSORY = ['DataUnit', 'TimeUnit', 'ChanUnit',
                        'DataFrame', 'TimeFrame', 'ChanFrame']

global _SUPPORTED_DATA_FRAMES
global _SUPPORTED_TIME_FRAMES
global _SUPPORTED_CAN_FRAMES

_SUPPORTED_DATA_FRAMES = ['G']
_SUPPORTED_TIME_FRAMES = ['MJD', 'UNIX']
_SUPPORTED_CAN_FRAMES = ['freq', 'vel']


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
                wf_desc = {'DataUnit' : 'amp',
                    'TimeUnit' : 's',
                    'ChanUnit' : 'Hz',
                    'DataFrame' : 'G',
                    'TimeFrame' : 'UNIX',
                    'ChanFrame' :'freq'}):
        """

        """

        self.data_array = data_array
        self.time_array = time_array
        self.chan_array = chan_array
        self.wf_desc = wf_desc

        for WF_DESC in _WF_DESC_COMPULSORY:
            if WF_DESC not in set(self.wf_desc.keys()):
                raise VAlueError('The key: {0:s} is not defined in the \
    WF description (wf_descr)!'.format(WF_DESC))


        if self.wf_desc['DataFrame'] not in _SUPPORTED_DATA_FRAMES:
            raise TypeError('The data frame: {0:s} is not supported!'.format(data_frame))

        if self.wf_desc['TimeFrame'] not in _SUPPORTED_TIME_FRAMES:
            raise TypeError('The data frame: {0:s} is not supported!'.format(time_frame))

        if self.wf_desc['ChanFrame'] not in _SUPPORTED_CAN_FRAMES:
            raise TypeError('The data frame: {0:s} is not supported!'.format(chan_frame))


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

        #Save to a numpy binary file
        if compressed:
            np.savez_compressed(outfile,
                    data_array = self.data_array,
                    time_array = self.time_array,
                    chan_array = self.chan_array,
                    wf_desc = wf_desc_array)
        else:            
            np.savez(outfile,
                    data_array = self.data_array,
                    time_array = self.time_array,
                    chan_array = self.chan_array,
                    wf_desc = wf_desc_array)

#=== Functions II ===
def load_WFdata(input_path, input_fname):
    """

    """
    #Append WFdata_name
    input_fname += '.npz'

    input_WFname = os.path.join(input_path, input_fname)

    with np.load(input_WFname, mmap_mode='r') as raw_WFdata:
        #Convert numpy array to dict
        wf_desc_dict = {}
        raw_wf_desc_array = raw_WFdata['wf_desc']
        for k, v in zip(raw_wf_desc_array[:,0], raw_wf_desc_array[:,1]):
            wf_desc_dict[k] = v


        wf_obj = WFdata(data_array = raw_WFdata['data_array'],
                        time_array = raw_WFdata['time_array'],
                        chan_array = raw_WFdata['chan_array'],
                        wf_desc = wf_desc_dict)

        return wf_obj


#=== MAIN ===

working_dir = './'
test_file_name = 'test'

datarr = np.array([[1.+1.j,0.+0.j],[0.+0.j,0.+0.j],[0.+0.j,1.+0.j]], dtype=complex)
tarr = np.array([0,1,2])
farr = np.array([0,1])

print(type(datarr[0,0]))

a = WFdata(data_array = datarr,
            time_array = tarr,
            chan_array = farr)


a.save_WFdata(working_dir,test_file_name,overwrite=True)

b = load_WFdata(working_dir, test_file_name)

print(type(b.data_array[0,0]))

print(b.data_array)

exit()


print(a.wf_desc)

b = np.array(list(a.wf_desc.items()))

print(type(b),b)


c = {}

for k, v in zip(b[:,0], b[:,1]):
    c[k] = v

print(b[:,0], b[:,1])
print(zip(b[:,0], b[:,1]))

print(c)