"""The core functions and definition of the WFdata class which holds the 
time-frequency data and the associated axes
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

	for WF_DESC in _WF_DESC_COMPULSORY:
		if WF_DESC not in set(self._mapping.keys()):
			raise VAlueError('The key: {0:s} is not defined in the \
WF description (wf_descr)!'.format(WF_DESC))

	self._data_array = data_array
	self._time_array = time_array
	self._chan_array = chan_array
	self._wf_des = wf_desc


	if self._wf_desc['DataFrame'] not in _SUPPORTED_DATA_FRAMES:
    	raise TypeError('The data frame: {0:s} is not supported!'.format(data_frame))

	if self._wf_desc['TimeFrame'] not in _SUPPORTED_TIME_FRAMES:
    	raise TypeError('The data frame: {0:s} is not supported!'.format(time_frame))

	if self._wf_desc['ChanFrame'] not in _SUPPORTED_CAN_FRAMES:
    	raise TypeError('The data frame: {0:s} is not supported!'.format(chan_frame))


	def save_WFdata(self, output_path, WFdata_name, overwrite=True):
		"""

		"""
		#Append WFdata_name
		WFdata_name += '.npy'

		outpath = os.path.join(output_path, WFdata_name)

        if os.path.isdir(outpath) and overwrite == False: 
            raise TypeError('WFdata file already exist, and the \
overwrite parameters is set to False!')


        #Save to a numpy binary file

