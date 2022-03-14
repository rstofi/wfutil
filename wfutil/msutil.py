"""Functions to interact with Measurement Sets (MS).

The main function of these functions is to extract data in the time-frequency
space from MS.
"""

__all__ = ['create_MS_object', 'close_MS_object', 'get_MS_subtable_path',
             'get_chann_array_from_MS', 'get_antenna_name_list_from_MS']

import numpy as np
import logging
import copy

from casacore import tables as casatables

import wfutil as wf

#=== GLOBALS ===

global _ACK #Enabling messages of successful interaction with the MS e.g. successful opening of a table
global _SUPPORTED_MS_COLUMNS

_ACK = False
_SUPPORTED_MS_COLUMNS =  ['DATA', 'FLAG']

#*******************************************************************************
#=== Functions ===
def create_MS_object(mspath, readonly=True, **kwargs):
    """This function aims to speed up other bits of this module, 
    by returning a ``casacore.tables.table.table`` object.
    The trick is, that the ``mspath`` argument can be either a string i.e. the path
    to the MS which will be read in and returned, **or** it can be already an
    in-memory ``casacore.tables.table.table`` object.
    This might not be the best solution, but I hope overall a check in a lot of cases will
    speed up code, rather than reading in the same MS again-and again. So ideally, only
    one reading in happens for each MS and all inside this function!

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    readonly: bool, optional
        If True, the tables of the MS can be read only, but if set to False one can modify the MS

    Returns
    =======
    MS: ``casacore.tables.table.table`` object
        The in-memory Measurement Set
    """
    #create an empty MS in-memory to check the object type: the working solution
    MS_type_tmp = casatables.table('',casatables.maketabdesc([casatables.makescacoldesc('DATA',0)]),memorytable=True,ack=False)

    #if type(mspath) == 'casacore.tables.table.table': #This approach does not work
    if type(mspath) == type(MS_type_tmp):
        MS_type_tmp.close()
        return mspath
    else:
        #log.debug('Open MS: {0:s}'.format(str(mspath))) #We know it is a string in this case
        MS = casatables.table(mspath, ack=_ACK, readonly=readonly)
        MS_type_tmp.close()
        return MS

def close_MS_object(mspath):
    """This bit of code should be called at the end of whan working with an MS.
    Basically only closes the MS opened in the beginning. Aims to prevent memory
    leaks and just implementing good practices...

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
 

    """
    if type(mspath) == 'casacore.tables.table.table':
        MS = wf.msutil.create_MS_object(mspath)

        MS.close()
    else:
        pass

def get_MS_subtable_path(mspath, subtable_name):
    """Subroutine to generate absolute paths for MS subtables using the right
    syntax (i.e. ::SUBTABLE instead of /SUBTABLE)

    See: https://casacore.github.io/python-casacore/casacore_tables.html

    NOTE: this piece of code only works on UNIX systems!

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    subtable_name: str
        FULL name of the subtable.

    Returns
    =======
    subtable_path: str
        Absolute path to the subtable

    """
    MS = wf.msutil.create_MS_object(mspath)

    #List of subtables:
    #print(MS.keywordnames())

    #Select all substrings containing `subtable_name` and select the first result
    #NOTE: only one table sould exist named `/subtable_name`
    subtable_path = [anttables_path for anttables_path in MS.getsubtables() if '/' + subtable_name in anttables_path][0]

    #Get the index of the dash using reverse
    subtable_dash_index = subtable_path.rindex("/")

    subtable_path = subtable_path[:subtable_dash_index] + \
                    "::" + subtable_path[subtable_dash_index+1:]

    return subtable_path


def get_chann_array_from_MS(mspath):
    """The MS is structured such a way, that the channel information is the same
    for all time chunks, and so we need to get it independently from the time
    sub-selection of the data.

    This subroutine gets the channel array from the MS.

    TO DO: get the channel frame and units as well

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    chan_array: <numpy.Ndarray>
        Numpy array containing the channel values

    """
    chan_table_path = wf.msutil.get_MS_subtable_path(MS,'SPECTRAL_WINDOW')

    chan_table = wf.msutil.create_MS_object(chan_table_path)

    chan_array = copy.deepcopy(chan_table.getcol('CHAN_FREQ')[0])

    wf.msutil.close_MS_object(chan_table)

    return chan_array


def get_antenna_name_list_from_MS(mspath):
    """Return the list of antenna names from the given MS
    
    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    ant_name_list: list of strings
        List containing the antenna names

    """
    anttable_path = wf.msutil.get_MS_subtable_path(MS,'ANTENNA')

    #Open `ANTENNA` table and read the list of antennas
    anttable = wf.msutil.create_MS_object(anttable_path)
    ant_name_list = copy.deepcopy(np.unique(anttable.getcol('NAME')))
    wf.msutil.close_MS_object(anttable)

    return ant_name_list


def get_baseline_data(mspath, ant1, ant2, qcolumn):
    """Core function for retrieveing time-frequency data from an MS. 

    The code uses a TAQL query to create a subtable selecting the given columns
    from the MAIN table and convert the data to a `<numpy.Ndarray>` obect.

    I assume using TAQL is the fastest solution available from the python casacore
    API...

    The baseline is defined by the antenna IDs that can be computed from the names
    using higher-level functions.

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    ant1: int
        ID of the first antenna

    ant2: int
        ID of the second antenna

    qcolumn: str
        Name of the column to query, i.e. the column we select from the MAIN table.
        The currently supported columns are defined in the `_SUPPORTED_MS_COLUMNS`
        global variable

    Returns
    =======
    wf_data: <>

    """
    #Check if the queried column is supported
    if qcolumn not in _SUPPORTED_MS_COLUMNS:
        raise ValueError('Query the column {0:s} is not supported!'.format(qcolumn))

    MS = wf.msutil.create_MS_object(mspath)

    #Query the data
    qtable = MS.query(query='ANTENNA1 == {0:d} AND ANTENNA2 == {1:d}'.format(ant1,ant2), 
                      columns='{0:s},TIME'.format(qcolumn))

    #Copy the results inot a data array and a time array
    tarray = copy.deepcopy(qtable.getcol('TIME'))
    qarray = copy.deepcopy(qtable.getcol(qcolumn))


    #Close the selection result table
    qtable.close()

    #No column including spectral info only TIME info...
    #TO DO: somehow check for this
    #print(MS.colnames())

    #c_shape = ???

    #NOTE the current code only works for two+ ploarisations!
    #TO DO: check for exact polarisation info

    p_shape = np.amin(np.shape(qarray))

    #TIME array can be exactly extracted from the data

    t_shape = np.size(tarray)

    #Reshape data until it is in the correct wfutil shape:
    # [pol x time x chan]
    
    #This is a dumb, but quick routine
    while np.shape(qarray)[0] != p_shape or np.shape(qarray)[1] != t_shape:

        qarray = qarray.transpose() #Swap first and last
        qarray = np.moveaxis(qarray,0,1) #Swap first and second

    #print(np.shape(qarray))    

    return qarray, tarray


def loop_baselines(mspath):
    """A wrapper around `get_baseline_data()`. It can be used to get some statisics
    on the baseline axis.


    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======





    """
    pass


#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    #pass

    #exit()


#*******************************************************************************
    #=== Quick and dirty testing ===
    working_dir = '/home/krozgonyi/Desktop/playground/'

    phase_plots = True

    #=== First test MS (small one) ===

    mspath = working_dir + '1630519596-1934_638_d0_5-allflag.ms'
    #mspath = working_dir + '1630519596-1934_638-allflag.ms'

    MS = wf.msutil.create_MS_object(mspath)
    #print(type(MS))

    #ant_name_list = get_antenna_name_list_from_MS(MS)
    #print(ant_name_list)

    flag_matrix, time_array = get_baseline_data(MS, 0, 40,'DATA')
    chan_array = get_chann_array_from_MS(MS)

    wf.msutil.close_MS_object(MS)

    pol_array = np.array(['XX', 'YY'])

    #Create a WF object
    WFd1 = wf.core.WFdata(data_array = flag_matrix,
                    time_array = time_array,
                    chan_array = chan_array,
                    pol_array = pol_array)

    #Quick&dirty plot
    wf.visutil.quick_and_dirty_plot(WFd1, plot_phase=phase_plots)

    #=== Second test MS (large one) ===

    mspath = working_dir + '1630519596-1934_638-allflag.ms'

    MS = wf.msutil.create_MS_object(mspath)

    flag_matrix, time_array = get_baseline_data(MS, 0, 40,'DATA')
    chan_array = get_chann_array_from_MS(MS)

    wf.msutil.close_MS_object(MS)

    pol_array = np.array(['XX', 'YY'])

    #Create a WF object
    WFd2 = wf.core.WFdata(data_array = flag_matrix,
                    time_array = time_array,
                    chan_array = chan_array,
                    pol_array = pol_array)

    #Quick&dirty plot
    wf.visutil.quick_and_dirty_plot(WFd2, plot_phase=phase_plots)

    #=== Combine the two WFdata ====

    WFd_to_merge_list = [WFd1, WFd2]

    merged_WFd = wf.core.merge_WFdata(WFd_to_merge_list)

    #Quick&dirty plot
    wf.visutil.quick_and_dirty_plot(merged_WFd, plot_phase=phase_plots)