"""Functions to interact with Measurement Sets (MS).

The main function of these functions is to extract data in the time-frequency
space from MS.
"""

__all__ = ['create_MS_object', 'close_MS_object', 'get_MS_subtable_path',
             'get_chann_array_from_MS', 'get_antenna_name_list_from_MS',
             'get_antname_and_ID_dict_from_MS', 'get_pol_array_from_MS',
             'get_fieldname_and_ID_list_dict_from_MS']

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


#Define the indices for the correlator types in MSv2.
#
#See e.g.: https://github.com/SKA-ScienceDataProcessor/algorithm-reference-library/
#blob/1b2c8d6079249202864abf8c60cdea40f0f123cb/processing_components/visibility/base.py#L644
#
# 1-4: Stokes polarization frame
# 5-8: Circular polarization frame
# 9-12: Linear polarization frame
#
_MSV2_CORR_TYPE_DICT = {1:'I',
                        2:'Q',
                        3:'U',
                        4:'V',
                        5:'RR',
                        6:'RL',
                        7:'LR',
                        8:'LL',
                        9:'XX',
                        10:'XY',
                        11:'YX',
                        12:'YY'}

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
    subtable_path = [subtables_path for subtables_path in MS.getsubtables() if '/' + subtable_name in subtables_path][0]

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
    chan_table_path = wf.msutil.get_MS_subtable_path(mspath,'SPECTRAL_WINDOW')

    chan_table = wf.msutil.create_MS_object(chan_table_path)

    chan_array = copy.deepcopy(chan_table.getcol('CHAN_FREQ')[0])

    wf.msutil.close_MS_object(chan_table)

    return chan_array

def get_pol_array_from_MS(mspath):
    """he MS is structured such a way, that the polarisation information is the
    same for all time chunks (this is an assumption), and so we need to get it
    independently from the time sub-selection of the data.

    This subroutine gets the polarisation array from the MS.

    TO DO: get the pol frame and units as well

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    pol_array: <numpy.Ndarray>
        Numpy array containing the polarisation values as strings

    """
    pol_table_path = wf.msutil.get_MS_subtable_path(mspath,'POLARIZATION')

    pol_table = wf.msutil.create_MS_object(pol_table_path)

    #Use the correlator type dictionary based on MSv2

    #NOTE if multiple fields with different correlator products are used the code only looks at the first row!
    corr_type = pol_table.getcol('CORR_TYPE')[0]
    N_pol = pol_table.getcol('NUM_CORR')[0]

    if np.size(corr_type) != N_pol:
        raise ValueError('Inconsistent number of polariastions and defined correlator types!')

    #Generate the human readable polarisation string from the IDs and dictionary
    pol_array = []
    for p in range(0,N_pol):
        pol_array.append(_MSV2_CORR_TYPE_DICT[corr_type[p]])

    wf.msutil.close_MS_object(pol_table)

    return np.array(pol_array)

def get_antname_and_ID_dict_from_MS(mspath):
    """Generate the antenna name -- ID pairs from an MS and return it as a dictionary.

    This is a key function as antenna naming and ID can be inconsistent between MS'

    e.g. ant named m004 can have an ID of 2 (e.g. whan ant m000 and m002 are missing
    from the MS)

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    antname_ID_dict: dict
        Dictionatry containing the antenna names and the corresponding antenna ID's

    """
    anttable_path = wf.msutil.get_MS_subtable_path(mspath,'ANTENNA')

    #Open `ANTENNA` table and read the list of antennas
    anttable = wf.msutil.create_MS_object(anttable_path)
    #ant_name_list = copy.deepcopy(np.unique(anttable.getcol('NAME')))

    #Set up the empty list
    antname_ID_dict = {}

    #The row number in the ANTENNA table corresponds to the antenna ID
    #See: https://casaguides.nrao.edu/index.php?title=Measurement_Set_Contents

    #Loop trough the rows in the ANTENNA table and build the dict
    for i in anttable.rownumbers():
        antname_ID_dict[anttable.getcol('NAME')[i]] = i

    wf.msutil.close_MS_object(anttable)

    return antname_ID_dict


def get_antenna_name_list_from_MS(mspath):
    """Return the list of antenna names from the given MS this is just a
    convienience function
    
    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    ant_name_list: list of strings
        List containing the antenna names

    """
    antname_ID_dict = get_antname_and_ID_dict_from_MS(mspath)

    return list(antname_ID_dict.keys())


def get_fieldname_and_ID_list_dict_from_MS(mspath,scan_ID=False):
    """Generate the field name -- ID list pairs or with `scan_ID` set to True,
    the scan IDs returned instead of the field IDs, from an MS as a dictionary.

    This is a key function as each field/scan could be queued by it's ID

    Thius is a core function in inspecting MS for field selection

    NOTE: a field can have different ID's attached i.e. the OTZFDUMMY field in
            the CNSS data set. The same applies to scans obviously

    NOTE: I am not 100% if this is how the MS works in terms of field ID's, but
            based on tests, seems legit

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    fieldname_ID_dict: dict
        Dictionatry containing the field names and the corresponding field ID's

    """
    fieldtable_path = wf.msutil.get_MS_subtable_path(mspath,'FIELD')

    #Open `FIELD` table and read the list of antennas
    fieldtable = wf.msutil.create_MS_object(fieldtable_path)

    #Set up the empty list
    fieldname_ID_dict = {}

    #The row number in the ANTENNA table corresponds to the field ID
    #See: https://casaguides.nrao.edu/index.php?title=Measurement_Set_Contents

    #I *assume* the same thing is true for the FIELD ID

    #Loop trough the rows in the ANTENNA table and build the dict
    for i in fieldtable.rownumbers():
        
        field_name = fieldtable.getcol('NAME')[i]
        
        if field_name not in list(fieldname_ID_dict.keys()):
            fieldname_ID_dict[field_name] = [i]
        else:
            fieldname_ID_dict[field_name].append(i)
            

    wf.msutil.close_MS_object(fieldtable)

    #The fieldname ID dict has to be creted to generate the scanname_ID_disct:
    if scan_ID:

        #NOTE this is a really slow sub-routine!

        MS = wf.msutil.create_MS_object(mspath)

        scanname_ID_disct = {}

        fieldnames = list(fieldname_ID_dict.keys())

        for field_name in fieldnames:

            field_selection_string = wf.miscutil.convert_list_to_string(list(fieldname_ID_dict[field_name]))

            scan_qtable = MS.query(query='FIELD_ID IN {0:s}'.format(
                        field_selection_string),
                        columns='SCAN_NUMBER')

            scanname_ID_disct[field_name] = list(np.unique(scan_qtable.getcol('SCAN_NUMBER')))

        wf.msutil.close_MS_object(MS)
    
        return scanname_ID_disct

    else:
        return fieldname_ID_dict

def get_field_name_list_from_MS(mspath):
    """Return the list of field names from the given MS this is just a
    convienience function
    
    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    Returns
    =======
    field_name_list: list of strings
        List containing the antenna names

    """
    fieldname_ID_dict = get_fieldname_and_ID_dict_from_MS(mspath)

    return list(fieldname_ID_dict.keys())

def get_baseline_data(mspath, ant1, ant2, qcolumn,
                    apply_flag=False,
                    field_ID_list=None,
                    scan_ID_list=None):
    """Core function for retrieveing time-frequency data from an MS. 

    The code uses a TAQL query to create a subtable selecting the given columns
    from the MAIN table and convert the data to a `<numpy.Ndarray>` obect.

    I assume using TAQL is the fastest solution available from the python casacore
    API...

    The baseline is defined by the antenna IDs that can be computed from the names
    using higher-level functions.

    A subset of scans/fields can be selected, however if one wants to select e.g.
    all fields from a scan and onla a single field from another scan, they must 
    provide both scan IDs and ALL field IDs, including from the one scan all fields
    are selected.

    That is, the scan/field selection is really flexible but has to be specified
    with great precision!

    NOTE: it is faster and better practice to select based solely on field ID,
        the option to select based on  scan ID is for testing purposes mostly...

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

    apply_flag: bool, opt
        If True, the flag columns are also read and the flag tables are applied to
        the data column.

        NOTE currently only works with `qcolumn` == 'DATA'

        if `qcolumn` == 'FLAG' this option is automatically disabled
    
    field_ID_list: list of int, optional
        A list of field IDs to select. If None, all fields are selected.

    scan_ID_list: list of int, optional
        A list of scan IDs to select. If None, all scans are selected.

    Returns
    =======
    wf_data: <numpy.Ndarray>
        The waterfall data matrix for the selected baseline

    tarray: <numpy.Ndarray>
        The time axis of the waterfall matrix for the selected baseline

    """
    #Check if the queried column is supported
    if qcolumn not in _SUPPORTED_MS_COLUMNS:
        raise ValueError('Query the column {0:s} is not supported!'.format(qcolumn))

    #Do not apply flags when the flag table is retriewed
    if qcolumn == 'FLAG' and apply_flag == True:
        apply_flag = False

    #The MS has the lower index antenna as ant1 (always I think), so swap the antennas,
    #if ant1 have larger index than ant2!
    if ant1 > ant2:
        antx = ant1
        ant1 = ant2
        ant2 = antx
        del antx

    MS = wf.msutil.create_MS_object(mspath)

    #Query the data
    def query_table(q_ant1, q_ant2, q_qcolumn, scan_ID_list,field_ID_list):
        """I am not sure if this is a good coding practice, but I want to code the
        table query only once and not duplicate it...
        
        Also I am using sligtly different variable names useful for debugging...
        """
        if scan_ID_list == None and field_ID_list == None:
            q_qtable = MS.query(query='ANTENNA1 == {0:d} AND ANTENNA2 == {1:d}'.format(
                         q_ant1,q_ant2), columns='{0:s},TIME'.format(q_qcolumn))

        else:
            #Get selection strings
            if field_ID_list != None:

                field_selection_string = ' AND FIELD_ID IN '
                field_selection_string += wf.miscutil.convert_list_to_string(field_ID_list)

            if scan_ID_list != None:

                scan_selection_string = ' AND SCAN_NUMBER IN '
                scan_selection_string += wf.miscutil.convert_list_to_string(scan_ID_list)

            #Actually query the data
            if scan_ID_list != None and field_ID_list == None:
                q_qtable = MS.query(query='ANTENNA1 == {0:d} AND ANTENNA2 == {1:d}{2:s}'.format(
                                    q_ant1,q_ant2,scan_selection_string),
                            columns='{0:s},TIME'.format(q_qcolumn))

            elif scan_ID_list == None and field_ID_list != None:
                q_qtable = MS.query(query='ANTENNA1 == {0:d} AND ANTENNA2 == {1:d}{2:s}'.format(
                                    q_ant1,q_ant2,field_selection_string),
                            columns='{0:s},TIME'.format(q_qcolumn))

            else:
                q_qtable = MS.query(query='ANTENNA1 == {0:d} AND ANTENNA2 == {1:d}{2:s}{3:s}'.format(
                                    q_ant1,q_ant2,field_selection_string,scan_selection_string),
                            columns='{0:s},TIME'.format(q_qcolumn))

        #Copy the results inot a data array and a time array
        q_tarray = copy.deepcopy(q_qtable.getcol('TIME'))
        q_qarray = copy.deepcopy(q_qtable.getcol(q_qcolumn))

        #Close the selection result table
        q_qtable.close()

        return q_qarray, q_tarray

    #Calling it for the input parameters
    qarray, tarray = query_table(q_ant1=ant1, q_ant2=ant2,
                                q_qcolumn=qcolumn,
                                field_ID_list=field_ID_list,
                                scan_ID_list=scan_ID_list)

    #Get the flags if needed
    if apply_flag == True:
        #Getting the flag from the table
        m_qarray, _ = query_table(q_ant1=ant1, q_ant2=ant2,
                                q_qcolumn='FLAG',
                                field_ID_list=field_ID_list,
                                scan_ID_list=scan_ID_list)
        
        #Applying the flag (NOTE that the axis are not WFdata ordered!)
        qarray = np.ma.masked_array(qarray, m_qarray)
    
    #No column including spectral info only TIME info...
    #TO DO: somehow check for this
    #print(MS.colnames())

    #NOTE this is only a problem if the time and channel shape is the same...

    #chan_array = wf.msutil.get_chan_array_from_MS(mspath)
    #c_shape = np.size(chan_array)

    #Check the polarisation array
    pol_array = wf.msutil.get_pol_array_from_MS(mspath)
    p_shape = np.size(pol_array)

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

def get_baseline_flag_fraction_data(mspath, sant=None, 
                                    echo_counter=False,
                                    continue_one_baseline_error=False,
                                    field_ID_list=None,
                                    scan_ID_list=None):
    """A wrapper around `get_baseline_data()`. In particular, this routine generates
    waterfall plots for the fraction of baselines flagged (0 = none, 1 = all)

    This script does not take autocorrelations, but only baselines into account.

    It depends on the inherent property of the MS being the same dimension for
    all baselines (i.e. `get_baseline_data` returns the same-sized quarr and the
    tarr are the same for any baseline). NOTE that the script does not check if the
    time stams are the same but onyl for the array shapes

    NOTE: this routine is not too fast...

    To derive the flag solutions for all baselines including a selected antenna,
    specify the antenna name via the 'sant' variable. This is the way to generate
    crosscorr flag tables that acan be applied to single-dish data!

    Parameters
    ==========
    mspath: str
        The input MS path or a ``casacore.tables.table.table`` object

    sant: string, opt
        If an antenna name is given only the flags from the baselines including
        the selected antenna will be returned. 

    echo_counter: bool, opt
        If True a counter of baselines processed is printed to the stdout

    continue_one_baseline_error: bool, opt
        If True, the code terminates even if not all baselines are processed

    field_ID_list: list of int, optional
        A list of field IDs to select. If None, all fields are selected.

    scan_ID_list: list of int, optional
        A list of scan IDs to select. If None, all scans are selected.

    Returns
    =======
    frac_qarr: <nump.Ndarray>
        The waterfall data matrix containing the fraction of flagged baselines

    frac_tarr: <numpy.Ndarray>
        The time axis for `frac_qarr`

    """
    MS = wf.msutil.create_MS_object(mspath)

    #Get the list of antennas and indices
    antname_ID_dict = get_antname_and_ID_dict_from_MS(MS)

    #Get unique antenna indices list
    unique_ant_name_list = list(antname_ID_dict.keys())
    unique_ant_ID_list = np.sort(list(antname_ID_dict.values())) #Sort the IDs for covieneient processing

    #If sant is given, chack if exist in the dict
    if sant != None:
        if sant not in unique_ant_name_list:
            raise ValueError('Provided sant is not a valid antenna name!')

        #Get sant ID
        sant_ID = antname_ID_dict[sant]
    else:
        sant_ID = None

    #Get the list of antenna ID's for ant1 and ant2 (can't get from the dict)
    ant1_ID_list = np.unique(MS.getcol('ANTENNA1'))
    ant2_ID_list = np.unique(MS.getcol('ANTENNA2'))

    #Check if this is the same as the union of ant1 and ant2 lists
    if np.all(unique_ant_ID_list) != np.all(np.unique([ant1_ID_list,ant2_ID_list])):
        raise ValueError('Antenna name missmach in abselines and ANTENNA table!')

    #Nuber of expected baselines => used only for loop counter
    if sant != None:
        N_b = len(unique_ant_name_list) - 1
    else:
        N_b = int((len(unique_ant_name_list) * (len(unique_ant_name_list) - 1)) / 2)

    #Get initial data sizes (assuming that this baseline exist!)
    qarr, tarr = get_baseline_data(MS, ant1=ant1_ID_list[0], ant2=ant1_ID_list[1],
                                    qcolumn='FLAG',
                                    field_ID_list=field_ID_list,
                                    scan_ID_list=scan_ID_list)

    if qarr.dtype != bool:
        raise TypeError('Error when reading in boolean mask from MS!')

    #Now create master qarr to average onto
    frac_qarr = np.zeros(np.shape(qarr))
    frac_tarr = copy.deepcopy(tarr)

    #I use a brute-force check for all baselines regardless how they organised
    #Thus I loop through all possible baseline configurations

    N_b_count = 0
    for i in unique_ant_ID_list:
        for j in unique_ant_ID_list:
            if i >= j:
                #Basically ignore possible autocorrelations 
                #Plus the already tried combinations
                continue
            else:
                #Check if baseline is valid
                if i in ant1_ID_list and j in ant2_ID_list:

                    #If sant is given check if it is in the baseline being processed
                    if sant != None:
                        #If not go to next baseline
                        if sant_ID not in [i,j]:
                            continue
                        else:
                            pass
                    else:
                        pass

                    #print(i,j)
                    qarr, tarr = get_baseline_data(MS,
                                ant1=i,
                                ant2=j,
                                qcolumn='FLAG',
                                field_ID_list=field_ID_list,
                                scan_ID_list=scan_ID_list)

                    #Check array sizes but not the individual time stamps
                    if np.shape(frac_qarr) != np.shape(qarr) or np.shape(frac_tarr) !=  np.shape(tarr):
                        raise ValueError('Retrieved baseline data shape mismatch!')

                    #qarr should be bool, so I convert it to zeros and ones and
                    #add it to the master array
                    frac_qarr = np.add(frac_qarr,qarr.astype(int))

                    N_b_count += 1

                    if echo_counter:
                        wf.miscutil.echo_for_loop_counter(0,N_b,N_b_count, 'Baselines processed')

                else:
                    pass 

    #Check if all baselines found and processed
    if not continue_one_baseline_error:
        if N_b != N_b_count:
            raise ValueError('Some baselines either missed by wfutil or missing from the MS!')

    #Now normalise with the number of baselines
    frac_qarr = np.divide(frac_qarr, N_b)

    return frac_qarr, tarr

#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass
