"""Collection of fuctions operation on the time array
"""

__all__ = ['convert_MJD_to_UNIX']


from astropy.time import Time

#=== Functions ===
def convert_MJD_to_UNIX(time_array):
    """Convert the MJD (Modified Julian Date) times to UNIX time

    MJD is the default time format in some MS, but it is more convienient to
    work with UNIX time

    Parameters:
    ===========

    time_array: <numpy.ndarray>
        The numpy array containing the time values in MJD format

    Returns:
    ========
    unix_time_array: <numpy.ndarray>
        The numpy array containing the time values in UNIX format
    """
    #Conversion
    unix_time_array = Time(time_array / 86400, format='mjd')
    unix_time_array.format = 'unix'

    #Do not return an Astropy Time object but a numpy array
    return unix_time_array.value

#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass
