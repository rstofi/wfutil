""" Miscellaneous utility functions useful for testing
"""

__all__ = ['echo_for_loop_counter', 'convert_list_to_string']

import os
import sys
import time

import wfutil as wf

#*******************************************************************************
#=== Functions ===
def echo_for_loop_counter(start,end,count,loop_string='Loop state'):
    """Simple routine to print the status of a for loop to the stdout and update
    it dynamically.
    
	It is useful for testing and running kinda long code on a screen session

	The code should be called within the loop, but it only works with for loops
	as one needs to know the start and end values of the loop.

	NOTE if it is logged every loop cycle will print a line to the logfile (!)

	NOTE this code is not compatible with the logging module
	
	NOTE no other print statements should be in the for loop (!)

    Parameters:
    ===========
    start: int
    	The start value of the loop
    end: int
    	The end value of the loop

	count: int
		The current vale of the loop

	loop string: str
		A string to name the loop. The counter is printed in the following format:

		`loop string`: [===.......] x% y/zzzz

	Returns:
	========
		Prints the loop status to stderr
    """

    #=== Compute percentage and object counter ===
    #Cause int() always round downwards == np.floor()
    loop_range = int(end - start)
    floor_percentage = int(100 * (count - start) / loop_range)

    #=== Build percentage bar string ===
    floor_10_percentage_bar = int(floor_percentage * 0.1)

    bar = ''.join(['=' for bar in range(floor_10_percentage_bar)])
    bar_left = ''.join(['.' for bar in range(10 - floor_10_percentage_bar)])
    floor_percentage_bar = bar + bar_left

    #=== Display ===
    print('{4:s}: [{0:s}] {1:d}% {2:d}/{3:d}'.format(
        floor_percentage_bar,floor_percentage,count - start,loop_range,
        loop_string),
        end='\r')

    if count == end:
    	print('{4:s}: [{0:s}] {1:d}% {2:d}/{3:d}'.format(
        floor_percentage_bar,floor_percentage,count - start,loop_range,
        loop_string))

def convert_list_to_string(list_to_convert):
    """Litarally convert a list to a string with format:

    [element_1,...,element_N]

    Parameters:
    ===========
    list_to_convert: list of anything
        The list to be converted
    
    Returns:
    ========
    literal_list_string: str
        The string created from the list
    """
    literal_list_string = '[' 
    listToStr = ','.join(map(str, list_to_convert))

    literal_list_string += listToStr + ']'

    return literal_list_string

#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass