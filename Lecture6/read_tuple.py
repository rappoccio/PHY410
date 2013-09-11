import math

def read_tuple(filename) :
    """This is a simple function that can be used to
    read in a "tuple" of data from a text file, which
    is separated by newline characters.
    It returns an array of values in python.

    Args : filename : name of the file

    Returns : an array of the values, or None if error. 
    
    """
    try: 
        file = open( filename, 'r')
        lines = file.readlines()
        x = []
        for line in lines :
            x.append( float(line) )
        return x
    except:
        print 'Error in read_tuple. Returning None'
        return None

