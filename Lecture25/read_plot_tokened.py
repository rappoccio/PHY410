import math

def read_plot_tokened(filename) :
    """This is a simple function that can be used to
    read in a plot of SETS of data from a text file, which
    is separated by newline characters. Each set itself
    is separated by a double newline character. 
    It returns the arrays as x,y

    Args : filename : name of the file

    Returns : x[] : x values
              y[] : y values
    
    """
    try: 
        file = open( filename, 'r')
        lines = file.readlines()
        x = [ [] ]
        y = [ [] ]
        itok = 0
        for line in lines :
            if line == '\n' :
                x.append( [] )
                y.append( [] )
                itok += 1
            else :                
                words = line.split()
                ix, iy = [float(s) for s in words]
                x[itok].append( ix )
                y[itok].append( iy )
        # This always botches the last one so remove it : 
        return x[0:len(x)-1],y[0:len(y)-1]  
    except:
        print 'Error in read_plot_tokened. Returning None'
        return None

