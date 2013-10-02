from cpt import *




# data downloaded from http:#neic.usgs.gov/neis/epic/epic_global.html
co2_data = "co2_mm_mlo.txt"

print " CO2 data "

x_values = []
y_values = []

# data downloaded from ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt
print ' C02 Data from Mauna Loa'
data_file_name = 'co2_mm_mlo.txt'
file = open(data_file_name, 'r')
lines = file.readlines()
file.close()
print ' read', len(lines), 'lines from', data_file_name

window = False

xvec = []
yvec = []

for line in lines :
    if line[0] != '#' :
        try:
            words = line.split()
            xval = float( words[2] )
            yval = float( words[4] )
            yvec.append( yval )
            xvec.append( xval )
        except ValueError :
            print 'bad data:',line

N = len(xvec)
sigma = [0.001] * N


POLYNOMIAL_ORDER = 3
# fit histogram to straight line
a_fit = Matrix(POLYNOMIAL_ORDER,POLYNOMIAL_ORDER)
sig_a = Matrix(POLYNOMIAL_ORDER,POLYNOMIAL_ORDER)
yy = [0.0] * N
chi2 = 0.0

[a_fit,sig_a,yy,chi2] = polyfit(xvec, yvec, sigma, POLYNOMIAL_ORDER)
print a_fit
print sig_a
print chi2



    
co2_polyfit_output = open("co2_polyfit_output_python.txt", 'w')

for i in xrange(N) :    
    s = "{0:6.2f} {1:6.2f} {2:6.2f}\n".format( xvec[i], yvec[i], yy[i] )
    co2_polyfit_output.write(s)


co2_polyfit_output.close()

