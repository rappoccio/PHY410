import math
import matplotlib.pyplot as plt

file_name = "co2_mm_mlo.txt"
file = open(file_name, "r")
lines = file.readlines()
file.close()

dates = []
data = []
for line in lines:
    if len(line) > 4:
        try:
            year = int(line[0:4])
            if year > 1957 and year < 2013:
                words = str.split(line)
                dates.append(words[2])
                ppm = float(words[4])
                if ppm > 0: data.append(ppm)
        except ValueError:
            pass

print " read", len(data), "values from", file_name


plt.plot( dates, data )

plt.show()
