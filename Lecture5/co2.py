import math
import matplotlib.pyplot as plt

def sine_transform(data):
    """Return Fourier sine transform of a real data vector"""
    N = len(data)
    transform = [ 0 ] * N
    for k in range(N):
        for j in range(N):
            angle = math.pi * k * j / N
            transform[k] += data[j] * math.sin(angle)
    return transform

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


transform_raw = sine_transform(data)
abs_transform = [ abs(x) for x in transform_raw]

freqs = [ float(i) for i in xrange(len(transform_raw))]


plt.subplot(2, 1, 1)
plt.plot( dates, data )

ax = plt.subplot(2, 1, 2)
plt.plot( freqs, abs_transform )
ax.set_yscale('log')

plt.show()
