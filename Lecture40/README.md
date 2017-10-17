# Calculate the differential production cross section of pp -> Z -> e+e-

This code computes the differential production cross section of pp -> Z -> e+e- using the [MSTW 2008 PDF set](https://mstwpdf.hepforge.org/code/code.html).
It also computes the total production cross section. This can be compared to the data collected by
the CMS experiment at 7 TeV [here](https://hepdata.net/record/ins921788). 



## Instructions for MSTW

- Download the [code](https://mstwpdf.hepforge.org/code/mstw2008code.tar.gz)
- Download the [PDFs](https://mstwpdf.hepforge.org/code/mstw2008grids.tar.gz)

Unpack:
```
tar -zxvf mstw2008code.tar.gz
tar -zxvf mstw2008grids.tar.gz
```

Compile the MSTW c++ file:
```
g++ -c -Wall -W mstwpdf.cc -o mstwpdf_cpp.o
```

Now compile our test example:
```
g++ mstw.cpp mstwpdf_cpp.o -o mstw
```
