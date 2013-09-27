Computational Physics Toolkit
=============================

Get and Install Current Version:

    Download: http://www.physics.buffalo.edu/phy410-505/tools/cpt.zip

    In PowerShell or Bash, cd to download directory

        ls  cpt.zip
        unzip  cpt.zip
        cd  cpt/
        make
        cd  ..

Compile and link your application code in file myapp.cpp

    Assuming  myapp.cpp  and  cpt.zip  are in the same directory

        c++  -Icpt  myapp.cpp  -Lcpt  -lcpt

    If  myapp.cpp  and  cpt.zip  are in different directories on Windows
        C:\cptfolder\cpt.zip
        C:\cptfolder\cpt\matrix.hpp
        C:\cptfolder\cpt\libcpt.a

        c++  -IC:\cptfolder\cpt  myapp.cpp  -LC:\cptfolder\cpt  -lcpt
