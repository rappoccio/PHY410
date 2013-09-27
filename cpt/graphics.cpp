#include "graphics.hpp"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
using namespace std;

namespace cpt {

// Gnuplot plots and movies

    void Gnuplot::animate(int steps)
    {
#       ifdef _WIN32
            string gnuplot = "env pgnuplot.exe ";
            string terminal = "windows";
            FILE *gnupipe = _popen(gnuplot.c_str(), "w");
#       else
            string gnuplot = "/usr/bin/env gnuplot ";
            string terminal = "x11";
            FILE *gnupipe = popen(gnuplot.c_str(), "w");
#       endif
        if (!gnupipe) {
            time_t start_time = clock();
            while (!gnupipe &&
                   (clock() - start_time) / double(CLOCKS_PER_SEC) < 5.0)
                ;
        }
        if (!gnupipe) {
            cerr << " Gnuplot::animate: cannot open pipe" << endl;
            exit(1);
        }
        double dt_frame = 1 / frame_rate;
        for (int step = 0; step < steps; step++) {
            time_t start_time = clock();
            step_callback();
            ostringstream os;
            for (int i = 0; i < commands.size(); i++) {
                string cmd = commands[i];
                string::size_type pos = cmd.find('\n');
                if (pos == string::npos)
                    cmd += '\n';
                os << cmd;
            }
            fprintf(gnupipe, "%s", os.str().c_str());
            fflush(gnupipe);
            while (true) {
                double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
                if (secs > dt_frame)
                    break;
            }
        }
        fclose(gnupipe);
    }

// OpenGL and GLUT

    // http://msdn.microsoft.com/en-us/library/ms680621%28VS.85%29.aspx
    // SetErrorMode(GetErrorMode () | SEM_NOGPFAULTERRORBOX);

} /* namespace cpt */
