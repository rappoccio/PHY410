/*
   Hopfield Neural Network Model
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

const int N = 10;                   // number of neurons
const int Mmax = 26;                // maximum number of memories
int M = 1;                          // number of memories
double P = 30;                      // error percentage
int t = 0;                          // time step

int J[N][N][N][N];                  // connection strength
int L[Mmax][N][N];                  //
int B[N][N];                        // bias term
int H[N];                           // Hamming distances
int X[N], Y[N];                     // visited neurons

#include "letters.hpp"              // 10x10 bitmaps of 26 letters

bool regularSweep;                  // if false do random sweep
bool bias;                          // use bias if true
int iNext = 0, jNext = -1;          // sweep indicies

void getInput() {
    cout << "Enter number of letters to memorize: ";
    cin >> M;
    cout << "Enter percentage errors: ";
    cin >> P;
}

void initialize() {
    int i, j, k, l, m;

    if (M < 1) M = 1;
    if (M > Mmax) M = Mmax;

    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            for (k = 0; k < N; ++k)
                for (l = 0; l < N; ++l)
                    J[i][j][k][l] = 0;

    for (m = 0; m < M; ++m)
        for (i = 0; i < N; ++i)
            for (j = 0; j < N; ++j)
                for (k = 0; k < N; ++k)
                    for (l = 0; l < N; ++l)
                        J[i][j][k][l] += Letter[m][i][j] * Letter[m][k][l];

    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            J[i][j][i][j] = B[i][j] = 0;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                for (l = 0; l < N; l++)
                    B[i][j] += J[i][j][k][l];

    for (m = 0; m < M; ++m)
        for (i = 0; i < N; ++i)
            for (j = 0; j < N; ++j)
                if (std_rand() > P / 100.0)
                    L[m][i][j] = Letter[m][i][j];
                else L[m][i][j] = -Letter[m][i][j];

    t = 0;
    iNext = 0;
    jNext = -1;
}

double sumH;                        // Hamming sum
void computeHammingDistances() {
    sumH = 0;
    for (int m = 0; m < M; ++m) {
        H[m] = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int d = Letter[m][i][j] - L[m][i][j];
                H[m] += d * d;
            }
        }
        sumH += H[m] / (double) (4 * N * N);
    }
}

double framesPerSec = 30;           // control animation rate if too fast

void takeStep() {
    clock_t start_time = clock();   // get CPU clock time

    if (++jNext == N) {
        jNext = 0;
        if (++iNext == N)
            iNext = 0;
    }

    for (int m = 0; m < M; ++m) {
        int i = 0;
        int j = 0;
        if (regularSweep) {
            i = iNext;
            j = jNext;
        } else {
            i = int(std_rand() * N);
            j = int(std_rand() * N);
        }
        X[m] = j;
        Y[m] = i;
        double E_flip = 0;
        for (int k = 0; k < N; ++k)
            for (int l = 0; l < N; ++l)
                E_flip += J[i][j][k][l] * L[m][k][l];
        if (bias)
            E_flip -= 0.5 * B[i][j];
        E_flip *= 2 * L[m][i][j] / double(M);
        if (E_flip < 0)
            L[m][i][j] = -L[m][i][j];
    }

    computeHammingDistances();
    ++t;

    // delay if animation is too fast
    while ((double(clock()) - start_time)/CLOCKS_PER_SEC < 1.0 / framesPerSec)
        ;

    glutPostRedisplay();
}

void drawText(const string& str, int x, int y) {
     glRasterPos2i(x, y);
     int len = str.find('\0');
     for (int i = 0; i < len; i++)
         glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]);
}

double Wdx = 10, Wdy = 11;           // width and height of world

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    int Lx = int(sqrt(double(M)));
    if (Lx * Lx != M)
        Lx += 1;
    int Ly = Lx;
    double letterSize = Wdx / Lx;
    double pixelSize = letterSize / N;
    double radius = letterSize / 4 / N;
    int m = 0;
    int x = 0;
    int y = 0;
    while (m < M) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                if (L[m][N - j - 1][i] == 1) {
                    if (H[m] == 0)
                        glColor3ub(0, 0, 0);      // black
                    else
                        glColor3ub(255, 0, 0);    // red
                } else
                    glColor3ub(255, 255, 0);      // yellow
                double x0 = x * letterSize + i * pixelSize;
                double y0 = y * letterSize + j * pixelSize;
                double dd = 0.9 * pixelSize;
                glRectd(x0, y0, x0 + dd, y0 + dd);
            }
        }

        if (++x == Lx) {
            ++y;
            x = 0;
        }
        ++m;
    }
    ostringstream os;
    os << "Timestep = ";
    glColor3ub(0, 0, 0);
    drawText(os.str(), 1, 1);

    glutSwapBuffers();
}

int VPx0 = 0;            // left pixel of viewport on window
int VPdx = 400;          // width of viewport on window in pixels
int VPy0 = 0;            // bottom pixel of viewport on window
int VPdy = 410;          // height of viewport on window in pixels

void reshape(int w, int h) {
    double aspect = (w * Wdy) / (h * Wdx);
    VPx0 = VPy0 = 0;
    VPdx = w;
    VPdy = h;
    if (aspect > 1) {
        VPx0 = int(w * (1 - 1 / aspect) / 2.0);
        VPdx = int(w / aspect);
    } else {
        VPy0 = int(h * (1 - aspect) / 2.0);
        VPdy = int(h * aspect);
    }
    glViewport(VPx0, VPy0, VPdx, VPdy);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, Wdx, 0, Wdy);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

bool running;

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
         if (running) {
            glutIdleFunc(NULL);
            running = false;
        } else {
            glutIdleFunc(takeStep);
            running = true;
        }
    }
}

void reset(int menuItem) {
    switch (menuItem) {
    case 1:
        if (running) {
            glutIdleFunc(NULL);
            running = false;
        }
        initialize();
        glutPostRedisplay();
        break;
    default:
        break;
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    getInput();
    initialize();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(VPdx, VPdy);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Hopfield Neural Network Model");
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(NULL);
    glutMouseFunc(mouse);
    glutCreateMenu(reset);              // create a popup menu
    glutAddMenuEntry("Reset", 1);
    glutAttachMenu(GLUT_RIGHT_BUTTON);  // attach it to right mouse button
    glutMainLoop();
}
