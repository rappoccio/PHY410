#include <cmath>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

struct Site {               // object to represent lattice site
    int x;                  // x coordinate
    int y;                  // y coordinate

    // strict weak-ordering comparison operator required by set
    bool operator< (const Site& s) const {
        const int maxDiff = 9999;   // must be larger than x or y difference
        return (x - s.x) + maxDiff * (y - s.y) < 0;
    }
};

deque<Site> snake;          // double-headed reptile
set<Site> occupiedSites;    // set of occupied sites

bool occupied(Site s) {     // return true if s is occupied
    return occupiedSites.find(s) != occupiedSites.end();
}

void clear() {              // remove all sites
    snake.clear();
    occupiedSites.clear();
}

void addBack(Site s) {      // add s to back of reptile
    snake.push_back(s);
    occupiedSites.insert(s);
}

void addFront(Site s) {     // add s to back of reptile
    snake.push_front(s);
    occupiedSites.insert(s);
}

void removeBack() {         // remove back end of reptile
    occupiedSites.erase(snake.back());
    snake.pop_back();
}

void removeFront() {        // remove front end of reptile
    occupiedSites.erase(snake.front());
    snake.pop_front();
}

const int EAST = 0, NORTH = 1, WEST = 2, SOUTH = 3, DIRECTIONS = 4;

const int                   // initial configuration choices
    STAIR = 1,              // East-North random staircase
    COIL  = 2,              // tight coil East, North, West, West, South, ...
    LINE  = 3;              // straight line East

void createSnake(           // make a snake with
    int steps,              // this number of segments
    int config = LINE)      // and this initial configuration
{
    clear();                // remove all sites

    Site s;
    s.x = s.y = 0;
    addFront(s);

    for (int step = 1; step <= steps; step++) {

        int stp = 0, dir = EAST;    // initialize variables to construct coil

        switch (config) {
        case STAIR:             // add randomly East or North
            std_rand() < 0.5 ? ++s.x : ++s.y;
            break;
        case COIL:              // add in sequence E,N,W,W,S,S,E,E,E,N,N,N,...
            while (stp < step)
                stp += ++dir / 2;
                switch ((dir + 2) % DIRECTIONS) {
                case EAST  :  ++s.x;  break;
                case NORTH :  ++s.y;  break;
                case WEST  :  --s.x;  break;
                case SOUTH :  --s.y;  break;
                }
                break;
        case LINE:              // add East
        default:                // also the default
            ++s.x;
        }

        addFront(s);
    }
}

Site randomAllowed(         // return a random allowed site
    Site head,              // adjacent to this head site
    Site neck)              // excluding this neck site
{
    deque<Site> allowed;

    // find the 3 allowed directions and add site to allowed deque

    for (int direction = EAST; direction < DIRECTIONS; direction++) {

        Site s;
        s.x = head.x;
        s.y = head.y;

        switch (direction) {
        case EAST:
            if ( !(head.x == neck.x - 1 && head.y == neck.y) ) {
                ++s.x;
                allowed.push_back(s);
            }
            break;
        case NORTH:
            if ( !(head.x == neck.x && head.y == neck.y - 1) ) {
                ++s.y;
                allowed.push_back(s);
            }
            break;
        case WEST:
            if ( !(head.x == neck.x + 1 && head.y == neck.y) ) {
                --s.x;
                allowed.push_back(s);
            }
            break;
        case SOUTH:
            if ( !(head.x == neck.x && head.y == neck.y + 1) ) {
                --s.y;
                allowed.push_back(s);
            }
            break;
        }
    }

    // choose and return a random allowed site
    return allowed[int(2.999999 * std_rand())];
}

bool reptate() {            // attempt random move and return true if succeeded

    if (snake.size() < 2)               // cannot reptate
        return false;

    Site head, neck, sNext;

    if (std_rand() < 0.5) {    // choose front end of snake
        head = snake[0];
        neck = snake[1];
        sNext = randomAllowed(head, neck);
        if (occupied(sNext))
            return false;
        removeBack();
        addFront(sNext);
    } else {                            // choose back end of snake
        int n = snake.size();
        head = snake[n - 1];
        neck = snake[n - 2];
        sNext = randomAllowed(head, neck);
        if (occupied(sNext))
            return false;
        removeFront();
        addBack(sNext);
    }

    return true;
}

double rSquared() {         // end-to-end size squared
    if (snake.size() < 2)
        return 0.0;
    double dx = snake.front().x - snake.back().x;
    double dy = snake.front().y - snake.back().y;
    return dx * dx + dy * dy;
}

int main() {

    cout << " Reptation Method for Self-Avoiding Walks on a Square Lattice\n"
         << " ------------------------------------------------------------\n"
         << " Enter maximum number of steps in walk: ";
    int n_steps, n_walks;
    cin >> n_steps;
    cout << " Enter number of walks to generate: ";
    cin >> n_walks;
    cout << " Enter initial configuration 1 = stair, 2 = coil, 3 = line: ";
    int config;
    cin >> config;

    ofstream file("reptation.data");
    cout << " Steps   <r^2>       Std. Dev.   Success%  CPU secs" << endl;

    for (int steps = 1; steps <= n_steps; steps++) {

        double r2sum = 0;
        double r4sum = 0;
        int success = 0;
        clock_t startTime = clock();

        createSnake(steps, config);

        for (int i = 0; i < n_walks; i++) {
            if (reptate())
                ++success;
            double r2 = rSquared();
            r2sum += r2;
            r4sum += r2 * r2;
        }

        clock_t endTime = clock();

        double r2av = r2sum / n_walks;
        double stdDev = sqrt(r4sum / n_walks - r2av * r2av);
        double successPercent = success / double(n_walks) * 100.0;
        double cpu = (endTime - startTime) / double(CLOCKS_PER_SEC);

        cout << right << setw(4)  << steps           << "     "
             << left  << setw(10) << r2av            << "  "
             << left  << setw(10) << stdDev          << "  "
             << left  << setw(8)  << successPercent  << "  "
            << left  << setw(8)  << cpu             << '\n';
        file << steps << '\t' << r2av << '\t' << stdDev
             << '\t' << successPercent << '\t' << cpu << '\n';
    }

    file.close();
    cout << " Data in file reptation.data" << endl;
}
