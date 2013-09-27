#ifndef CPT_GRAPHICS_HPP
#define CPT_GRAPHICS_HPP

#include <string>
#include <vector>

namespace cpt
{

    // Model-View-Controller classes

    class Model {
    public:
        Model() { }
    private:
    };

    class View {
    public:
        View() { }
    private:
    };

    class Controller {
    public:
        Controller() { }
    private:
    };

    // Gnuplot plots and movies

    class Gnuplot {
    public:
        Gnuplot() : frame_rate(24) { }
        void set_step_callback(void (*func)()) { step_callback = func; }
        void add_command(std::string cmd) { commands.push_back(cmd); }
        void clear_commands() { commands.clear(); }
        void animate(int steps);
        void set_frame_rate(double rate) { frame_rate = rate; }

    private:
        void (*step_callback)();
        std::vector<std::string> commands;
        double frame_rate;                      // frames per second

    };

} /* namespace cpt */

#endif /* CPT_GRAPHICS_HPP  */
