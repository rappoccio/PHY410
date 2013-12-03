import math
import random
from Tkinter import *
import matplotlib.pyplot as plt

class Disk:

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
        self.apply_pbcs()

    def apply_pbcs(self):
        while self.x >= 1: self.x -= 1
        while self.x <  0: self.x += 1
        while self.y >= 1: self.y -= 1
        while self.y <  0: self.y += 1

class System:

    def __init__(self, nu=0, N=224, rows=16, cols=14):
        self.disks = [ Disk() for i in range(N) ]
        d = 1.0 / float(cols)
        self.d_0 = d * (1 - 2**(nu-8.0))
        self.alpha = d - self.d_0
        dx = d
        dy = d * math.sqrt(3.0) / 2.0
        for j in range(rows):
            for i in range(cols):
                n = i + j * cols
                self.disks[n].x = i * dx + j % 2 * dx / 2.0
                self.disks[n].y = j * dy

    def overlap(self, disk1, disk2):
        dx = disk1.x - disk2.x
        dy = disk1.y - disk2.y
        if abs(dx) > 0.5:
            dx *= 1 - 1 / abs(dx)
        if abs(dy) > 0.5:
            dy *= 1 - 1 / abs(dy)
        return dx**2 + dy**2 < self.d_0**2

    def move_allowed(self, disk):
        x_try = disk.x + self.alpha * random.uniform(-1, 1)
        y_try = disk.y + self.alpha * random.uniform(-1, 1)
        moved_disk = Disk(x_try, y_try)
        for disk2 in self.disks:
            if disk2 is not disk and self.overlap(disk2, moved_disk):
                return False
        disk.x = moved_disk.x
        disk.y = moved_disk.y
        return True

    def cycle(self):
        accepts = 0
        for disk in self.disks:
            if self.move_allowed(disk):
                accepts += 1

    def pair_histogram(self, K=1.5, N_bins=64):
        histogram = [ 0.0 ] * N_bins
        DeltaA2 = (K**2 - 1) * math.pi * self.d_0**2 / float(N_bins)
        n_disks = len(self.disks)
        for i in range(n_disks - 1):
            for j in range(i + 1, n_disks):
                dx = self.disks[i].x - self.disks[j].x
                if abs(dx) > 0.5:
                    dx *= 1 - 1 / abs(dx)
                dy = self.disks[i].y - self.disks[j].y
                if abs(dy) > 0.5:
                    dy *= 1 - 1 / abs(dy)
                bin = int(math.pi *
                          (dx**2 + dy**2 - self.d_0**2) / DeltaA2)
                if bin >= 0 and bin < N_bins:
                    histogram[bin] += 1
        return histogram

class Animation(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        Pack.config(self)
        self.create_widgets()
        self.after(10, self.animation_step)

    def create_widgets(self):
        self.QUIT = Button(self, text='QUIT', command=self.quit)
        self.QUIT.pack(side=BOTTOM)
        self.draw = Canvas(self, width='400', height='400')
        self.draw.pack(side=TOP)

    def draw_disks(self):
        self.draw.delete("all")
        self.draw.create_rectangle(50, 50, 350, 350)
        for disk in system.disks:
            x = 50 + disk.x * 300
            y = 50 + disk.y * 300
            r = system.d_0 / 2.0 * 300
            self.draw.create_oval(x - r, y - r, x + r , y + r, fill='red')

    def animation_step(self):
        self.draw_disks()
        system.cycle()
        self.after(10, self.animation_step)

print(" Monte Carlo Simulation of Hard Disk Gas")
nu = float(input(" Enter value of nu [0..7]: "))
system = System(nu)

animate = False
try:
    response = input(' Run animation? 1 (yes) or 0 (no): ')
    if response :
       animate = True
except:
    pass
if animate:
    print 'Animating...'
    root = Tk()
    anim = Animation(master=root)
    anim.mainloop()
    root.destroy()
else:
    K = 1.5
    N_bins = 64
    histogram = [ 0.0 ] * N_bins
    xvals = [ i for i in xrange(N_bins) ]
    mc_cycles = 100
    print(" Running " + str(mc_cycles) + " Monte Carlo cycles ...")
    for step in range(mc_cycles):
        ph = system.pair_histogram(K, N_bins)
        for bin in range(N_bins):
            histogram[bin] += ph[bin]
        system.cycle()

    plt.scatter(xvals, histogram)

    plt.ylim(0, 4000)
    plt.xlim(0, N_bins)
    
    plt.show()

