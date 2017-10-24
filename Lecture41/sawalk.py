import math
import random

print(" Self-Avoiding Walks on a Square Lattice")
print(" ---------------------------------------")
n_steps = int(input(" Enter number of steps in walk: "))
n_walks = int(input(" Enter number of walks to generate: "))

walks = 0
failed_walks = 0
r2_av = 0
r4_av = 0

# generate walks
while walks < n_walks:

    sites = []          # list of occupied lattice sites
    x, y = [0, 0]       # starting site
    sites.append([x, y])
    walk_failed = False

    # loop over desired number of steps
    for step in range(n_steps):

        # take a random step
        d = random.random()
        if d < 0.25:
            x += 1      # step East
        elif d < 0.50:
            y += 1      # step North
        elif d < 0.75:
            x -= 1      # step West
        else:
            y -= 1      # step South

        # check whether the site is occupied
        if [x, y] in sites:
            walk_failed = True
            failed_walks += 1
            break
        else:
            sites.append([x, y])

    if walk_failed:
        continue        # try again

    r2 = x**2 + y**2
    r2_av += r2
    r4_av += r2**2
    walks += 1

r2_av /= n_walks
r4_av /= n_walks
std_dev = math.sqrt(r4_av - r2_av**2)
total_walks = n_walks + failed_walks
failed_percent = failed_walks / total_walks * 100.0
print(" Mean square distance <r^2> =", r2_av)
print(" Standard deviation         =", std_dev)
print(" Number of walks attempted  =", total_walks)
print(" Percentage failed walks    =", failed_percent)
