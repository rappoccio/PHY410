import math
import random
import sys
import timeit
from collections import deque

def create_snake(steps, config='coil'):
    global snake
    snake = deque()         # double-ended priority queue

    # initialize with head and tail at x=0, y=0
    x, y = [0, 0]
    snake.append([x, y])

    # choice of 'coil', 'stair', or straight line
    if config == 'coil':
        # tight coil East, North, West^2, South^2, E^3, N^3, ...
        directions = [ ['East', 'North'], ['West', 'South'] ]
        ds = directions[0]  # step E then N
        n_ds = 1            # number of directed steps in each direction
        step = 0            # number of step done
        while step < steps:
            for i in range(2):
                for n in range(n_ds):
                    if ds[i] == 'East':
                        x += 1
                    elif ds[i] == 'North':
                        y += 1
                    elif ds[i] == 'West':
                        x -= 1
                    else:
                        y -= 1
                    if step < steps:
                        snake.append([x, y])
                    step += 1
            n_ds += 1       # increase number of directed steps
    elif config == 'stair':
        # step randomly East or North
        for step in range(steps):
            if random.randrange(2) == 0:
                x += 1
            else:
                y += 1
            snake.append([x, y])
    else:
        # straight line East
        for step in range(steps):
            x += 1
            snake.append([x,y])

def random_allowed(head, neck):
    # given x,y coordinates of head and neck
    # return a random site for the head to move to

    # construct a list of 4 nearest neighbors of the head and remove the neck
    x, y = head
    neighbors = [ [x + 1, y], [x - 1, y], [x, y + 1], [x, y -1] ]
    neighbors.remove(neck)

    # return a random allowed site
    return neighbors[random.randrange(3)]

def reptate_succeeded():    # attempt random move and return true if succeeded
    global snake

    # choose one end of the chain at random to be the head by tossing a coin
    if random.randrange(2) == 0:
        snake.reverse()     # switch head and tail

    # choose a random head move
    next = random_allowed(snake[0], snake[1])

    # check whether it is an interior site other than the neck
    tail = snake[-1]
    if next in snake and next != tail:
        return False        # move failed, retain current configuration

    # repate by removing the tail and advancing the head
    snake.remove(tail)
    snake.appendleft(next)
    return True

print(" Reptation Method for Self-Avoiding Walks on a Square Lattice")
print(" ------------------------------------------------------------")
n_steps = int(input(" Enter maximum number of steps in walk: "))
n_walks = int(input(" Enter maximum number walks to generate: "))
config = input(" Enter initial configuration (coil, stair, line): ")

file = open("reptation.data", "w")
print(" Steps   <r^2>       Std. Dev.   Success%  CPU secs")

for steps in range(1, n_steps + 1):
    r2_sum = 0.0
    r4_sum = 0.0
    success = 0
    start_time = timeit.default_timer()

    create_snake(steps, config)
    for i in range(n_walks):
        if reptate_succeeded():
            success += 1
            head = snake[0]
            tail = snake[-1]
            r2 = (head[0] - tail[0])**2 + (head[1] - tail[1])**2
            r2_sum += r2
            r4_sum += r2**2

    end_time = timeit.default_timer()
    r2_av = r2_sum / n_walks
    std_dev = math.sqrt(r4_sum / n_walks - r2_av**2)
    success_percent = success / float(n_walks) * 100.0
    cpu = (end_time - start_time)

    result = ('  '
            + str(steps).ljust(4) + '   '
            + str(r2_av)[:10].ljust(10) + '  '
            + str(std_dev)[:10].ljust(10) + '  '
            + str(success_percent)[:8].ljust(8) + '  '
            + str(cpu))
    print(result)
    file.write(result + '\n')

file.close()
print(" Data in file reptation.data")
