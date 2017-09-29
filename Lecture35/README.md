# Using N-body starter code by Hut and Makino for N-body simulation

Full details are [here](https://www.ids.ias.edu/~piet/act/comp/algorithms/starter).

## Download, compile, execute figure 8

### Get the code:

- Linux:
```
wget https://www.ids.ias.edu/sites/ids.ias.edu/files/imported/act/comp/algorithms/starter/nbody_sh1.tar.gz
```
- Mac:
```
curl -O https://www.ids.ias.edu/sites/ids.ias.edu/files/imported/act/comp/algorithms/starter/nbody_sh1.tar.gz
```

### Unpack the code:
```
tar -zxvf nbody_sh1.tar.gz
```

### Compile:
```
g++ nbody_sh1.C -o nbody_sh1
```

### Execute:
```
./nbody_sh1 -o 0.01 < figure8.in > figure8.out
```

## Generate random:

### 

## Animate:
```
python nbody_animate.py
```
