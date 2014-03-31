Segment a 3D affinity graph using various modified watershed algorithms.  Code thanks to Aleks Zletaski of MIT, my work has been to prepare watershed for simpler compilation/execution.

# Getting Started

Get the dependencies:
* g++4.7 (compiler for C++11)
* boost (C++ library)

Compile via:
```
g++-4.7 -Wall -std=c++0x zi/watershed/watershed.cpp -I/usr/local/lib/boost -I. -o bin/watershed
```

Execute via:
```
bin/watershed [INPUT_FILE] [INPUT_FILE_DIMENSIONS] [TYPE]
```
* INPUT_FILE: Affinity graph as a raw binary file of floats.  Example files can be downloaded at [http://microglia.mit.edu/data/].
* INPUT_FILE_DIMENSIONS: Dimensions of the affinity graph, assuming a cube.
* TYPE: Type of watershed algorithm to run.  See comments in zi/watershed/watershed.cpp for details.

Once compiled and executed, the results can be viewed in MATLAB with the
following commands:

```
f = fopen('./aleks-160-0.out', 'r');
v = fread(f, 'uint32');
v = reshape(v, [160 160 160]);
fclose(f);
flyOver(v);
```
