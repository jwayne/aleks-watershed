aleks-watershed
========

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
bin/watershed [INPUT_FILE] [XSIZE] [YSIZE] [ZSIZE] [TYPE]
```
* INPUT_FILE: Affinity graph as a raw binary file of floats.  Example files can be downloaded at [http://microglia.mit.edu/data/].
* XSIZE
* YSIZE
* ZSIZE
* TYPE: Type of watershed algorithm to run.  All algorithms first run connected components with threshold T_h=.99, and then watershed to expand the seeds up to T_l=0.3.  Then, the types diverge:
  * 0: nothing else
  * 1: merge pairs of neighboring segments where some segment has size < watershed::limit_fn_avg(avg edge weight among edges straddling that pair of segments) [this does exactly what zi/watershed/just_watershed_avg.cpp does]
  * 2: merge pairs of neighboring segments where some segment has size < watershed::limit_fn_bup(max edge weight among edges straddling that pair of segments) [this does exactly what zi/watershed/just_watershed_bup.cpp does]
  * 3: merge pairs of neighboring segments where the max edge weight among edges straddling that pair of segments > 0.1 [this does exactly what zi/watershed/just_watershed.cpp does]
  * Note that the algorithm in Aleks' 2011 master's thesis (merge pairs of neighboring segments where some segment has size < T_s and the max straddling edge > T_e) is not available through this code.  However, types 1-2 are intended to be updates of that algorithm.  Furthermore, that algorithm was not available in this package before I modified it either.  That said, if desired, this should be easy to implement.

Once compiled and executed, the results can be viewed using my other package [https://github.com/jwayne/iw-seung].  Alternatively, the results can be viewed in MATLAB with the following commands:

```
f = fopen('./aleks-160-0.out', 'r');
v = fread(f, 'uint32');
v = reshape(v, [160 160 160]);
fclose(f);
flyOver(v);
```
