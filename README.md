# Getting Started

Compile via..?

Download the data files at [http://microglia.mit.edu/data/]

Execute via..?

Once compiled and executed, the results can be viewed in MATLAB with the
following command:

```
f = fopen('./vout.out', 'r');
v = fread(f, 'uint32');
v = reshape(v, [160 160 160]);
fclose(f);
flyOver(v);
```
