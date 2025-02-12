## Instructions for Building *HSCA*

```sh
sh build.sh
```

By executing this script file, users can build both *SamplingCA* and *HSCA*. Note that both *SamplingCA* and *HSCA* should be built on a 64-bit GNU/Linux operating system.

**Tip**: Both *SamplingCA* and *HSCA* depend on *MiniSAT*, but *MiniSAT* may not compile successfully when using specific versions of gcc. In this case, users may seek for solutions on the Internet or in the Github page of [MiniSAT](https://github.com/niklasso/minisat). Specifically, in our experiments, the version of `GCC` is 9.4.0, and the version of `Zlib` is 1.2.11.

We also provide a Docker image for running *HSCA*.

To build and run *HSCA* with Docker, users can use the following commands:

1. Load the Docker image from the tar file:
    ```
    docker load -i hsca-img.tar
    ```

2. Run the Docker container using the loaded image:
    ```
    docker run -d hsca:latest
    ```

## Quick Installation Testing

To check whether the compilation is successful or not, the user may run the following command in the root directory:

```
python3 run.py benchmarks/real-world/spins_4wise.model benchmarks/real-world/spins.constraints -seed 1
```

For Docker users:：

```
# Note: the benchmarks provided in this repo are already copied to the docker image,
# so you can directly use the following command to run *HSCA* with docker. (no more volume mapping is needed)
docker run --rm hsca:latest benchmarks/real-world/spins_4wise.model benchmarks/real-world/spins.constraints -seed 1
```

The command above calls *HSCA* to solve the instance `benchmarks/real-world/spins_4wise.model` with default hyper-parameter setting (L = 5,000). And here *HSCA* use the random seed of 1.

The console output is expected to be similar with the following text:
```
...
...  <-- possibly some huge amount of output by Coprocessor
...
1: 561
2: 1056
3: 1359
4: 1594
5: 1721
6: 1816
7: 1882
8: 1945
9: 1993
10: 2032
11: 2057
12: 2080
13: 2100
14: 2116
15: 2131
16: 2143
17: 2154
18: 2162
19: 2169
20: 2175
21: 2181
22: 2186
23: 2190
24: 2193
25: 2195
26: 2197
27: 2198

c Clear final: fix-up all remaining tuples ...
c All possible 2-tuple number: 2198
c 2-tuple coverage: 1
c Generate testcase set finished, containing 27 testcases!
c CPU time cost by generating testcase set: 0.024 seconds
c Testcase set saved in tmp/04d1e6fe89c7cf16978e33a918de9c50dcc7ce879c4238b7cf3a4b70c82fa2b9.out_boolean_CA.out
c 2-tuple number of generated testcase set: 2198
...
...  <-- possibly some huge amount of output by Coprocessor
...
Expandor init success
opt_method = 1
------------------------ strengh = 3 ------------------------
num_combination_all_possible_ = 47872
covered valid 3-wise tuple nums: 11932
uncovered valid 3-wise tuple nums: 2861
all valid 3-wise tuple nums: 14793
invalid 3-wise tuple nums: 519
all tuples: 15312
c current test suite size: 28, current uncovered valid 3-wise tuple nums: 2694
c current test suite size: 29, current uncovered valid 3-wise tuple nums: 2534
c current test suite size: 30, current uncovered valid 3-wise tuple nums: 2397
c current test suite size: 31, current uncovered valid 3-wise tuple nums: 2263
c current test suite size: 32, current uncovered valid 3-wise tuple nums: 2139
...
... <-- Process of Generating 3-wise CA
...
c current test suite size: 122, current uncovered valid 3-wise tuple nums: 2
c current test suite size: 123, current uncovered valid 3-wise tuple nums: 1
c current test suite size: 124, current uncovered valid 3-wise tuple nums: 0
uncovered_nums: 0
add new testcase: 97
uncovered_nums: 0
------------------------ strengh = 4 ------------------------
num_combination_all_possible_ = 742016
covered valid 4-wise tuple nums: 128884
uncovered valid 4-wise tuple nums: 13118
all valid 4-wise tuple nums: 142002
invalid 4-wise tuple nums: 9694
all tuples: 151696
c current test suite size: 125, current uncovered valid 4-wise tuple nums: 12820
c current test suite size: 126, current uncovered valid 4-wise tuple nums: 12546
c current test suite size: 127, current uncovered valid 4-wise tuple nums: 12278
c current test suite size: 128, current uncovered valid 4-wise tuple nums: 12032
...
... <-- Process of Generating 4-wise CA
...
c current test suite size: 474, current uncovered valid 4-wise tuple nums: 1
c current test suite size: 475, current uncovered valid 4-wise tuple nums: 0
uncovered_nums: 0
add new testcase: 351
uncovered_nums: 0
the first pass of HSCA's optimization approach
Optimizer init success
c current remove 88
c current 4-wise CA size: 474, step #0
c current remove 138
c current 4-wise CA size: 472, step #0
c current remove 139
c current 4-wise CA size: 470, step #0
c current remove 198
c current 4-wise CA size: 468, step #0
c current 4-wise CA size: 467, step #0
c current 4-wise CA size: 466, step #1
c current 4-wise CA size: 465, step #2
c current 4-wise CA size: 464, step #6
c current 4-wise CA size: 463, step #9
...
... <-- The process of HSCA's first optimization pass
...
c current 4-wise CA size: 369, step #5314
c current 4-wise CA size: 368, step #6138
c current 4-wise CA size: 367, step #6532
final 4-wise CA size is: 367
End
the second pass of HSCA's optimization approach
c current 4-wise CA size: 367, step #0, time #0.051086
c current 4-wise CA size: 366, step #18, time #0.052989
c current 4-wise CA size: 365, step #35, time #0.055852
c current 4-wise CA size: 364, step #58, time #0.057852
...
... <-- The process of HSCA's first optimization pass
...
c current 4-wise CA size: 309, step #117017, time #21.1948
c current 4-wise CA size: 308, step #117169, time #21.2256
```

We finally note that due to potential differences in the random number generation mechanism over different versions of g++, the console output and the size of the generated 4-wise CA on the user's machine may be slightly different from the one presented here.
