# Cylinder Flow Simulation

## Pre-requisites

- `Eigen` Library：[libeigen / eigen · GitLab](https://gitlab.com/libeigen/eigen)。
- `SYCL` Library：、`icpx` Compiler：[Intel® oneAPI Base Toolkit: Essential oneAPI Tools & Libraries](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/base-toolkit.html)。

Need to generate the calculation grid (like `gen_gird.cpp`) to `X.txt` and `Y.txt`.

## Build & Run

### On Linux

```sh
icpx -fsycl src/main.cpp
./a.out
```

If you encounter the following problems:

```
error: Double type is not supported on this problem.
in kernel: 'xxxxxx'
error: backend compiler failed build.
```

Enter the following command on the command line before building the project:

```bash
export OverrideDefaultFP64Settings=1 
export IGC_EnableDPEmulation=1
```

### In Intel® DevCloud

```sh
chmod 755 q; chmod 755 run.sh;if [ -x "$(command -v qsub)" ]; then ./q run.sh; else ./run.sh; fi
```


## Params

```c++
const int DIV = 1000; // The number of divisions per second
const int SEC = 200; // Total simulation time (seconds)

const double Re = 200; // Reynolds number of the fluid
const double err = 0.005; // Minimum error value required for psi iterations
```

## Output

Output to the `./output/` folder in the same directory, including `PSI`, `U`, `V` three kinds of data, the number contained in the file name represents the number of seconds of simulation.

You can use `Matlab` or `Python` for data analysis, such as using `plot_flt.m` and `rainbow_figure.m` in this repository.

The following is a schematic diagram of the speed when `Re = 200, T = 100s` (the color depth represents the speed):

![schematic diagram of the speed](./img/spd_cl.jpg)