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

## Technical Details

The specific optimization is reflected in the three functions of the vorticity FTCS iteration `push()`、the flow function iteration `psi_iteration()`、and the velocity solution `velocity()`, because these three parts perform a large number of matrix operations, and the amortized complexity of each iteration is $\mathcal O(n^2)$. Therefore, consider using [SYCL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/data-parallel-c-plus-plus.html) for heterogeneous computing. The other function `boundary()` in the same iteration has a complexity of only $\mathcal O(n)$, no need for special optimization.

Here is the iterated discrete form of the stream function $\psi$，Since the values of $\alpha, \beta, \gamma, \delta, \varepsilon$ are only determined by the grid information，Therefore, preprocessing can be performed at the beginning of the program. Although this approach increases the addressing time, it reduces the number of floating-point calculations.

$\begin{matrix} \psi_{i, j}^{k+1}&=\frac{1}{2(\alpha+\gamma)}\left[\left(\alpha-\frac{\delta}{2}\right) \psi_{i-1, j}^{k}+\left(\alpha+\frac{\delta}{2}\right) \psi_{i+1, j}^{k}+\left(\gamma-\frac{\varepsilon}{2}\right) \psi_{i, j-1}^{k}+\left(\gamma+\frac{\varepsilon}{2}\right) \psi_{i, j+1}^{k}\right] \\  &+\frac{1}{2(\alpha+\gamma)}\left[2 \beta\left(\psi^{k}{ }_{i+1, j+1}-\psi^{k}{ }_{i+1, j-1}-\psi^{k}{ }_{i-1, j+1}+\psi^{k}{ }_{i-1, j-1}\right)-\omega_{i, j}\right]\end{matrix}$

`MatrixXd` defaults to column priority, which leads to discontinuous memory in each row when converted to an array by `MatrixXd::data()`. I exchange the summation symbol and transform the formula to traverse by column to increase addressing efficiency.