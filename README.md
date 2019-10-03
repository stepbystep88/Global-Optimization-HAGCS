# Global-Optimization-HAGCS

## Introduction
A new global optimization algorithm developed by using gradient infromation, multiple search strategies, and parameters self-adaptation. 

**HAGCS** is desgined for the objective functions whose gradient is readily available. Essentially, the **HAGCS** algorithm is based on the Cuckoo search [CS, 1] algorithm.

## Installation
addpath(genpath("./WZDOM"));

## Test

1. View the animation of optimization process (run test_global_optimization/testShowAnimation.m )

> Compare `CS` with `HAGCS` on the Drop Wave test function:
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Drop%20Wave_25.gif)

> Compare `CS` with `HAGCS` on the Drop Wave test function (bigger population size):
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Drop%20Wave_100.gif)

> Compare `CS` with `HAGCS` on the Rastrigin test function:
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Rastrigin_25.gif)

> Compare `CS` with `HAGCS` on the Rastrigin test function (smaller population size):
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Rastrigin_5.gif)

> Compare `CS` with `HAGCS` on the Schwefel test function:
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Schwefel.gif)

> Compare `CS` with `HAGCS` on the Schwefel's P2.22 test function:
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Schwefel's%20P2.22.gif)

> Compare `CS` with `HAGCS` on the Ackley test function:
> ![alt text](https://github.com/stepbystep88/Global-Optimization-HAGCS/blob/master/test_global_optimization/frames/Cmp_CS_HAGCS_Ackley.gif)

2. To compare the performance of finding global minimum, convergence speed of different algorithms, please run ./test_global_optimization/testCompareMethodsFinal.m

## References
[1] Xin-She Yang and Suash Deb. Cuckoo search via levy flights. In 2009 World Congress on Nature Biologically Inspired Computing (NaBIC), pages 210-214, Dec 2009.
