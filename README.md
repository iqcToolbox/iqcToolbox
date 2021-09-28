![GitHub](https://img.shields.io/github/license/iqcToolbox/iqcToolbox?color=green)
[![Tests](https://github.com/iqcToolbox/iqcToolbox/actions/workflows/ci.yml/badge.svg)](https://github.com/iqcToolbox/iqcToolbox/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/iqcToolbox/iqcToolbox/branch/develop/graph/badge.svg?token=C3E2E3V80K)](https://codecov.io/gh/iqcToolbox/iqcToolbox)
![PRs Welcome](https://img.shields.io/badge/PRs-welcome-green)
![GitHub commits since latest release (by date) for a branch](https://img.shields.io/github/commits-since/iqcToolbox/iqcToolbox/latest/develop) (on `develop` branch)

# iqcToolbox

MATLAB toolbox for modeling uncertain/nonlinear systems and using integral quadratic constraint theory for conducting robustness analysis (stability, performance, and reachability).

## Installation

Either clone this repository (`git clone https://github.com/iqcToolbox/iqcToolbox.git`)
or download and unzip the [zipball](https://github.com/iqcToolbox/iqcToolbox/zipball/master) into a directory.

In MATLAB simply run the function `installIqcToolbox` located in `name_of_toolbox_directory/scripts/`.

This function will request permission to download and install the required dependencies [yalmip](https://yalmip.github.io/), [SDPT3](https://github.com/Kim-ChuanToh/SDPT3/), and [LPSOLVE](http://lpsolve.sourceforge.net/5.5/).  You may benefit from having other SDP solvers besides `sdpt3`.  This might include [SeDuMi](https://github.com/sqlp/sedumi), [mosek](https://www.mosek.com/), or [csdp](https://github.com/coin-or/Csdp).  The latter two solvers are capable of leveraging parallel computing resources; MOSEK is free for academic use, CSDP is free for general use (EPL 2.0 license).

## Tests

After installation, all the tests may be run in MATLAB with `tests.run_tests_script`.

A specific test in the `+test` directory may be run by executing the first section in `run_tests_script.m`, and then calling
```
result = runner.run(matlab.unittest.TestSuite.fromClass(?tests.SUBPACKAGE.TESTCLASS));
```
For example, to run the tests defined by the `testUlft` class:
```
result = runner.run(matlab.unittest.TestSuite.fromClass(?tests.unit.testUlft));
```

## Purpose and capabilities of toolbox

This toolbox is capable of both modeling and analyzing uncertain/nonlinear systems which can be modeled as an interconnection of a linear time-varying system `M` (describing the nominal system) and a structured set of bounded operators `Delta` (describing the system's uncertainties and nonlinearities). For example, to determine if the system described by the dynamical equations
```
x(k+1) = del * x(k) + u(k)
y(k)   = x(k)
```
is stable for all `del \in [-0.9, 0.9]`, we can build the uncertain system (linear fractional transformation, LFT) with the commands
```
del_dimension = 1;
del = toLft(DeltaSlti('del', del_dimension, -0.9, 0.9)); % Creating the uncertainty "del"
a = del; 
b = 1; 
c = 1; 
d = 0;
timestep = -1;
g = toLft(a, b, c, d, timestep); % Creating the uncertain system 
```
and determine if the uncertain system is robustly stable with the command
```
result = iqcAnalysis(g)
```
Note, if `del \in [-1, 1]`, the uncertain system would not be robustly stable, and IQC analysis would return an "infeasible problem" result.

A strength in IQC analysis is the flexibility afforded in analyzing many different types of uncertainties, with an extensive library of so-called "IQC multipliers" which are used to characterize different types of uncertainties.  For example, there is oftentimes uncertainty in the body-relative location of a UAV's center-of-gravity, which could potential destabilize the aircraft.  A far more challenging problem is if the UAV's center-of-gravity is not just uncertain, but also changing in flight (when, for example, fuel sloshes or is consumed, or when payloads are repositioned or ejected). IQC analysis can leverage properties of uncertainties (such as time-invariance, or frequency-domain characteristics) to reduce conservatism when conducting worst-case analysis of a system.

For example, the system
```
  ┌─     ─┐   ┌─                 ─┐ ┌─   ─┐   ┌─ ─┐
  │x1(k+1)│   │    0     3/4 + del│ │x1(k)│   │ 1 │
  │       │ = │                   │ │     │ + │   │ u(k)
  │x2(k+1)│   │3/4 - del     0    │ │x2(k)│   │ 0 │
  └─     ─┘   └─                 ─┘ └─   ─┘   └─ ─┘

                                    ┌─   ─┐
                          ┌─     ─┐ │x1(k)│
     y(k)   =             │ 1   0 │ │     │
                          └─     ─┘ │x2(k)│
                                    └─   ─┘
```
where `del \in [-1/2, 1/2]`, can be proven to be stable for any admissible instance of `del`, as long as it is time-invariant.  This can be checked with the following commands:
```
del = toLft(DeltaSlti('del', 1, -1/2, 1/2)); % Create an LFT object modeling the uncertainty del
a = [0, 3/4 + del; 3/4 - del, 0];
b = [1; 0];
c = [1, 0];
d = 0;
timestep = -1; % Discrete-time system with unspecifed timestep
g = toLft(a, b, c, d, -1);
result = iqcAnalysis(g);
result.valid
```
Indeed, the eigenvalues of the a-matrix are `sqrt(9/16 - del^2)`, which are strictly within the unit circle (meaning that the a-matrix is "Schur"), implying stability of the system for any `del \in [-1/2, 1/2]`.

However, if `del` could possibly be time-varying, robust stability cannot be ensured.
```
a.delta.deltas{1} = DeltaSltv('del', 2, -1/2, 1/2);
g_tv = toLft(a, b, c, d, timestep);
result_tv = iqcAnalysis(g_tv);
result_tv.valid
```
This is expected. If, for example, `del(k) = 1/2 for odd k` and `del(k) = -1/2 for even k`, then the aforementioned system is unstable (example inspired by [1]). IQC analysis is capable of leveraging the time-invariance property of `del`to prove robust stability, and rightly fails to prove robust stability if `del` were a member of the more general set of time-varying parameters. This illustrates the benefit of using IQC multipliers to characterize as tightly as possible the sets of uncertainties and disturbances that affect a system.

Much more complicated systems may be analyzed, with time-varying nominal systems, and a mixture of numerous uncertainties/nonlinearities (rate-bounded time-varying parameters, sector-bounded nonlinearities, time delays, and more) [2].  IQC analysis can also consider different performance measures, limit the admissible disturbances to specific sets (frequency-banded signals, white noise, periodically increasing/decreasing signals, etc.), incorporate uncertain initial conditions, and conduct reachability analysis [3, 4].  

As the theoretical body regarding IQC analysis is large, we welcome contributions that enable the incorporation of uncertainties, disturbances, and analytical capabilities that are currently absent from this toolbox.


## Contributing

Please refer to the [contributing guidelines](https://github.com/iqcToolbox/iqcToolbox/blob/master/CONTRIBUTING.md).

## References

[1] Elaydi, Saber N. "An Introduction to Difference Equations,(1996)." Springer-Verlag, New York.(10) M. Hirsch, J. Palis, G. Pugh and M. Shub, Neighborhoods of hyperbolic sets, Invent. Math 9 (1970): 121-13

[2] Megretski, Alexandre, and Anders Rantzer. "System analysis via integral quadratic constraints." IEEE Transactions on Automatic Control 42, no. 6 (1997): 819-830.

[3] Veenman, Joost, Carsten W. Scherer, and Hakan Köroğlu. "Robust stability and performance analysis based on integral quadratic constraints." European Journal of Control 31 (2016): 1-32.

[4] Fry, J. Micah, Dany Abou Jaoude, and Mazen Farhood. "Robustness analysis of uncertain time‐varying systems using integral quadratic constraints with time‐varying multipliers." International Journal of Robust and Nonlinear Control 31, no. 3 (2021): 733-758.

## License Header

Copyright (C) 2021 Massachusetts Institute of Technology

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

A full license statement is provided [here](https://github.com/iqcToolbox/iqcToolbox/blob/release/LICENSE.md)

## Distribution Statement

DISTRIBUTION STATEMENT A. Approved for public release: distribution unlimited.

© 2021 MASSACHUSETTS INSTITUTE OF TECHNOLOGY

Subject to FAR 52.227-11 – Patent Rights – Ownership by the Contractor (May 2014)
SPDX-License-Identifier: GPL-2.0

The software/firmware is provided to you on an As-Is basis
