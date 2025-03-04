# Finite Element Method - Part1

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/mrahimi74/FEM/graph/badge.svg?token=RqAIkRNOJH)](https://codecov.io/gh/mrahimi74/FEM)
[![tests](https://github.com/Lejeune-Lab-Graduate-Course-Materials/bisection-method/actions/workflows/tests.yml/badge.svg)](https://github.com/Lejeune-Lab-Graduate-Course-Materials/bisection-method/actions)
---

### Table of Contents
* [Getting Started](#gs)
* [FEM algorithm](#algo)
* [Conda environemnt, installation, and testing](#install)
* [Tutorial](#tutorial)

---

### Getting Started

To be written

---

### Bisection Method Algorithm <a name="algo"></a>

The **Bisection Method** is a numerical technique to find roots of a continuous function f(x). The method works by repeatedly dividing an interval [a, b] in half and selecting the subinterval in which the root lies. The algorithm for the Bisection Method is as follows:

1. **Choose an interval $[a, b]$**:
   - choose $[a, b]$ such that $f(a) \cdot f(b) < 0$, which means the function has a root in the interval.
2. **Compute the midpoint**:
   - $c = \frac{a + b}{2}$
3. **Check the sign of $f(c)$**:
   - If $f(c) = 0$, then $c$ is the root.
   - If $f(a) \cdot f(c) < 0$, the root lies in $[a, c]$. Update $b = c$.
   - Otherwise, the root lies in $[c, b]$. Update $a = c$.
4. **Repeat** until the interval size $|b - a|$ is smaller than the desired tolerance or $f(c)$ is smaller than the desired tolerance.

**Advantages of the Bisection Method**:
1. **Simplicity**: The method is easy to understand and implement.
2. **Guaranteed Convergence**: If the function $f(x)$ is continuous and $f(a) \cdot f(b) < 0$, the method is guaranteed to converge to a root.
3. **Robustness**: It works well for a wide range of functions without requiring derivatives or complex calculations.
4. **Predictable Behavior**: The error decreases by approximately half in each iteration, providing predictable convergence.

**Limitations of the Bisection Method**:
1. **Slow Convergence**: Compared to other methods like Newton-Raphson or secant, the bisection method can converge more slowly.
2. **Requires Bracketing**: The method requires an interval $[a, b]$ where $f(a) \cdot f(b) < 0$, which may not always be easy to identify.
3. **Limited Precision**: The method converges linearly, which may not be efficient for high-precision requirements.
4. **Not Suitable for Multiple Roots**: The method may fail or behave inconsistently if the function has multiple roots within the interval.
5. **Function Continuity Required**: It assumes $f(x)$ is continuous in $[a, b]$, and any discontinuities can cause issues.

---

### Conda environment, install, and testing <a name="install"></a>

To install this package, please begin by setting up a conda environment (mamba also works):
```bash
conda create --name bisection-method-env python=3.12
```
Once the environment has been created, activate it:

```bash
conda activate bisection-method-env
```
Double check that python is version 3.12 in the environment:
```bash
python --version
```
Ensure that pip is using the most up to date version of setuptools:
```bash
pip install --upgrade pip setuptools wheel
```
Create an editable install of the bisection method code (note: you must be in the correct directory):
```bash
pip install -e .
```
Test that the code is working with pytest:
```bash
pytest -v --cov=bisectionmethod --cov-report term-missing
```
Code coverage should be 100%. Now you are prepared to write your own code based on this method and/or run the tutorial. 

If you would like, you can also open python and check to make sure that the import works properly:
```bash
(bisection-method-env) $ python
Python 3.12.8 | packaged by Anaconda, Inc. | (main, Dec 11 2024, 11:37:13) [Clang 14.0.6 ] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from bisectionmethod import bisection_method as bim
>>> bim.hello_world()
'hello world'
```
If you are using VSCode to run this code, don't forget to set VSCode virtual environment to bisection-method-env.

If you would like the open `tutorial.ipynb` located in the `tutorials` folder as a Jupyter notebook in the browser, you might need to install Jupyter notebook in your conda environment as well:
```bash
pip install jupyter
```
```bash
cd tutorials/
```
```bash
jupyter notebook tutorial.ipynb
```
---

### Tutorial <a name="tutorial"></a>

#### **What Does the Function Do?**

The `run_bisection_method` function takes:
- A continuous function $f(x)$,
- Two bounds $a$ and $b$ where $f(a)$ and $f(b)$ have opposite signs, and
- A tolerance and a maximum number of iterations.

It iteratively computes the root of the function using the bisection method and returns the root along with detailed iteration data.

---

#### **Inputs and Outputs**

#### **Inputs**
1. **`fcn`**: A Python callable (function) that represents $f(x)$.
2. **`a`**: The lower bound of the interval.
3. **`b`**: The upper bound of the interval. Must satisfy $f(a) \times f(b) < 0$ (i.e., the root lies between $a$ and $b$.
4. **`tol_input`**: (Optional) The tolerance for the interval size. Default is $10^{-9}$.
5. **`tol_output`**: (Optional) The tolerance for the function output. Default is $10^{-30}$.
6. **`max_num_iter`**: (Optional) The maximum number of iterations to perform. Default is $1000$.

#### **Outputs**
The function returns a dictionary with the following keys:
- **`solution`**: The computed root.
- **`num_iter`**: The number of iterations performed.
- **`all_a`**: A list of all the intermediate $a$ values.
- **`all_fcn_a`**: A list of function values corresponding to $a$.
- **`all_b`**: A list of all the intermediate $b$ values.
- **`all_fcn_b`**: A list of function values corresponding to $b$.

---

### **Summary of Errors and Their Causes**

The `run_bisection_method` function includes several checks to ensure valid input and proper conditions for the bisection method. 

| **Error Message**                                                                 | **Cause**                                                                                  | **Solution**                                                                                  |
|-----------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| `Invalid input: {a} is greater than {b}.`                                         | $a \geq b$, invalid interval.                                                         | Ensure $a < b$.                                                                           |
| `a and b are not guaranteed to contain a root of the continuous function provided` | $f(a)$ and $f(b)$ have the same sign, root not guaranteed in the interval.        | Choose bounds where $f(a) \times f(b) < 0$.                                               |
| `Maximum number of iterations ({max_iter}) reached without convergence`           | Root not found within the specified maximum number of iterations.                         | Increase the `max_num_iter` parameter or check the function for potential issues.             |
| `The function evaluations must have one positive and one negative value.`         | The function is discontinuous, or bounds do not guarantee a root.                        | Ensure the function is continuous and the root lies between $a$ and $b$.              |

By handling these errors carefully, you can debug and ensure proper usage of the bisection method for finding roots in a wide variety of scenarios.

---

#### **Examples**

After following the installation instructions above, it will be possible to run the tutorial examples contained in the `tutorials` folder.

---

### More information <a name="more"></a>
More information can be found here:
* https://en.wikipedia.org/wiki/Bisection_method
* https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html
* https://en.wikipedia.org/wiki/Root-finding_algorithm