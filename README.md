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

### FEM Algorithm <a name="algo"></a>

The **FEM** is a numerical method to solve partial differential equations.

1. **Node class**:
   - You should enter the coordinates, boundary conditions and loading conditions.
   - Note that BCs is a 1*6 array containing True and False. True is for known and False for unknown DoFs.
   - Loads is also a 1*6 array, the first three of which are loads and the other are moments.
   - You should enter the id of that node which is an integer. (i.e 1 or 2 ...)
   - Finally, this class gives you the coordinates, Bcs, Loads and id of that specific node.
2. **Element class**:
   - First, you should enter the mechanical properties of the element.
   - You should also enter an array containing the id of both nodes in the element.
   - Using Element.el_info() method, you can get the stiffness matrix and the number of that element.
3. **Fem class$**:
   - Depending on how many elements you have, you should make an array for stiffness matrices of the elements, boundary conditions and loads of the nodes.
   - Finally, you need to creat an id array for the connectivity of elements.
   - The class gives the assembled stiffness matrix, BCs and Loading condition arrays of whole nodes.
4. **Repeat**:

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