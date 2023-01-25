# Numerical-Methods-For-Solving-A-System-Of-Linear-Equations
Numerical methods for solving a system of linear equations

## User Manual
Run GUI.py
choose the method for solving
enter each equation seperately then press ok
enter the required parameters then solve
![img](images/gui1.jpg)

or you can put the data in a file
Reading from files will follow the following template:
1st line: number of equations
2nd line: Method Name (e.g. ‘Gaussian-elimination)
3rd line --> nth line: equations
last line: Space separated initial points (e.g. 1.1 2)

A sample file in this case will be: 
3
Gaussian-elimination 
3*a + 2*b + c - 6 
2*a + 3*b - 7 
2*c - 4 
![img](images/gui2.jpg)

### Analysis for the behavior of different examples:
❖ Gaussian Elimination
Example 1) 4*a-b+c-4
a+6*b+2*c-9
-1*a-2*b+5*c-2
![img](images/gauss1.jpg)

❖ LU decomposition
Example) 1*a+1*b+1*c-1
3*a+1*b-3*c-5
1*a-2*b-5*c-1
![img](images/lu1.jpg)

❖ Gaussian Jordan
Example) 1x+3y-7
3x+4z-11
![img](images/gj1.jpg)

❖ Gauss Seidel
Example) 4*a-b+c-4
a+6*b+2*c-9
-1*a-2*b+5*c-2
![img](images/gs1.jpg)

❖ Jacobi
Example) 25a+5b-1c-106.8
-4.8b+2c+96.21
0.7c-0.735
![img](images/jacobi.jpg)

❖ All methods at once
Example) 12x+3y−5z-1
x+5y+3z-28
3x+7y+13z-76
Initial val = [1 0 1]
![img](images/all1.jpg)
![img](images/all2.jpg)
![img](images/all_roots.jpg)

