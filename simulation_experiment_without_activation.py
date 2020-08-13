import numpy as np;
import csv;
import scipy;
from scipy import stats;
import matplotlib.pyplot as plt;
from scipy.optimize import fsolve, root, broyden1;
from sympy import Symbol,nsolve;

import reactorbioquim;

K=10;
k=1;
kd=0.022;
e0=18;
Si=200;
Xopt=0.9;
nlotes=6;
kr=0.5;
V=1;
Mbiocat=1;
elimit=0.25

efinish, t_reaccion, product, productivity_specif, Productividad_especifica_global = reactorbioquim.inactivity(K, k, kd, e0, Si, Xopt, nlotes, V, Mbiocat, elimit);

print('Producto',product);
print('Productividad reacción',productivity_specif);
print('Productividad_especifica_global',Productividad_especifica_global);

print(efinish);
print(t_reaccion);
print(product);

x = np.arange(0,nlotes+1,1)
plt.plot(x,efinish);
plt.xlabel('lotes')
plt.title('Actividad Específica')
plt.show();

plt.plot(x,productivity_specif);
plt.xlabel('lotes')
plt.title('Productividad Específica')
plt.show();


plt.plot(x,t_reaccion);
plt.xlabel('lotes');
plt.ylabel('t')
plt.title('Tiempo de activación para búsqueda de energía activación')
plt.show();

plt.plot(x[1:nlotes],product[1:nlotes]);
plt.xlabel('lotes');
plt.ylabel('producto')
plt.title('Producto')
plt.show();

