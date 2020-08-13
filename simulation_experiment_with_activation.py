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
nlotes=8;
beta=0.95;
kr=0.1;
V=1;
Mbiocat=1;
dbeta=0.1;
xlimit=0.95;
elimit=0.25;


Xopt=0.9;

et_activity, betaf, Gamma, e0_activity, t_reaccion, t_reactivacion, product, Productividad, Productividadr = reactorbioquim.reactivity(K, k, kd, e0, Si, Xopt, nlotes, beta, kr, V, Mbiocat, dbeta, xlimit, elimit);

print('Producto',product);
print('Productividad reacción',Productividad);
print('Productividad_especifica_global',Productividadr);

print(Gamma);
print(betaf);
print(e0_activity);
print(et_activity);
print(t_reaccion);
print(t_reactivacion);

x = np.arange(0,nlotes+1,1)
plt.plot(x,et_activity);
plt.xlabel('lotes')
plt.title('Actividad Específica en t final')
plt.show();

plt.plot(x,e0_activity);
plt.xlabel('lotes')
plt.title('Actividad Específica Reactivación en t=0')
plt.show();

plt.plot(x,Gamma);
plt.xlabel('lotes')
plt.title('Gamma')
plt.show();


plt.plot(x,t_reaccion);
plt.xlabel('lotes');
plt.ylabel('t')
plt.title('Tiempo de activación para búsqueda de energía activación')
plt.show();

plt.plot(x[0:nlotes],t_reactivacion[0:nlotes]);
plt.xlabel('lotes');
plt.ylabel('t')
plt.title('Tiempo de reactivación para búsqueda de energía activación deseada')
plt.show();

