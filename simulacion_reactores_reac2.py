import numpy as np;
import csv;
import scipy;
from scipy import stats;
import matplotlib.pyplot as plt;
from scipy.optimize import fsolve, root, broyden1;
from sympy import Symbol,nsolve;

import reactorbioquim;

km=10;
kcat=1;
kd=0.03;
e0=18;
Si=200;
Xopt=0.9;
nlotes=6;
beta=0.95;
kr=0.1;
V0=1;
Mcat=1;
dbeta=0.1;

e0_st, Gamma, e_ref_actividad, t_reaccion, t_reactivacion, product, productivity_specif, Productividad_especifica_global = reactorbioquim.reaccionactivada(km, kcat, kd, e0, Si, Xopt, nlotes, beta, kr, V0, Mcat, dbeta);


print('Productividad_especifica_global',Productividad_especifica_global);

print(e0_st);
print(Gamma);
print(e_ref_actividad);
print(t_reaccion);
print(t_reactivacion);
print(product);

x = np.arange(0,nlotes+1,1)
plt.plot(x,e0_st);
plt.xlabel('lotes')
plt.title('Actividad Específica')
plt.show();

plt.plot(x,e_ref_actividad);
plt.xlabel('lotes')
plt.title('Actividad Específica Reactivación')
plt.show();

plt.plot(x,productivity_specif);
plt.xlabel('lotes')
plt.title('Productividad Específica Reactivación')
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

plt.plot(x[1:nlotes],t_reactivacion[1:nlotes]);
plt.xlabel('lotes');
plt.ylabel('t')
plt.title('Tiempo de reactivación para búsqueda de energía activación deseada')
plt.show();

plt.plot(x[1:nlotes],product[1:nlotes]);
plt.xlabel('lotes');
plt.ylabel('producto')
plt.title('Producto')
plt.show();
