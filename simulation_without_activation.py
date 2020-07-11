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
kd=0.022;
e0=18;
Si=200;
Xopt=0.1;
nlotes=6;
kr=0.5;
V0=1;
Mcat=1;
elimit=0.25

e0_st, t_reaccion, product, productivity_specif, Productividad_especifica_global = reactorbioquim.reaccion(km, kcat, kd, e0, Si, Xopt, nlotes, V0, Mcat, elimit);

print('Producto',product);
print('Productividad reacción',productivity_specif);
print('Productividad_especifica_global',Productividad_especifica_global);

print(e0_st);
print(t_reaccion);
print(product);

x = np.arange(0,nlotes+1,1)
plt.plot(x,e0_st);
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

e1 = e0_st[0];
curva_without_activation=np.zeros((10,nlotes+2),dtype=float);
for i in range(0,nlotes):
    if (i==0):
        curva_without_activation[0:10,0:3] = reactorbioquim.curva(km, kcat, kd, e1, Si, V0, Mcat);
        e1=e1*np.exp(-kd*curva_without_activation[9,2]);
    else:
        if (e1 >= elimit*e0_st[0]):
           temp=reactorbioquim.curva(km, kcat, kd, e1, Si, V0, Mcat);
           curva_without_activation[0:10,i+2] =temp[0:10,2];
           e1=e1*np.exp(-kd*temp[9,2]);

fich2=open("plot_curve_whithout_activation.txt","w");
fich2.write("Point,Xopt,");
for p in range(0,nlotes):
    fich2.write("lotes_" +str(p+1) + ",");
fich2.write("\n");
for i in range(0,10):
   for p in range(0,nlotes+2):
      fich2.write(str(curva_without_activation[i,p]) + ",");
   fich2.write("\n");  

fich2.close();

plt.plot(curva_without_activation[0:9,2],curva_without_activation[0:9,1]);
plt.xlabel('time');
plt.ylabel('Activity')
plt.title('Plot of Activity verus time')
plt.show();
