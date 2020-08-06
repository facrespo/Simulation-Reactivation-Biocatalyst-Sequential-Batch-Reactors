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
nlotes=8;
beta=0.95;
kr=0.1;
V0=1;
Mcat=1;
dbeta=0.1;
xlimit=0.95;
elimit=0.25;

Xopt=0.9;

e0_st, betaf, Gamma, e_ref_actividad, t_reaccion, t_reactivacion, product, Productividad, Productividadr = reactorbioquim.reaccionactivada(km, kcat, kd, e0, Si, Xopt, nlotes, beta, kr, V0, Mcat, dbeta, xlimit, elimit);

print('Producto',product);
print('Productividad reacción',Productividad);
print('Productividad_especifica_global',Productividadr);

print(e0_st);
print(Gamma);
print(betaf);
print(e_ref_actividad);
print(t_reaccion);
print(t_reactivacion);


e1 = e0;
curva_without_activation=np.zeros((10,2*(nlotes+1)),dtype=float);
for i in range(0,nlotes):
    if (i==0):
        curva_without_activation[0:10,0:4] = reactorbioquim.curva_activada(km, kcat, kd, e0, e0_st[i], e_ref_actividad[i+1], Si, betaf[i+1], kr, V0, Mcat, Gamma[i], i+1, xlimit);
    else:
        temp=reactorbioquim.curva_activada(km, kcat, kd, e0, e0_st[i], e_ref_actividad[i+1], Si, betaf[i+1], kr, V0, Mcat, Gamma[i+1], i+1, xlimit);
        curva_without_activation[0:10,(2*i+2):(2*i+4)] =temp[0:10,2:4];


fich2=open("plot_curve_with_activation.txt","w");
fich2.write("Point,Xopt,");
for p in range(0,(nlotes)):
    fich2.write(",,lotes_" +str(p+1) + ", ,");
fich2.write("\n");
for i in range(0,10):
   for p in range(0,2*(nlotes)):
      fich2.write(str(curva_without_activation[i,p]) + ",");
      if (p>=3) and ((p%2)):
        fich2.write(" ,");  
   fich2.write("\n");  

fich2.close();

plt.plot(curva_without_activation[0:9,2],curva_without_activation[0:9,1]);
plt.xlabel('time');
plt.ylabel('Activity')
plt.title('Plot of Activity verus time')
plt.show();