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


e1 = e0;
curva_with_activation=np.zeros((10,2*(nlotes+1)),dtype=float);
for i in range(0,nlotes):
    if (i==0):
        curva_with_activation[0:10,0:4] = reactorbioquim.curv_activity(K, k, kd, e0, et_activity[i], e0_activity[i+1], Si, betaf[i+1], kr, V, Mbiocat, Gamma[i], i+1, xlimit);
    else:
        temp=reactorbioquim.curv_activity(K, k, kd, e0,  et_activity[i], e0_activity[i+1], Si, betaf[i+1], kr, V, Mbiocat, Gamma[i+1], i+1, xlimit);
        curva_with_activation[0:10,(2*i+2):(2*i+4)] =temp[0:10,2:4];


fich2=open("plot_curve_with_activation.txt","w");
fich2.write("Point,Xopt,");
for p in range(0,(nlotes)):
    fich2.write(",,lotes_" +str(p+1) + ",");
fich2.write("\n");
for i in range(0,10):
   for p in range(0,2*(nlotes)):
      fich2.write(str(curva_with_activation[i,p]) + ",");
      if (p>=3) and ((p%2)):
        fich2.write(" ,");  
   fich2.write("\n");  

fich2.close();

plt.plot(curva_with_activation[0:9,2],curva_with_activation[0:9,1]);
plt.xlabel('time');
plt.ylabel('Activity')
plt.title('Plot of Activity verus time')
plt.show();