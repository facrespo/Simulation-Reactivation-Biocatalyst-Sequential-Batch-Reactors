import numpy as np;
import csv;
import scipy;
from scipy import stats;
import matplotlib.pyplot as plt;
from scipy.optimize import fsolve, root, broyden1;
from sympy import Symbol,nsolve;

import reactorbioquim;

nlotes=9;

naleatorios=5000;

km=np.random.normal(10,0,naleatorios);
kcat=np.random.uniform(1,1,naleatorios);
kd=np.random.uniform(0.021,1.307,naleatorios);
e0=np.random.normal(18,0,naleatorios);
Si=np.random.normal(200,0,naleatorios);
Xopt=0.9;
beta=np.random.normal(0.95,0,naleatorios);
kr=np.random.uniform(0.1,0.37,naleatorios);
V0=np.random.normal(1,0,naleatorios);
Mcat=np.random.normal(1,0,naleatorios);
dbeta=0.1;
xlimit=0.95;
elimit=0.25;

fich=open("resultadosimulacionactivada_n"  +str(nlotes) + ".txt","w");
fich.write("Model,km,kcat,kd,e0,Si,Xopt,nlotes,beta,kr,V0,Mcat,");
for p in range(0,nlotes+1):
    fich.write("e0_st" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("Gamma" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("e_ref_actividad" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("t_reaccion" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("t_reactivacion" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("producto" +str(p) + ",");
for p in range(0,nlotes+1):
    fich.write("productivity_specif" +str(p) + ",");
fich.write("Productividad_especifica_global");
fich.write("\n");  
fich.close();

fich2=open("resultadosimulacionsin_n"  +str(nlotes) + ".txt","w");
fich2.write("Model,km,kcat,kd,e0,Si,Xopt,nlotes,V0,Mcat,");
for p in range(0,nlotes+1):
    fich2.write("e0_st" +str(p) + ",");
for p in range(0,nlotes+1):
    fich2.write("t_reaccion" +str(p) + ",");
for p in range(0,nlotes+1):
    fich2.write("producto" +str(p) + ",");
for p in range(0,nlotes+1):
    fich2.write("productivity_specif" +str(p) + ",");
fich2.write("Productividad_especifica_global");
fich2.write("\n");  
fich2.close();



for i in range(0,naleatorios):
    e0_st, Gamma, e_ref_actividad, t_reaccion, t_reactivacion, product, productivity_specif, Productividad_especifica_global  = reactorbioquim.reaccionactivada(km[i], kcat[i], kd[i], e0[i], Si[i], Xopt, nlotes, beta[i], kr[i], V0[i], Mcat[i], dbeta, xlimit, elimit);
    fich=open("resultadosimulacionactivada_n"  +str(nlotes) + ".txt","a");
    fich.write(str(i) + "," + str(km[i]) + "," + str(kcat[i]) + "," + str(kd[i]) + "," + str(e0[i]) + "," + str(Si[i]) + "," + str(Xopt) + "," + str(nlotes) + "," + str(beta[i]) + ","  + str(kr[i]) + "," + str(V0[i]) + "," + str(Mcat[i]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(e0_st[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(Gamma[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(e_ref_actividad[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(t_reaccion[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(t_reactivacion[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(product[p]) + ",");
    for p in range(0,nlotes+1):
        fich.write(str(productivity_specif[p]) + ",");   
    fich.write(str(Productividad_especifica_global));    
    fich.write("\n");  
    fich.close();
    e0_st, t_reaccion, product, productivity_specif, Productividad_especifica_global = reactorbioquim.reaccion(km[i], kcat[i], kd[i], e0[i], Si[i], Xopt, nlotes, V0[i], Mcat[i], elimit);
    fich2=open("resultadosimulacionsin_n"  +str(nlotes) + ".txt","a");
    fich2.write(str(i) + "," + str(km[i]) + "," + str(kcat[i]) + "," + str(kd[i]) + "," + str(e0[i]) + "," + str(Si[i]) + "," + str(Xopt) + "," + str(nlotes) + "," +  str(V0[i]) + "," + str(Mcat[i]) + "," );
    for p in range(0,nlotes+1):
        fich2.write(str(e0_st[p]) + ",");
    for p in range(0,nlotes+1):
        fich2.write(str(t_reaccion[p]) + ",");
    for p in range(0,nlotes+1):
        fich2.write(str(product[p]) + ",");
    for p in range(0,nlotes+1):
        fich2.write(str(productivity_specif[p]) + ",");   
    fich2.write(str(Productividad_especifica_global));    
    fich2.write("\n");  
    fich2.close();