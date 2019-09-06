import numpy as np;
import csv;
import scipy;
from scipy import stats;
import matplotlib.pyplot as plt;
from scipy.optimize import fsolve, root, broyden1;
from sympy import Symbol,nsolve;


def reactor(t, Xopt, Si, km, kcat, kd, e0):
    final=np.zeros(1,dtype=float)
    final[0]=np.matrix(((Si/km)*Xopt-np.log(1-Xopt))-(((kcat*e0)/(km*kd))*(1-np.exp(-kd*t))));
    return final;

def actividad(t, beta, gamma, kr, e0):
    finale=np.zeros(1,dtype=float)
    finale[0]=e0*((gamma-beta)*(np.exp(-kr*t))+beta);
    return finale;

def eactividad(t, beta, gamma, kr, e0, en):
    finale=np.zeros(1,dtype=float)
    finale[0]=actividad(t, beta, gamma, kr, e0)-en;
    return finale;

def eqreactor(observacion, data):
    t = observacion;
    X=data[0];
    Si=data[1];
    km=data[2];
    kcat=data[3];
    kd=data[4];
    e0=data[5];
    return reactor(t, X, Si, km, kcat, kd, e0);

def solvereactor(Xopt, Si, km, kcat, kd, e0):
    x0=0;
    resolucion=scipy.optimize.fsolve(reactor,x0,args=(Xopt, Si, km, kcat, kd, e0))
    return resolucion;

def solveactividad(Eopt, beta, gamma, kr, e0):
    x0=0;
    resolucion=scipy.optimize.fsolve(eactividad,x0,args=(beta, gamma, kr, e0, Eopt))
    return resolucion;

def reaccionactivada(km, kcat, kd, e0, Si, Xopt, nlotes, beta, kr, V0, Mcat, dbeta):
    e0_st=np.zeros(nlotes+1,dtype=float);
    betat=np.zeros(nlotes+1,dtype=float);
    e_ref_actividad=np.zeros(nlotes+1,dtype=float);
    Gamma=np.zeros(nlotes+1,dtype=float);
    t_reaccion=np.zeros(nlotes+1,dtype=float);
    t_reactivacion=np.zeros(nlotes+1,dtype=float);
    V=np.zeros(nlotes+1,dtype=float);
    X=np.zeros(nlotes+1,dtype=float);
    product=np.zeros(nlotes+1,dtype=float);
    productivity_specif=np.zeros(nlotes+1,dtype=float);
    betat[0]=beta;
    betat[1]=beta;
    e0_st[0]=e0;
    e_ref_actividad[0]=e0;
    Gamma[0]=1;
    t_reactivacion[0]=0;
    t_reaccion[0]=0;
    V[0]=V0;
    X[0]=1;
    product[0]=0;
    tiempo=0;
    Productividad=0;
    Productividadr=0;
    for i in range(1,nlotes+1):
        en=e0;
        if (i>=2):
           betat[i]=betat[i-1]-dbeta;
           if (betat[i]<=0):
               betat[i]=0;
        if i>=1:
           t_reaccion[i] = solvereactor(Xopt, Si, km, kcat, kd, e_ref_actividad[i-1]);
           e0_st[i]=e_ref_actividad[i-1]*np.exp(-kd*t_reaccion[i]);
           Gamma[i]=e0_st[i]/e_ref_actividad[i-1];
           Eopt=e0*(0.95)*(0.95-0.1*(i-1));
           e_ref_actividad[i]=Eopt;
           t_reactivacion[i]=solveactividad(e_ref_actividad[i], betat[i], Gamma[i], kr, e0);
           tiempo=tiempo+ t_reaccion[i]+t_reactivacion[i-1];
           X[i]=Xopt;
           V[i]=V0;
           product[i]=product[i-1]+(((Si*V[i]*X[i])/Mcat));
           productivity_specif[i]=product[i]/tiempo;
       
    Productividad_especifica_global=sum(product)/(tiempo);
    return e0_st, Gamma, e_ref_actividad, t_reaccion, t_reactivacion, product, productivity_specif, Productividad_especifica_global;

def reaccion(km, kcat, kd, e0, Si, Xopt, nlotes, V0, Mcat):
    e0_st=np.zeros(nlotes+1,dtype=float);
    e_ref_actividad=np.zeros(nlotes+1,dtype=float);
    t_reaccion=np.zeros(nlotes+1,dtype=float);
    V=np.zeros(nlotes+1,dtype=float);
    X=np.zeros(nlotes+1,dtype=float);
    product=np.zeros(nlotes+1,dtype=float);
    productivity_specif=np.zeros(nlotes+1,dtype=float);
    e0_st[0]=e0;
    t_reaccion[0]=0;
    V[0]=V0;
    X[0]=1;
    product[0]=0;
    tiempo=0;
    Productividad=0;
    Productividadr=0;
    for i in range(1,nlotes+1):
        en=e0;
        if i>=1:
           t_reaccion[i] = solvereactor(Xopt, Si, km, kcat, kd, e0_st[i-1]);
           if t_reaccion[i]>0:
               e0_st[i]=e0_st[i-1]*np.exp(-kd*t_reaccion[i]);
               tiempo=tiempo+ t_reaccion[i];
               X[i]=Xopt;
               V[i]=V0;
               product[i]=product[i-1]+(((Si*V[i]*X[i])/Mcat));
               productivity_specif[i]=product[i]/tiempo;
           else:
               e0_st[i]=0;
               tiempo=tiempo;
               X[i]=0;
               V[i]=0;
               product[i]=product[i-1];
               productivity_specif[i]=0;
       
    Productividad_especifica_global=sum(product)/(tiempo);
    return e0_st, t_reaccion, product, productivity_specif, Productividad_especifica_global;
