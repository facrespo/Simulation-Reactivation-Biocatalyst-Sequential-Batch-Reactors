import numpy as np;
import scipy;
from scipy import stats;
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

def solveactivity(Eopt, beta, gamma, kr, e0):
    #x0=0;
    x0=0;
    resolucion=scipy.optimize.fsolve(eactividad,x0,args=(beta, gamma, kr, e0, Eopt))
    #resolucion=-(1/kr)*np.log((Eopt-beta*e0)/(e0*(gamma-beta)));
    return resolucion;

def reactivity(km, kcat, kd, e0, Si, Xopt, nlotes, beta, kr, V0, Mcat, dbeta, xlimit, elimit):
    efinish=np.zeros(nlotes+1,dtype=float);
    betat=np.zeros(nlotes+1,dtype=float);
    e0_activity=np.zeros(nlotes+1,dtype=float);
    Gamma=np.zeros(nlotes+1,dtype=float);
    t_reaccion=np.zeros(nlotes+1,dtype=float);
    t_reactivacion=np.zeros(nlotes+1,dtype=float);
    V=np.zeros(nlotes+1,dtype=float);
    X=np.zeros(nlotes+1,dtype=float);
    product=np.zeros(nlotes+1,dtype=float);
    productivity_specif=np.zeros(nlotes+1,dtype=float);
    betat[0]=0;
    betat[1]=1;
    betat[2]=beta;
    efinish[0]=e0;
    e0_activity[0]=e0;
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
        if (i>=3) and (efinish[i-1] >= elimit*efinish[0]):
           betat[i]=betat[i-1]-dbeta;
           if (betat[i]<=0):
               betat[i]=0;
        if (efinish[i-1] >= elimit*efinish[0]):
           if (i==1):
               e0_activity[i]=efinish[0];
           else:
               e0_activity[i]=xlimit*efinish[0]*betat[i];
           if (i==1):
               t_reaccion[i] = solvereactor(Xopt, Si, km, kcat, kd, e0_activity[i-1]);
           else:
               t_reaccion[i] = solvereactor(Xopt, Si, km, kcat, kd, e0_activity[i]);
           if t_reaccion[i]>0:
              efinish[i]=e0_activity[i]*np.exp(-kd*t_reaccion[i]);
              Gamma[i]=efinish[i]/efinish[0];
              if (i==1):
                  t_reactivacion[i]=0;
              else:
                  t_reactivacion[i]=solveactivity(e0_activity[i], betat[i], Gamma[i], kr, efinish[0]);
              tiempo=tiempo+ t_reaccion[i]+t_reactivacion[i];
              X[i]=Xopt;
              V[i]=V0;
              product[i]=(((Si*V[i]*X[i])/Mcat));
              productivity_specif[i]=product[i]/(t_reaccion[i]+t_reactivacion[i]);
       
    Productividad_especifica_global=sum(product)/(tiempo);
    return efinish, betat, Gamma, e0_activity, t_reaccion, t_reactivacion, product, productivity_specif, Productividad_especifica_global;

def inactivity(km, kcat, kd, e0, Si, Xopt, nlotes, V0, Mcat, elimit):
    efinish=np.zeros(nlotes+1,dtype=float);
    e_ref_actividad=np.zeros(nlotes+1,dtype=float);
    t_reaccion=np.zeros(nlotes+1,dtype=float);
    V=np.zeros(nlotes+1,dtype=float);
    X=np.zeros(nlotes+1,dtype=float);
    product=np.zeros(nlotes+1,dtype=float);
    productivity_specif=np.zeros(nlotes+1,dtype=float);
    efinish[0]=e0;
    t_reaccion[0]=0;
    V[0]=V0;
    X[0]=1;
    product[0]=0;
    tiempo=0;
    Productividad=0;
    Productividadr=0;
    for i in range(1,nlotes+1):
        en=e0;
        if (i>=1) and (efinish[i-1] >= elimit*efinish[0]):
           t_reaccion[i] = solvereactor(Xopt, Si, km, kcat, kd, efinish[i-1]);
           if t_reaccion[i]>0:
               efinish[i]=efinish[i-1]*np.exp(-kd*t_reaccion[i]);
               tiempo=tiempo+ t_reaccion[i];
               X[i]=Xopt;
               V[i]=V0;
               #product[i-1]+
               product[i]=(((Si*V[i]*X[i])/Mcat));
               productivity_specif[i]=product[i]/t_reaccion[i];
           else:
               e0_st[i]=0;
               tiempo=tiempo;
               X[i]=0;
               V[i]=0;
               product[i]=0;
               productivity_specif[i]=0;
       
    Productividad_especifica_global=sum(product)/(tiempo);
    return efinish, t_reaccion, product, productivity_specif, Productividad_especifica_global;

def curv(km, kcat, kd, e0, Si, V0, Mcat):
    Xopt=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], dtype=float);
    proyeccion=np.zeros((10,3),dtype=float);
    proyeccion[0,0]=1;
    for i in range(1,proyeccion.shape[0]):
        proyeccion[i,0]=proyeccion[0,0]+i;
        proyeccion[i,1]=Xopt[i-1];
        proyeccion[i,2] = solvereactor(Xopt[i-1], Si, km, kcat, kd, e0);
                   
    return proyeccion;

def curv_activity(km, kcat, kd, e0, e0f, eref, Si, beta, kr, V0, Mcat, gamma, nlote, xlimit):
    Xopt=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], dtype=float);
    proyeccion=np.zeros((10,4),dtype=float);
    proyeccion[0,0]=1;
    for i in range(1,proyeccion.shape[0]):
        t2=0;
        proyeccion[i,0]=proyeccion[0,0]+i;
        proyeccion[i,1]=Xopt[i-1];
        if (nlote==1):
           t1 = solvereactor(Xopt[i-1], Si, km, kcat, kd, e0);
        else:
           t1 = solvereactor(Xopt[i-1], Si, km, kcat, kd, eref);
        if (i==9):
           if (nlote==1):
               t2=0;
           else:
               t2=solveactivity(eref, beta, gamma, kr, e0);
           if (t2<=0):
               t2=0;
        proyeccion[i,2]=t1;
        proyeccion[i,3]=t2;
    
    return proyeccion;



