setwd("E:/investigacion/nadia_guajardo");

library(stats);
library(ggplot2);
library(Rmisc);
library(ggpubr);
library(psych);

ntotalr=10;

for (i in 5:ntotalr){
  eval(parse(text=noquote(paste('sin', as.character(i), '<- read.table(file=\"resultadosimulacionsin_n', as.character(i), '.txt\", header = T, dec = \".\", sep = \",\");', sep=""))));
  eval(parse(text=noquote(paste('con', as.character(i), '<- read.table(file=\"resultadosimulacionactivada_n', as.character(i), '.txt\", header = T, dec = \".\", sep = \",\");', sep=""))));
}

ndim=dim(con5);


datos <- matrix(0,2*(ntotalr-4)*ndim[1],4);

l=1;
for (i in 5:ntotalr){
  for (j in 1:ndim[1]){
    datos[l,1]="reactivation";
    datos[l,2]=i;
    datos[l,3]=(eval(parse(text=noquote(paste('con', as.character(i), '$Productividad_especifica_global[ ', as.character(j) ,']', sep="")))));
    datos[l,4]=(eval(parse(text=noquote(paste('con', as.character(i), '$producto',  as.character(i), '[', as.character(j), ']', sep="")))));
    l=l+1;
  }
  for (j in 1:ndim[1]){
    datos[l,1]="normal";
    datos[l,2]=i;
    datos[l,3]=(eval(parse(text=noquote(paste('sin', as.character(i), '$Productividad_especifica_global[ ', as.character(j) ,']', sep="")))));
    datos[l,4]=(eval(parse(text=noquote(paste('sin', as.character(i), '$producto',  as.character(i), '[', as.character(j), ']', sep="")))));
    l=l+1;
  }
}

colnames(datos) <- c("Type","Nbatch","Productividad_especifica_global","Producto");

datos = as.data.frame(datos);

str(datos)

datos$Nbatch = as.numeric(as.character(datos$Nbatch));
datos$Productividad_especifica_global = as.numeric(as.character(datos$Productividad_especifica_global));
datos$Producto = as.numeric(as.character(datos$Producto));

summary(datos);

#Productividad

sum = summarySE(datos, 
                measurevar="Productividad_especifica_global", 
                groupvars=c("Type","Nbatch"));
sum

pd = position_dodge(.1);

ggplot(sum, aes(x=Nbatch, 
                y=Productividad_especifica_global, 
                color=Type)) + 
  geom_errorbar(data=sum,aes(ymin=Productividad_especifica_global-se, 
                    ymax=Productividad_especifica_global+se), 
                width=.9, size=0.7, position=pd) +
  geom_point(shape=15, size=4, position=pd) +
  theme_bw() +
  theme(
    axis.title.y = element_text(vjust= 1.8),
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  scale_color_manual(values = c("black", "blue"));

ggboxplot(datos, x = "Nbatch", y = "Productividad_especifica_global", 
          color = "Type", ylab = "Productividad_especifica_global", xlab = "Batch");

boxplot(Productividad_especifica_global ~ Nbatch:Type,
        data = datos,
        xlab = "Batch x Type",
        ylab = "Productividad_especifica_global");


model = lm(Productividad_especifica_global ~ Type + Nbatch + Type:Nbatch, data=datos);

Anova(model, type="II");

anova(model);

summary(model);

hist(residuals(model), col="darkgray");

plot(fitted(model),residuals(model));

res.aov <- aov(Productividad_especifica_global ~ Type:Nbatch, data = datos);

summary(res.aov);

TukeyHSD(res.aov);

plot(res.aov, 1);

leveneTest(Productividad_especifica_global ~ Type:Nbatch, data = datos);

oneway.test(Productividad_especifica_global ~ Type:Nbatch, data = datos);

plot(res.aov, 2);

aov_residuals <- residuals(object = res.aov );

shapiro.test(x = aov_residuals );

kruskal.test(Productividad_especifica_global ~ Type, data = datos);

kruskal.test(Productividad_especifica_global ~ Nbatch, data = datos);

describeBy(datos$Productividad_especifica_global,datos$Type);

describeBy(datos$Productividad_especifica_global,datos$Type);



-#Producto
  
sum = summarySE(datos, 
                  measurevar="Producto", 
                  groupvars=c("Type","Nbatch"));
sum

pd = position_dodge(.1);

ggplot(sum, aes(x=Nbatch, 
                y=Producto, 
                color=Type)) + 
  geom_errorbar(data=sum,aes(ymin=Producto-se, 
                             ymax=Producto+se), 
                width=.9, size=0.7, position=pd) +
  geom_point(shape=15, size=4, position=pd) +
  theme_bw() +
  theme(
    axis.title.y = element_text(vjust= 1.8),
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  scale_color_manual(values = c("black", "blue"));

ggboxplot(datos, x = "Nbatch", y = "Producto", 
          color = "Type", ylab = "Producto", xlab = "Batch");

boxplot(Producto ~ Nbatch:Type,
        data = datos,
        xlab = "Batch x Type",
        ylab = "Producto");


model = lm(Producto ~ Type + Nbatch + Type:Nbatch, data=datos);
           
Anova(model, type="II");
           
anova(model);
           
summary(model);
        
hist(residuals(model), col="darkgray");
           
plot(fitted(model),residuals(model));
           
res.aov <- aov(Producto ~ Type:Nbatch, data = datos);
           
summary(res.aov);
           
TukeyHSD(res.aov);
           
plot(res.aov, 1);
           
leveneTest(Producto ~ Type:Nbatch, data = datos);
           
oneway.test(Producto ~ Type:Nbatch, data = datos);

by(datos$Producto ~ datos$Type:Nbatch, shapiro.test);
           
plot(res.aov, 2);
           
aov_residuals <- residuals(object = res.aov );
           
shapiro.test(x = aov_residuals );
           
kruskal.test(Producto ~ Type, data = datos);
           
kruskal.test(Producto ~ Nbatch, data = datos);
           

#Manova

res.man <- manova(cbind(Productividad_especifica_global, Producto) ~ Type:Nbatch, data = datos);

summary(res.man);

summary.aov(res.man);
