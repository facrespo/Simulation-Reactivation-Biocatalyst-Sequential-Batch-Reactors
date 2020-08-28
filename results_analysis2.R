setwd("E:/investigacion/nadia_guajardo");

library(stats);
library(ggplot2);
library(Rmisc);
library(ggpubr);
library(psych);
library(fields);
library(scatterplot3d);
library(mvnormtest);
library(tidyverse)
library(ggpubr)
library(rstatix)
library(car)
library(broom)

ntotalr=9;

interleave <- function(l, how = c('cbind', 'rbind')) {
  how <- match.arg(how)
  if (how %in% 'rbind')
    do.call(how, l)[order(sequence(sapply(l, nrow))), ]
  else do.call(how, l)[, order(sequence(sapply(l, ncol))), ]
}


for (i in ntotalr:ntotalr){
  eval(parse(text=noquote(paste('sin', as.character(i), '<- read.table(file=\"resultadosimulacionwithoutactivate_n', as.character(i), '.txt\", header = T, dec = \".\", sep = \",\");', sep=""))));
  eval(parse(text=noquote(paste('con', as.character(i), '<- read.table(file=\"resultadosimulacionactivated_n', as.character(i), '.txt\", header = T, dec = \".\", sep = \",\");', sep=""))));
}

plot(sin9$kd,sin9$Productividad_especifica_global,xlab = "k_d", ylab = "Global_specific_productivity", main="Case without reactivation");

scatterplot3d(x=con9$kd, y=con9$kr,z=con9$Productividad_especifica_global,  xlab="k_d", ylab="k_r",  zlab = "Global_specific_productivity", main="Case With reactivation");

mm=max(sin9$Productividad_especifica_global);
i=which.max(sin9$Productividad_especifica_global);
sin9[i,]$kd

mm=max(con9$Productividad_especifica_global);
i=which.max(con9$Productividad_especifica_global);
con9[i,]



ndim1=dim(con9);
ndim2=dim(sin9);

ntot <- 1;
datos <- as.data.frame(matrix(0,2*(ntot)*ndim1[1],12));

l=1;
for (i in ntotalr:ntotalr){
  for (j in 1:ndim1[1]){
    datos[l,1]="with reactivation";
    datos[l,2]=i;
    datos[l,3:12]=con9[j,74:83];
    l=l+1;
  }
  for (j in 1:ndim2[1]){
    datos[l,1]="without reactivation";
    datos[l,2]=i;
    datos[l,3:12]=sin9[j,42:51];
    l=l+1;
  }
}

colnames(datos)[3:12] <- colnames(con9[,74:83]);
colnames(datos)[1:2] <- c("Type","Nbatch");
colnames(datos)[12] <- c("Global_specific_productivity");

datos = as.data.frame(datos);

#str(datos)

#datos$Nbatch = as.numeric(as.character(datos$Nbatch));
#datos$Productividad_especifica_global = as.numeric(as.character(datos$Productividad_especifica_global));
#datos$Producto = as.numeric(as.character(datos$Producto));

summary(datos);

#Productividad

sum = summarySE(datos, 
                measurevar="Global_specific_productivity", 
                groupvars=c("Type","Nbatch"));
sum

colnames(sum)[4] <- c("Global_specific_productivity");

sum

write.table(sum,file="resumen_productividad_global.txt",quote=FALSE,row.names=FALSE,sep=';');


pd = position_dodge(.1);

ggplot(sum, aes(x=Nbatch, 
                y=Global_specific_productivity, 
                color=Type)) + 
  geom_errorbar(data=sum,aes(ymin=Global_specific_productivity-se, 
                    ymax=Global_specific_productivity+se), 
                width=.9, size=0.7, position=pd) +
  geom_point(shape=15, size=4, position=pd) +
  theme_bw() +
  theme(
    axis.title.y = element_text(vjust= 1.8), xlab = "Batch",
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  scale_color_manual(values = c("black", "blue"));

ggboxplot(datos, x = "Nbatch", y = "Global_specific_productivity", 
          color = "Type", ylab = "Global_specific_productivity", xlab = "Batch");

boxplot(Global_specific_productivity ~ Type,
        data = datos,
        xlab = "Type",
        ylab = "Global specific productivity");


model = lm(Global_specific_productivity ~ Type, data=datos);

Anova(model, type="II");

anova(model);

summary(model);

hist(residuals(model), col="darkgray");

plot(fitted(model),residuals(model));

res.aov <- aov(Productividad_especifica_global ~ Type, data = datos);

resumen <- (summary(res.aov));

resumen_global <- describeBy(datos[,12],datos$Type);

resumen_global <- interleave(resumen_global,'r');

write.table(resumen_global,file="analysis_anova_global_productivity.txt",quote=FALSE,sep=';');


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





summary(datos[,4:6])



measurevar=colnames(datos[,3:8]);

testshapiro1 <- datos[1:5000,] %>%
  select(productivity_specif1,productivity_specif2,productivity_specif3,productivity_specif4,productivity_specif5,productivity_specif6) %>%
  mshapiro_test()

testshapiro2 <- datos[5001:10000,] %>%
  select(productivity_specif1,productivity_specif2,productivity_specif3,productivity_specif4) %>%
  mshapiro_test()

res.manova <- manova(cbind(productivity_specif1,productivity_specif2,productivity_specif3,productivity_specif4,productivity_specif5,productivity_specif6)
                  ~ Type, data = datos);

shapiro1 <- datos %>%
              group_by(Type) %>%
              shapiro_test(productivity_specif1, productivity_specif2,productivity_specif3,productivity_specif4)

shapiro2 <- datos %>%
  shapiro_test(productivity_specif1, productivity_specif2,productivity_specif3,productivity_specif4)




summary(res.manova)

resumenprod <- describeBy(datos[,3:11],datos$Type);

productivity <- interleave(resumenprod,'r');

write.table(productivity[1:12,],file="analysis_productivity.txt",quote=FALSE,sep=';');


summary.aov(res.manova);

ggboxplot(
  datos, x = "Type", y = c("productivity_specif1", "productivity_specif2", 
                           "productivity_specif3", "productivity_specif4",
                           "productivity_specif5", "productivity_specif6"), 
  merge = TRUE, palette = "grey"
)

