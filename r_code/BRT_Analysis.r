###########################################################################
# Boosted Regression Tree analysis
###########################################################################

setwd(".\\data")

library(gbm)
library(dismo)

###########################################################################
# for case-control fire occurrence analysis 
###########################################################################
Fire.1 = data.frame(read.csv("Fire.Occ.csv")[,-1])
Fire.1$LOshp = as.factor(Fire.1$LOshp)
Variable.list1.exclu = c("H90Mn","T90Mn","HPAo")
Fire.1 = Fire.1[ , -which(names(Fire.1) %in% Variable.list1.exclu)]

# select best BRT settings based on deviance and AUC
tc = c(1,2,3); lr = c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1); nt = c(10, 20, 30, 40,50)
model.perf = c()						
for (i in 1:length(tc)){
 for (j in 1:length(lr)){
  for (k in 1:length(nt)){
  gbm1=gbm.step(data = Fire.1, gbm.x = 1:(ncol(Fire.1)-2),
						gbm.y = ncol(Fire.1)-1, 
						family = "bernoulli",
						tree.complexity = tc[i], 
						n.trees = nt[k], learning.rate = lr[j], bag.fraction = 0.75)
						
preds = predict.gbm(gbm1, Fire.1, n.trees = gbm1$gbm.call$best.trees, type="response")
dev = calc.deviance(obs=Fire.1$Status, pred=preds, calc.mean=TRUE)
	
d <- cbind(Fire.1$Status, preds)
pres <- d[d[,1]==1, 2]
abse <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abse)
model.perf.t = c(tc[i],lr[j],nt[k],dev, e@auc, gbm1$gbm.call$best.trees)
model.perf = rbind(model.perf, model.perf.t)

print(paste("Finish Fitting tree complexity = ", tc[i], "; learning.rate = ", lr[j]," n.trees = ", nt[k], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  }
 }
}
write.csv(model.perf, "model.perf.occ.csv")

#best BRT setting: tree.complexity = 3; n.trees = 40; learning.rate = 0.075,with best tree names at 2720;  AUC = 0.97

gbm1=gbm.step(data = Fire.1, gbm.x = 1:(ncol(Fire.1)-2),
						gbm.y = ncol(Fire.1)-1, 
						family = "bernoulli",
						tree.complexity = 3, 
						n.trees = 40, learning.rate = 0.075, bag.fraction = 0.75)
###########################################################################
# for fire size and severity analysis
###########################################################################
Fire1Var.1 = as.data.frame(read.csv("Fire.SizeSeverity.csv")[,-1])
Fire1Var.1$LOshp = factor(Fire1Var.1$LOshp)
Fire1Var.1$Calacre_ha = log(Fire1Var.1$Calacre_ha)
Variable.list2.exclu = c("T90Mn","H90Ao","HPAo","HWAo")
Fire1Var.1 = Fire1Var.1[ , -which(names(Fire1Var.1) %in% Variable.list2.exclu)]

# select best BRT settings for fire size model
model.perf = c()						
for (i in 1:length(tc)){
 for (j in 1:length(lr)){
  for (k in 1:length(nt)){
gbm.fs=gbm.step(data = Fire1Var.1, gbm.x = 1:(ncol(Fire1Var.1)-4),
						gbm.y = ncol(Fire1Var.1)-3, 
						family = "gaussian",
						tree.complexity = tc[i], 
						n.trees = nt[k], learning.rate = lr[j], bag.fraction = 0.75)
						
r <- Fire1Var.1$Calacre_ha - predict.gbm(gbm.fs, Fire1Var.1,n.trees = gbm.fs$gbm.call$best.trees, type="response")
dev = 1-var(r)/var(Fire1Var.1$Calacre_ha) 
	
model.perf.t = c(tc[i],lr[j],nt[k],dev, gbm.fs$gbm.call$best.trees)
model.perf = rbind(model.perf, model.perf.t)

print(paste("Finish Fitting tree complexity = ", tc[i], "; learning.rate = ", lr[j]," n.trees = ", nt[k], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  }
 }
}
write.csv(model.perf, "model.perf.size.csv")
# best setting for fire size model, tree.complexity = 3; n.trees = 40; learning.rate = 0.05,with best tree names at 1960;  variance explained = 0.76

gbm.fs=gbm.step(data = Fire1Var.1, gbm.x = 1:(ncol(Fire1Var.1)-4),
						gbm.y = ncol(Fire1Var.1)-3, 
						family = "gaussian",
						tree.complexity = 3, 
						n.trees = 25, learning.rate = 0.05, bag.fraction = 0.75)
# select best BRT settings for fire severity model
model.perf = c()						
for (i in 1:length(tc)){
 for (j in 1:length(lr)){
  for (k in 1:length(nt)){
gbm.p=gbm.step(data = Fire1Var.1, gbm.x = 1:(ncol(Fire1Var.1)-4),
						gbm.y = ncol(Fire1Var.1)-1, 
						family = "gaussian",
						tree.complexity = tc[i], 
						n.trees = nt[k], learning.rate = lr[j], bag.fraction = 0.75)
						
r <- Fire1Var.1$FSH_PRO - predict.gbm(gbm.p, Fire1Var.1,n.trees = gbm.p$gbm.call$best.trees, type="response")
dev = 1-var(r)/var(Fire1Var.1$FSH_PRO) 
	

model.perf.t = c(tc[i],lr[j],nt[k],dev, gbm.p$gbm.call$best.trees)
model.perf = rbind(model.perf, model.perf.t)

print(paste("Finish Fitting tree complexity = ", tc[i], "; learning.rate = ", lr[j]," n.trees = ", nt[k], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  }
 }
}
write.csv(model.perf, "model.perf.severity.csv")

# best setting for percent high severity model, tree.complexity = 3; n.trees = 20; learning.rate = 0.1,with best tree names at 600;  variance explained = 0.60
gbm.per=gbm.step(data = Fire1Var.1, gbm.x = 1:(ncol(Fire1Var.1)-3),
						gbm.y = ncol(Fire1Var.1), 
						family = "gaussian",
						tree.complexity = 3, 
						n.trees = 20, learning.rate = 0.1, bag.fraction = 0.75)


###########################################################################
# plot results
###########################################################################
best.iter=gbm.perf(gbm1,method="OOB")
a=summary(gbm1,n.trees=best.iter,cBars=5,las=1,cex.axis=1.5,cex.lab=1.6,mar=c(5,8,4,2)+0.1,axes=FALSE,axisnames=FALSE, plotit=F)
				
best.iter=gbm.perf(gbm.fs,method="OOB")
b = summary(gbm.fs,n.trees=best.iter,cBars=5,las=1,cex.axis=1.5,cex.lab=1.6,mar=c(5,8,4,2)+0.1,axes=FALSE,axisnames=FALSE, plotit=F)

best.iter=gbm.perf(gbm.per,method="OOB")
b2 = summary(gbm.per,n.trees=best.iter,cBars=5,las=1,cex.axis=1.5,cex.lab=1.6,mar=c(5,8,4,2)+0.1,axes=FALSE,axisnames=FALSE, plotit=F)

#set color scheme
Variable.list1.c = c(rep("lightgreen", 13), rep("lightcoral", 2), rep("skyblue", 3), rep("yellow3", 13), "lightcoral", "lightcoral" , "lightcoral")
Variable.list1.c = cbind(Variable.list1, Variable.list1.c)

Variable.list2.c = c(rep("lightgreen", 13), rep("lightcoral", 2), rep("skyblue", 3), rep("yellow3", 17), "lightcoral", "lightcoral" , "lightcoral", "lightcoral")
Variable.list2.c = cbind(Variable.list2, Variable.list2.c)

#plot
png(file = ".data\\VariInfluence.png", width = 2000, height = 2000, units = "px", res = 300)
par(mfrow=c(2,2),mar=c(0,3,1,0),oma=c(0,3,0,0))

a.col = merge(a, Variable.list1.c, by.x = "var", by.y = "Variable.list1")
rownames(a.col) = a.col$var
a.col = a.col[rownames(a), ]
a.col = as.character(a.col[,3])
mar=c(3,0,0,0)
barplot(height=a$rel.inf[1:5],horiz=FALSE,xlab="",
		#xlim = c(, 15),
		ylim = c(-4.9, 15),
		cex.axis=1.5,
		cex.lab=1.5,
		col = a.col)
#for(i in 1:nrow(a)){text(0.5+1.2*(i-1),-2.5, labels=a$var[i],cex=1.5, srt=45)}
for(i in 1:5){text(0.5+1.2*(i-1),-2.5, labels=a$var[i],cex=1.5, srt=45)}

text(0.5+1.2*(i-1),14,labels="(a)",cex=2)
text(0.7,13,labels="17.3",cex=1.5, srt=90)
a.sign = c("-", "-", "+", "+", "-", "+", "", "+","","","","","","","","","+", "","-","","+","+")
for(i in 1:20){text(0.7+1.2*(i-1),2.5,labels=a.sign[i],cex=2)}

b.col = merge(b, Variable.list2.c, by.x = "var", by.y = "Variable.list2")
rownames(b.col) = b.col$var
b.col = b.col[rownames(b), ]
b.col = as.character(b.col[,3])
mar=c(3,0,0,0)
barplot(height=b$rel.inf[1:6],horiz=FALSE,xlab="",
		#xlim = c(, 15),
		ylim = c(-4.9, 15),
		cex.axis=1.5,
		cex.lab=1.5,
		col = b.col)
for(i in 1:6){text(0.5+1.2*(i-1),-2.5, labels=b$var[i],cex=1.5, srt=45)}
text(0.5+1.2*(i-1),14,labels="(b)",cex=2)
b1.sign = c("+", "-","+", "+", "+", "-", "","","","","","","","-","","-","","-","-","+" )
for(i in 1:20){text(0.7+1.2*(i-1),2.5,labels=b1.sign[i],cex=2)}

b2.col = merge(b2, Variable.list2.c, by.x = "var", by.y = "Variable.list2")
rownames(b2.col) = b2.col$var
b2.col = b2.col[rownames(b2), ]
b2.col = as.character(b2.col[,3])
mar=c(3,0,0,0)
barplot(height=b2$rel.inf[1:5],horiz=FALSE,xlab="",
		#xlim = c(, 15),
		ylim = c(-4.9, 15),
		cex.axis=1.5,
		cex.lab=1.5,
		col = b2.col)
for(i in 1:5){text(0.5+1.2*(i-1),-2.5, labels=b2$var[i],cex=1.5, srt=45)}
text(0.5+1.2*(i-1),14,labels="(c)",cex=2)
text(0.7,13,labels="42.3",cex=1.5, srt=90)
b2.sign = c("+", "+", "-", "+", "+", "","","","-","","","+","+","","","+","+")
for(i in 1:20){text(0.7+1.2*(i-1),2.5,labels=b2.sign[i],cex=2)}

plot.new()

#add legend
legend(0.05,0.75, 
       legend = c("Vegetation", "Human", "Topography","Short-term climate" ), 
       fill = c("lightgreen", "lightcoral", "skyblue", "yellow3"), cex = 1.5, , bty = "n")

mtext("Relative importance (%)", side = 2, cex = 1.5,outer=TRUE,padj = -0.8,adj = 0.5) 

dev.off()

	