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
# plot results (full set of variables)
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

###########################################################################
# plot partial plot within the barplot (only for variable > 5%)
###########################################################################

######## plot fire occurrence probability
preds.dat.occ = c()
for(i in 1:ncol(Fire.1)){preds.dat.occ = cbind(preds.dat.occ, rep(mean(Fire.1[,i], na.rm = TRUE), nrow(Fire.1)))}
preds.dat.occ = data.frame(preds.dat.occ)
colnames(preds.dat.occ) <- colnames(Fire.1)
preds.dat.occ$LOshp = Fire.1$LOshp

pred.list = list()
for (i in 1:5){
idx = which(colnames(Fire.1) == a$var[i])
preds.dat.occ1 = preds.dat.occ
preds.dat.occ1[,idx] = Fire.1[,idx]
preds = predict.gbm(gbm1, preds.dat.occ1, n.trees = gbm1$gbm.call$best.trees, type="response")
pred.list[[i]] <- data.frame(preds.dat.occ1[,idx], preds)
}

######## for fire size and severity
preds.dat = c()
for(i in 1:ncol(Fire1Var.1)){preds.dat = cbind(preds.dat, rep(mean(Fire1Var.1[,i], na.rm = TRUE), nrow(Fire1Var.1)))}
preds.dat = data.frame(preds.dat)
colnames(preds.dat) <- colnames(Fire1Var.1)
preds.dat$LOshp = Fire1Var.1$LOshp

pred.list.fs = list()
for (i in 1:6){
idx = which(colnames(Fire1Var.1) == b$var[i])
preds.dat1 = preds.dat
preds.dat1[,idx] = Fire1Var.1[,idx]
preds = predict.gbm(gbm.fs, preds.dat1, n.trees = gbm.fs$gbm.call$best.trees, type="response")
pred.list.fs[[i]] <- data.frame(preds.dat1[,idx], preds)
}

pred.list.per = list()
for (i in 1:5){
idx = which(colnames(Fire1Var.1) == b2$var[i])
preds.dat1 = preds.dat
preds.dat1[,idx] = Fire1Var.1[,idx]
preds = predict.gbm(gbm.per, preds.dat1, n.trees = gbm.per$gbm.call$best.trees, type="response")
pred.list.per[[i]] <- data.frame(preds.dat1[,idx], preds)
}

### plot starts here
png(file = ".\\VariInfluence2.png", width = 2500, height = 1500, units = "px", res = 300)

par(mfrow=c(2,2),mar=c(0,3,1,0),oma=c(0,3,0,0))
# plot fire occurrece probability
a.col = merge(a, Variable.list1.c, by.x = "var", by.y = "Variable.list1")
rownames(a.col) = a.col$var
a.col = a.col[rownames(a), ]
a.col = as.character(a.col[,3])
mar=c(3,0,0,0)
barplot(height=a$rel.inf[1:5],
		horiz=FALSE,xlab="",
		#xlim = c(, 15),
		ylim = c(-4.9, 15),
		cex.axis=1.5,
		cex.lab=1.5,
		col = a.col)
#for(i in 1:nrow(a)){text(0.5+1.2*(i-1),-2.5, labels=a$var[i],cex=1.5, srt=45)}
for(i in 1:5){text(0.5+1.2*(i-1),-2.5, labels=a$var[i],cex=1.5, srt=45)}
text(0.5+1.2*(i-1),14,labels="(a)",cex=2)
text(0.7,13,labels="17.3",cex=1.5, srt=90)

xposition = c(0.08, 0.17, 0.25, 0.335, 0.415)
#plot margial effects
for (i in 1:5){
par(fig = c(xposition[i]-0.06, xposition[i]+0.06, 0.65, 0.77), new = TRUE)
mar=c(0,0,0,0)
if (i == 3){
plot(pred.list[[i]][,1], pred.list[[i]][,2],col = a.col[i],xaxt='n',yaxt='n',xlim = c(0,0.2),axes=FALSE)
lines(lowess(pred.list[[i]][,1], pred.list[[i]][,2]), lwd = 2, xlim = c(0,0.2),col = "red")
} else {
plot(pred.list[[i]][,1], pred.list[[i]][,2],  col = a.col[i],xaxt='n',yaxt='n',axes=FALSE)
lines(lowess(pred.list[[i]][,1], pred.list[[i]][,2]), lwd = 2, col = "red")
}
}

#plot fire size
par(fig = c(0.5, 1, 0.5, 1), new = TRUE)

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
#for(i in 1:nrow(b)){text(0.5+1.2*(i-1),-2.5, labels=b$var[i],cex=1.5, srt=45)}
for(i in 1:6){text(0.5+1.2*(i-1),-2.5, labels=b$var[i],cex=1.5, srt=45)}
text(0.5+1.2*(i-1),14,labels="(b)",cex=2)

xposition = c(0.075, 0.15, 0.215, 0.287, 0.35, 0.43)+0.5
#plot margial effects
for (i in 1:6){
par(fig = c(xposition[i]-0.055, xposition[i]+0.055, 0.65, 0.77), new = TRUE)
mar=c(0,0,0,0)
if(i == 2 | i == 6){
plot(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2],  col = b.col[i],xaxt='n',yaxt='n',xlim = c(0.1,1), axes=FALSE)
lines(lowess(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2]), lwd = 2, ,xlim = c(0.1,1),col = "red")
}else if (i == 3 | i == 5){
plot(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2],  col = b.col[i],xaxt='n',yaxt='n',xlim = c(0,0.2), axes=FALSE)
lines(lowess(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2]), lwd = 2, ,xlim = c(0,0.2),col = "red")
} else {
plot(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2],  col = b.col[i],xaxt='n',yaxt='n', axes=FALSE)
lines(lowess(pred.list.fs[[i]][,1], pred.list.fs[[i]][,2]), lwd = 2,,col = "red")
}
}

#plot fire severity
par(fig = c(0, 0.5, 0, 0.5), new = TRUE)
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
#for(i in 1:nrow(b2)){text(0.5+1.2*(i-1),-2.5, labels=b2$var[i],cex=1.5, srt=45)}
for(i in 1:5){text(0.5+1.2*(i-1),-2.5, labels=b2$var[i],cex=1.5, srt=45)}
text(0.5+1.2*(i-1),14,labels="(c)",cex=2)
text(0.7,13,labels="42.3",cex=1.5, srt=90)

xposition = c(0.075, 0.16, 0.245, 0.33, 0.415)
#plot margial effects
for (i in 1:5){
par(fig = c(xposition[i]-0.055, xposition[i]+0.065, 0.13, 0.25), new = TRUE)
mar=c(0,0,0,0)
if(i == 2){
plot(pred.list.per[[i]][,1], pred.list.per[[i]][,2],  col = b2.col[i],xaxt='n',yaxt='n',xlim = c(0,0.2), axes=FALSE)
lines(lowess(pred.list.per[[i]][,1], pred.list.per[[i]][,2]), lwd = 2,xlim = c(0,0.2), col = "red")

} else if (i == 3)
{
plot(pred.list.per[[i]][,1], pred.list.per[[i]][,2],  col = b2.col[i],xaxt='n',yaxt='n',xlim = c(0,20),axes=FALSE)
lines(lowess(pred.list.per[[i]][,1], pred.list.per[[i]][,2]), lwd = 2,xlim = c(0,20), col = "red")
} else {

plot(pred.list.per[[i]][,1], pred.list.per[[i]][,2],  col = b2.col[i],xaxt='n',yaxt='n',axes=FALSE)
lines(lowess(pred.list.per[[i]][,1], pred.list.per[[i]][,2]), lwd = 2, col = "red")
}
}

par(fig = c(0.5, 1, 0, 0.5), new = TRUE)
plot.new()

#add legend
legend(0.05,0.75, 
       legend = c("Vegetation", "Human", "Topography","Short-term climate" ), 
       fill = c("lightgreen", "lightcoral", "skyblue", "yellow3"), cex = 1.5, , bty = "n")

mtext("Relative importance (%)", side = 2, cex = 1.5,outer=TRUE,padj = -0.8,adj = 0.5) #http://stat.ethz.ch/R-manual/R-devel/library/graphics/html/mtext.html

dev.off()
		
