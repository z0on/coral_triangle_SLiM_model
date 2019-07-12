library(vegan)
library(ggplot2)

mig=as.matrix(read.csv("Amill_CRA_Mmat.csv"))
files=as.character(read.table("resfiles")[,1])

meta=do.call(rbind,strsplit(files,"[_.]"))[,3:5]
row.names(meta)=meta[,1]
meta=meta[,-1]
ids=paste(meta[,1],meta[,2],sep=".")
str(meta)
run=c("A","B","C","D","E")
rdas=list();sdss=list();gvar=list();stdms=list();pw=list();rr=list();relrr=list();ages=list()
par(mfrow=c(3,6))
for (a in 1:length(files)) { 
	load(files[a])
	rws=data.frame(raws[run[1]])
	hss=data.frame(hs[run[1]])
	for (r in run[2:length(run)]) {
		rws=rws+data.frame(raws[r])
		hss=hss+data.frame(hs[r])
	}
	rws=rws/length(run)
	hss=hss/length(run)
	gvar[[a]]=hss
	pw[[a]]=apply(rws[,1:20],1,mean)
	direction=as.numeric(pw[[a]]<apply(rws[,21:40],1,mean))
	direction[direction==0]=(-1)
	rr[[a]]=direction*apply(rws,1,function(x) { quantile(x,0.9)-quantile(x,0.1) })
	means=apply(rws,1,mean)
	sds=apply(rws,1,sd)
	stdm=(rws-means)/sds
	stdms[[a]]= stdm
	rdas[[a]]=rda(stdm~1,scale=T)
	ages[[a]]=ageall
	relrr[[a]]=rr[[a]]/pw[[a]]
	plot(rdas[[a]],main=ids[a])
}

#----- plotting genetic variation tracks 

par(mfrow=c(3,6))
for (a in 1:length(ids)){
	plot((t(gvar[[a]][100,]))~seq(-20,19,1),type="l",col="white",ylim=c(0,0.6),mgp=c(2.3,1,0), xlab="decade",ylab="SD of breeding value",main=ids[a])
	for (i in sample(1:680,680)) { lines(t(gvar[[a]][i,])~seq(-20,19,1),type="l",col=rgb(0,0,0,alpha=0.03)) }
}

#----- plotting age tracks 

plots=list()
for (a in 1:length(ids)){
	plots[[a]]=ggplot(ages[[a]],aes(x=dec,y=age,group=reef_id))+geom_line(alpha=0.01)+theme_bw()+theme(legend.position="none")+ylim(3,11)+xlab("decade")+ylab("adult age")+ggtitle(ids[a])
}

library(gridExtra)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],plots[[16]],plots[[17]],plots[[18]],nrow=3)

#------ realized heritability: 
e=0.5
for (a in 1:length(ids)){
 h2=mean(t(gvar[[a]][,19]),na.rm=T)^2/(mean(t(gvar[[a]][,19]),na.rm=T)^2+e^2)
 message(paste(ids[a],round(h2,2)))
}
	# 4.LM 0.13
	# 4.10L 0.12
	# 4.LP 0.1
	# 4.LE 0.09
	# 4.NV 0.15
	# 4.RND 0.15
	# 4.SIN 0.15
	# 4.LJM 0.15
	# 8.LM 0.13
	# 8.PS2 0.15
	# 8.10L 0.12
	# 8.LP 0.1
	# 8.LE 0.09
	# 8.NV 0.15
	# 8.RND 0.15
	# 8.SIN 0.15
	# 8.HH 0.39  << with e=0.25
	# 8.LJM 0.15

# ------- magnitudes, prewarms, and relative responses

mags=do.call(cbind, rr)[1:680,]
prs=do.call(cbind, pw)[1:680,]
rls=do.call(cbind,relrr)[1:680,]
rls[rls>1]=1.01
head(mags)
colnames(mags)=ids
colnames(rls)=ids
colnames(prs)=ids

panel.1to1=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(0,1, 
            col = "red")
}

#ordering
rlso=rls[,rev(order(apply(rls,2,mean)))]

# decline magnitudess
rlsa=abs(rlso)
rlsa[rlso>0]=NA

meta=meta[rev(order(apply(rls,2,mean))),]
prso= prs[,rev(order(apply(rls,2,mean)))]

color.gradient <- function(x, data.range=c(quantile(x,0.1),quantile(x,0.9)),colors=c("red","yellow","blue"), colsteps=100) {
		x[x<data.range[1]]= data.range[1]
		x[x>data.range[2]]= data.range[2]
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

colorScaleBar=function(x,width=4,cexx=0.5,step=0.05,colors=c("blue","yellow","red")) {
	xx=seq(min(x),max(x),length.out=100)
	yy=rep(0,100)
	plot(yy~xx,col=color.gradient(xx,colors=colors),pch=15,cex=cexx,bty="n",yaxt="n",ylab="")
	yy=yy-step
	points(yy~xx,col=color.gradient(xx,colors=colors),pch=15,cex= cexx,bty="n")
	yy=yy-step
	for (i in 1:width) {
		points(yy~xx,col=color.gradient(xx,colors=colors),pch=15,cex= cexx,bty="n")
		yy=yy-step
	}
}
#install.packages("viridis")
library(viridis)

load("sst_modified.RData")
names(sst)


# big pairs plot of response ranges
sel=c(3,which(meta[,1]=="8" & !(meta[,2] %in% c("RND","NV"))))
sel=c(which(meta[,1]=="4" & !(meta[,2] %in% c("RND","NV"))),11)
pdf(file="pairs_relResp_allGenetics_declines.pdf",height=11,width=11)
#pairs(rlsa[,sel],col=rgb(0,0,0,alpha=0.25),pch=16,cex=0.8,upper.panel=panel.1to1,lower.panel=panel.1to1)
pairs(rlsa[,sel],col=color.gradient(sst$meanT,colors=viridis(100)),pch=16,cex=0.8,upper.panel=panel.1to1,lower.panel=panel.1to1)
pairs(rlsa[,sel],col=color.gradient(sst$income.temp,colors=c("blue","yellow","red")),pch=16,cex=0.8,upper.panel=panel.1to1,lower.panel=panel.1to1)
dev.off()
colorScaleBar(sst$meanT,colors=viridis(100),cexx=0.3,step=0.03,width=20)


# correcting prewarms cover for 2x popsize setting
prso[,which(meta[,2]=="PS2")]=prso[,which(meta[,2]=="PS2")]/2
pdf(file="pairs_prewarms_allGenetics.pdf",height=11,width=11)
pairs(prso[,sel],col=color.gradient(sst$income.temp,colors=viridis(100)),pch=16,cex=0.5,upper.panel=panel.1to1,lower.panel=panel.1to1)
dev.off()
colorScaleBar(sst$income.temp,cexx=0.3,step=0.03,width=20)


# big pairs plot of response ranges
sel=which(meta[,2] %in% c("RND","NV","SIN"))
pdf(file="pairs_relResp_typeEnv_abs.pdf",height=8,width=8)
pairs(abs(rlso[,sel]),col=rgb(0,0,0,alpha=0.25),pch=16,cex=0.8,upper.panel=panel.1to1,lower.panel=panel.1to1)
#pairs(rlso[,sel],col=rgb(0,0,0,alpha=0.25),pch=16,cex=0.8,upper.panel=panel.1to1,lower.panel=panel.1to1)
dev.off()


#---------- income T and relresp
library(ggplot2)
load("sst_modified.RData")
names(sst)
sel=c(which(meta[,1]=="4" & meta[,2]=="SIN"))
sst$rr=rls[,sel]
sst$pw=prs[,sel]
sst=sst[1:680,]
summary(MASS::stepAIC(lm(rr~income.temp+meanT+income+logarea,sst)))


std=function(x) { return((x-mean(x))/sd(x)) }
sst$mTiT=std(sst$meanT)*std(sst$income.temp)

names(sst)
ggplot(sst,aes(x=meanT,y=rr))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("relative\nresponse")+xlab("")
ggplot(sst,aes(x=loginhots,y=rr))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("relative\nresponse")+xlab("")
ggplot(sst,aes(x=inhotp,y=rr))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("relative\nresponse")+xlab("")
ggplot(sst,aes(x=mTiT,y=rr))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("mean T : Tin")+xlab("")+xlim(-5,1)

ggplot(sst,aes(x=meanT,y=income.temp))+geom_point(alpha=0.2,cex=0.8)+theme_bw()
summary(lm(income.temp~meanT,sst))

sst$li=log(sst$inhots+0.001,10)
ggplot(sst,aes(x=meanT,y=loginhots))+geom_point(alpha=0.2,cex=0.8)+theme_bw()
ggplot(sst,aes(x=meanT,y=li))+geom_point(alpha=0.2,cex=0.8)+theme_bw()
ggplot(sst,aes(x=meanT,y=inhots))+geom_point(alpha=0.2,cex=0.8)+theme_bw()
summary(lm(loginhots~meanT,sst))
summary(lm(loginhots~poly(meanT,2),sst))


sst$it=sst$income.temp+rnorm(nrow(sst),0,0.3)
summary(lm(rr~income.temp,sst))
plot(income.temp~it,sst)
ggplot(sst,aes(x=meanT,y=pw))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("pre-warming\ncoral cover")+xlab("")
ggplot(sst,aes(x=mTiT,y=pw))+geom_point(alpha=0.2,cex=0.8)+theme_bw()+ylab("pre-warming\ncoral cover")+xlab("")



lm(rr~income.temp*meanT*income*logarea,sst)


summary(lm(rr~income.temp,sst))
# 4.5 : b= 0.49 R2 0.32
# 8.5 : b= 0.6  R2 0.30
summary(lm(rr~meanT,sst))
# 4.5 : b= -0.05  R2 0.11
# 8.5 : b= -0.14  R2 0.59

summary(lm(rr~income.temp*meanT,sst))
# 4.5 :  R2 0.34
# 8.5 :  R2 0.63
summary(lm(rr~loginhots*meanT,sst))
# 4.5 :  R2 0.38
# 8.5 :  R2 0.71

summary(aov(rr~loginhots*meanT,sst))

summary(aov(rr~income.temp+meanT+income,sst))
# 4.5
                   # Df Sum Sq Mean Sq F value   Pr(>F)    
# income.temp         1 13.502  13.502 330.997  < 2e-16 ***
# meanT               1  0.121   0.121   2.965 0.085548 .  
# income.temp:meanT   1  0.524   0.524  12.847 0.000362 ***
# Residuals         676 27.576   0.041                     

                   # Df Sum Sq Mean Sq F value Pr(>F)    
# income.temp         1 20.007  20.007 556.783 <2e-16 ***
# meanT               1 21.880  21.880 608.898 <2e-16 ***
# income.temp:meanT   1  0.185   0.185   5.152 0.0235 *  
# Residuals         676 24.291   0.036                   

ssq=summary(aov(rr~loginhots*meanT,sst))[[1]][,"Sum Sq"]
varexp=ssq/sum(ssq)
varexp
# 4: [1] 0.3452657049 0.0004326677 0.0392024146 0.6150992128
# 8: [1] 0.3719270 0.2788354 0.0649375 0.2843001

summary(lm(rr~mTiT,sst))


#------------ envfit

ll=load("~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/sst_mig.RData")
#save(sst,mig,mig2,file="sst_mig.RData")

# calculating genetic thermal tolerance of next generation (incl immigrants) relative to local temp

mig=as.matrix(read.csv("~/Dropbox/Documents/coral_triangle_SLiM/Amill_CRA_Mmat.csv"))
colnames(sst)
sst$income=1-diag(mig)
sums=apply(mig,1,sum)
intemp=mig;inmigs=mig;inhotp=c()
for (i in 1:680) {
	intemp[,i]=mig[,i]*sst$MEAN_meanT
	inmigs[,i]=mig[,i]*sst$km2[i]
	warmer=(sst$MEAN_meanT>sst$MEAN_meanT[i])
	table(warmer)
inhotp=append(inhotp,sum(mig[warmer,i]*sst$MEAN_meanT[warmer]))
}

hist(inhotp)
inhots=c();prophots=c();hot=0
for (i in 1:680) {
	loct=sst$meanT[i]
	hots=which(sst$meanT-loct > hot)
	message(length(hots))
	inhots[i]=sum(inmigs[hots,i])
	prophots[i]=inhots[i]/sst$km2[i]
}
hist(inhots)
sst$inhots=inhots
sst$inhotp=inhotp
summary(inhots)
hist(log(inhots+0.001,10))
sst$loginhots=log(inhots+0.001,10)
sst$inhot1=as.numeric(inhots>=1)
inmigs=sst$income*sst$km2
plot(inhots~inmigs,log="xy")
abline(0,1,col="red")


plot(inc~sst$income)
mts=apply(intemp,2,sum)
plot(mts~sst$meanT)
mtsd=mts-sst$meanT
plot(mtsd~sst$meanT)
plot(mtsd~sst$income)
sst$income.temp=mtsd
plot(mtsd~sst$income)
names(sst)=gsub("MEAN_","",names(sst))
save(sst,file="sst_modified.RData")

load("../../sst_modified.RData")

library(ggplot2)
str(sst)

a=which(meta$R=="45" & meta$V=="sin" & meta$H=="low" & meta$D=="hi" & meta$M=="hiPL")
sst$pc1=rots[[a]][1:680,"PC1"]
sst$pw=pw[[a]][1:680]
sst$pc1bi=as.numeric(sst$pc1>0)
sst$prs=prs[1:680,a]
sst$aTin=abs(sst$income.temp)
ggplot(sst,aes(x=loginhots,y=prs))+geom_point(alpha=0.2,size=1)+theme_bw()+xlab("Tin")+ylab("pre-warming\ncoral cover")
ggplot(sst,aes(x=aTin,y=prs))+geom_point(alpha=0.25)+theme_bw()
ggplot(sst,aes(x=loginhots,y=pc1))+geom_point(alpha=0.25)+theme_bw()
ggplot(sst,aes(x=income.temp,y=pc1))+geom_point(alpha=0.25)+theme_bw()
ggplot(sst,aes(x=loginhots,y=pc1))+geom_point(alpha=0.25)+theme_bw()
ggsave(paste(a,"_pc1_vs_incomeTemp.pdf",sep=""))
ggplot(sst,aes(x=income.temp,y=pw))+geom_point(alpha=0.25)+theme_bw()
ggsave(paste(a,"_prewarm_vs_incomeTemp.pdf",sep=""))
ggplot(sst,aes(y=income.temp,x=meanT))+geom_point(alpha=0.25)+theme_bw()
ggsave(paste(a,"_incomeTemp_vs_meanT.pdf",sep=""))
ggplot(sst,aes(x=meanT,y=pc1))+geom_point(alpha=0.3)+theme_bw() # +geom_smooth(method="lm")
ggsave(paste(a,"_pc1_vs_meanT.pdf",sep=""))

library(MASS)

a=which(meta$R=="45" & meta$V=="rand" & meta$H=="low" & meta$D=="hi" & meta$M=="hi")
sst$pc1=rots[[a]][1:680,"PC1"]
sst$pw=pw[[a]][1:680]
sst$pc1bi=as.numeric(sst$pc1>0)
sst$pc1=scales::rescale(rots[[a]][1:680,"PC1"])
plot(pc1~income.temp,sst,pch=16, cex=0.5,col=rgb(0,0,0,alpha=0.3),mgp=c(2.3,1,0),xlab="income T")
g=glm(pc1bi~income.temp,family="binomial",sst)
curve(scales::rescale(predict(g,data.frame(income.temp=x),type="resp"),to=c(min(sst$pc1),max(sst$pc1))),add=T)
abline(v=0,lty=3)
a=which(meta$R=="85" & meta$V=="rand" & meta$H=="low" & meta$D=="hi" & meta$M=="hi")
sst$pc1=rots[[a]][1:680,"PC1"]
sst$pw=pw[[a]][1:680]
sst$pc1bi=as.numeric(sst$pc1>0)
sst$pc1=scales::rescale(rots[[a]][1:680,"PC1"])
points(pc1~income.temp,sst,pch=16, cex=0.5,col=rgb(1,0,0,alpha=0.3))
g=glm(pc1bi~income.temp,family="binomial",sst)
curve(scales::rescale(predict(g,data.frame(income.temp=x),type="resp"),to=c(min(sst$pc1),max(sst$pc1))),add=T,col="red")

curve(predict(g,data.frame(bodysize=x),type="resp"))
summary(glm(pc1bi~meanT+income.temp+DT_RCP45,family="binomial",sst))
install.packages("popbio")
library(popbio)
logi.hist.plot(sst$income.temp, sst$pc1bi,type="hist")

stepAIC(glm(pc1bi~meanT+income.temp+DT_RCP45,family="binomial",sst))
stepAIC(lm(pc1~meanT+income.temp+DT_RCP45+logarea,sst))

summary(lm(pc1~income.temp,sst))
summary(lm(pc1~inhot1,sst))
summary(lm(pc1~loginhots,sst))

mean(sst$income.temp)
sst$itc=as.factor(sst$income.temp>0.06)
summary(lm(pc1~itc+meanT,sst))
ggplot(sst,aes(x=itc,y=pc1))+geom_boxplot()+theme_bw()+geom_vline(xintercept=0.062,lty=3,colour="grey80",lwd=1)
summary(lm(pc1~meanT,sst))

summary(lm(income.temp~income,sst))
            # Estimate Std. Error t value Pr(>|t|)  
# (Intercept) 0.001432   0.020000   0.072   0.9429  
# income      0.062842   0.034896   1.801   0.0722 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2841 on 678 degrees of freedom
# Multiple R-squared:  0.00476,	Adjusted R-squared:  0.003292 


str(rots)
for(i in 1:length(rots)) {  colnames(rots[[i]])=c("PC1","PC2") } 

a=which(meta$R=="45" & meta$V=="sin" & meta$H=="low" & meta$D=="hi" & meta$M=="hi10L")

#ef0=envfit(rots[["45sin"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP45","income")],perm=3999)
ef0=envfit(rots[[a]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP45","income")],perm=3999)
ef0
                 # PC1      PC2     r2  Pr(>r)    
# meanT        0.99137  0.13109 0.1087 0.00025 ***
# income.temp -0.96225 -0.27217 0.2001 0.00025 ***
# loginhots   -0.99997  0.00802 0.2822 0.00025 ***
# inhot1      -0.99967 -0.02552 0.2587 0.00025 ***
# logarea      0.96451 -0.26403 0.0147 0.00825 ** 
# DT_RCP45     0.32202 -0.94673 0.0253 0.00050 ***
# income      -0.73684  0.67607 0.1414 0.00025 ***
#ef00=envfit(rots[["452sin"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP45","income")],perm=3999)
ef00=envfit(rots[["452sin"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP45","income")],perm=3999)
ef00
                 # PC1      PC1     r2  Pr(>r)    
# meanT        0.93823 -0.34602 0.1584 0.00025 ***
# income.temp -0.99875 -0.04997 0.2296 0.00025 ***
# loginhots   -0.99662  0.08209 0.2657 0.00025 ***
# inhot1      -0.99976 -0.02182 0.2085 0.00025 ***
# logarea      0.99554  0.09439 0.0032 0.34600    
# DT_RCP45     0.22830 -0.97359 0.0076 0.07950 .  
# income      -0.78022  0.62550 0.0677 0.00025 ***

#ef1=envfit(rots[["85sin"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP85","income")],perm=3999)
ef1=envfit(rots[["85sin"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP85","income")],perm=3999)
ef1
                 # PC1      PC2     r2  Pr(>r)    
# meanT        0.87843  0.47787 0.1590 0.00025 ***
# income.temp -0.97387 -0.22711 0.2249 0.00025 ***
# loginhots   -0.99803 -0.06272 0.1977 0.00025 ***
# inhot1      -0.99895 -0.04575 0.1522 0.00025 ***
# logarea      0.94671  0.32209 0.0156 0.00375 ** 
# DT_RCP85     0.72520 -0.68854 0.0327 0.00025 ***
# income      -0.74513  0.66692 0.1386 0.00025 ***
#ef11=envfit(rots[["852sin"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP85","income")],perm=3999)
ef11=envfit(rots[["852sin"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP85","income")],perm=3999)
ef11
                 # PC1      PC1     r2  Pr(>r)    
# meanT        0.99160 -0.12932 0.1764 0.00025 ***
# income.temp -0.99998  0.00670 0.2754 0.00025 ***
# loginhots   -0.88873 -0.45843 0.1352 0.00025 ***
# inhot1      -0.84510 -0.53461 0.0992 0.00025 ***
# logarea      0.96916 -0.24642 0.0124 0.01600 *  
# DT_RCP85     0.64690  0.76258 0.0089 0.05350 .  
# income      -0.99910  0.04241 0.0420 0.00025 ***


#ef02=envfit(rots[["45sin2"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP45","income")],perm=3999)
ef02=envfit(rots[["45sin2"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP45","income")],perm=3999)
ef02

                 # PC1      PC2     r2  Pr(>r)    
# meanT        0.88740 -0.46100 0.1392 0.00025 ***
# income.temp -0.99118 -0.13249 0.2185 0.00025 ***
# loginhots   -0.91350  0.40683 0.3280 0.00025 ***
# inhot1      -0.90885  0.41713 0.2985 0.00025 ***
# logarea      0.96514 -0.26175 0.0117 0.01850 *  
# DT_RCP45     0.13047 -0.99145 0.0707 0.00025 ***
# income      -0.78231  0.62289 0.0995 0.00025 ***

#ef2=envfit(rots[["85sin2"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP85","income")],perm=3999)
ef2=envfit(rots[["85sin2"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP85","income")],perm=3999)
ef2
                 # PC1      PC2     r2  Pr(>r)    
# meanT        0.75304 -0.65797 0.3640 0.00025 ***
# income.temp -0.99254  0.12193 0.3193 0.00025 ***
# loginhots   -0.99646 -0.08401 0.2056 0.00025 ***
# inhot1      -0.99217 -0.12493 0.1641 0.00025 ***
# logarea      0.91328 -0.40733 0.0203 0.00100 ***
# DT_RCP85     0.98761 -0.15692 0.0165 0.00400 ** 
# income      -0.99984  0.01762 0.0598 0.00025 ***
#ef22=envfit(rots[["852sin"]][1:680,],sst[,c("meanT","income.temp","loginhots","inhot1","logarea","DT_RCP85","income")],perm=3999)
ef22=envfit(rots[["8522sin"]][1:680,],sst[,c("meanT","income.temp","logarea","DT_RCP85","income")],perm=3999)
ef22
                 # PC1      PC2     r2  Pr(>r)    
# meanT        0.99160 -0.12932 0.1764 0.00025 ***
# income.temp -0.99998  0.00670 0.2754 0.00025 ***
# loginhots   -0.88873 -0.45843 0.1352 0.00025 ***
# inhot1      -0.84510 -0.53461 0.0992 0.00025 ***
# logarea      0.96916 -0.24642 0.0124 0.01250 *  
# DT_RCP85     0.64690  0.76258 0.0089 0.04625 *  
# income      -0.99910  0.04241 0.0420 0.00025 ***

# how much variation is explained by meanT and income.temp?
for (a in files) {
	sst$pc1=rots[[a]][1:680,"PC1"]
	message(paste(a,"\t",round(summary(lm(pc1~meanT+income.temp,sst))$adj.r.squared,3),sep=""))
}
# 45	0.249
# 45sin	0.201
# 45nv	0.224
# 45sin2	0.228
# 85	0.279
# 85sin	0.232
# 85nv	0.268
# 85sin2	0.355


angle0 <- function(M,N=c(1,0)){
  aa=acos( sum(M*N) / ( sqrt(sum(M*M)) * sqrt(sum(N*N)) ) )
  if (M[2]<0) { aa=2*pi-aa }
  len=sqrt(sum(M*M))
  return (c(aa,len))
}

a="452sin";ef=ef00
datt=data.frame(rots[[a]])
names(datt)=c("PC1","PC2")
offset=mean(datt[c(683,687),1])
datt[,1]=datt[,1]-offset
d2=c()
for(i in 1:680) { d2=data.frame(rbind(d2,angle0(datt[i,]))) }
d2[,1]=d2[,1]*4/(2*pi)

plot(d2[,2]~rls[1:680,a])
resp.type=d2[1:680,1]
sst$resp.type=d2[1:680,1]
ggplot(sst,aes(x=income.temp,y=resp.type))+geom_point(alpha=0.25)+theme_bw()

label_offset=0.17 # offset for environment labels
scale.env.arrows=2
scale.datapoints=12
color_arrows="red"
color_sites="red"
color_points=rgb(0,0,0,alpha=0.3)

ef.df<-as.data.frame(ef$vectors$arrows*scale.env.arrows*sqrt(ef$vectors$r))
ef.df2<-as.data.frame(ef$vectors$arrows*scale.env.arrows*sqrt(ef$vectors$r+label_offset^2))
ef.df2$species<-rownames(ef.df2)
dats=data.frame(apply(datt[1:680,],2,function(x) { x*scale.datapoints}))
str(dats)

plot(sst$pc1~d2[,2])
ggplot()+
	geom_point(data=dats, aes(PC1,PC2,colour=resp.type),alpha=0.3)+
#	geom_point(data=dats, aes(PC1,PC2,colour=d2[,2]),alpha=0.3)+
	scale_colour_gradientn(colours=c("plum","cyan3","green3","coral","plum"))+
	geom_segment(data=ef.df[ef$vectors$pvals<0.05,],aes(x=0,xend=PC1,y=0,yend=PC2),arrow = arrow(length = unit(0.2, "cm")),colour=color_arrows)+
#	geom_text(data=ef.df2[ef$vectors$pvals<0.05,],aes(PC1,PC2,label=species),col=color_arrows,size=3)+
	ggtitle(a)+
	theme(legend.position="none")+
	guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)+
	theme_void()
ggsave(paste("pca_envfit_",a,"_noarrowlabels.pdf",sep=""),device="pdf")
#	geom_point(data=data.frame(dats[681:688,]), aes(PC1,PC2),pch=1,size=5,col="grey70")


#-------------- map

ll=load("sst_modified.RData")

angle0 <- function(M,N=c(1,0)){
  aa=acos( sum(M*N) / ( sqrt(sum(M*M)) * sqrt(sum(N*N)) ) )
  if (M[2]<0) { aa=2*pi-aa }
  len=sqrt(sum(M*M))
  return (c(aa,len))
}

a=which(meta$R=="85" & meta$run=="PS2")
#plot(rots[[a]])
dat=rots[[a]]
dat[,1]=dat[,1]-mean(dat[c(683,687),1]) 
d2=c()
for(i in 1:nrow(dat)) { d2=data.frame(rbind(d2,angle0(dat[i,]))) }
d2[,1]=d2[,1]*4/(2*pi)
ssr=sst
ssr$resp.type=d2[1:680,1]
ssr$pc1=dat[1:680,1]
ssr$magnitude=mags[1:680,a]
ssr$relmag=rls[1:680,a]
ssr=ssr[order(ssr$magnitude),]
ssr$income.temp[ssr$income.temp>1]=1
ssr$income.temp[ssr$income.temp<(-1)]=(-1)

ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=income.temp,cex=km2), alpha=0.5)+scale_colour_gradientn(colours=c("cyan3","coral")) + theme_void()+ggtitle(ids[a])

ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=meanT,cex=km2), alpha=0.7)+scale_colour_gradientn(colours=c("cyan3","coral")) + theme_void()+ggtitle(ids[a])

ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=resp.type,cex=relmag), alpha=0.7)+scale_colour_gradientn(colours=c("plum","cyan3","green3","coral","plum")) + theme_void()+ggtitle(ids[a])
# +scale_size(limits=c(0,1),range=c(1,5))
#ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=resp.type,cex=magnitude), alpha=0.7)+scale_colour_gradientn(colours=c("plum","cyan3","green3","coral","plum")) + theme_void()+ggtitle(a)
ggsave(paste("map_",ids[a],".pdf",sep=""),device="pdf")



