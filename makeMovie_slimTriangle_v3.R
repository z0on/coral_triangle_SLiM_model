#------------

library(tidyverse)
library(Rmisc)
library(vegan)
library(pheatmap)
library(viridis)
library(gridExtra)

setwd("~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/slim_triangle_results_apr2019")
load("~/Dropbox/Documents/coral_triangle_SLiM/Amill_CRA_OutMmat_sst_modified_may2019.RData")

#install.packages('gpclib', type='source')
library(gpclib)
if (!rgeosStatus()) gpclibPermit()

mapfile="../../HiRes_shoreline/gshhs_l.b"
sf1 <- getRgshhsMap(mapfile, xlim = c(min(sst$lon)-1, max(sst$lon)+1), ylim = c(min(sst$lat)-5, max(sst$lat)+3)) 

files=c("triRE85_3A.clean" ,"triRE853A.clean" ,"triRE852sinA.clean","triRE8522sinA.clean")
#files=c("triRE453A.clean" ,"triRE853A.clean")

names(sst)
readslim=function(filename) {
	g85=read.table(filename,sep="\t",skip=19)
	names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
	g85$reef_id=g85$reef_id+1
	g85$reef_id=as.factor(g85$reef_id)
	g85=subset(g85,g>5300 & g<=5700)
	g85=merge(g85,sst[,c(1:5,35:37)],by="reef_id",all.x=T)
	gg=g85 %>% mutate(cover=nad/km2)
	gg=gg %>% arrange(g)
	return(gg)
}

ress=list()
for (f in files) {
	message(f)
	ress[[f]]=readslim(f) 
	}

rescaler = function(x, to = c(0, 1), from = NULL) {
    			ifelse(x<1.2, 
    			scales::rescale(x, to = to, from = c(min(x, na.rm = TRUE), 1.2)), 
    			1)
}

i=5451
str(sst)
rcp=c("realRfSiz","TinyJuv","8.5","8.5")
mu=c("1e-5","1e-5","1e-5","1e-6")


mxx=c();mnn=c();minn=c();i=1
for (i in 1:length(files)) { mxx=append(mxx,max(ress[[i]]$cover[1:20])) }
for (i in 1:length(files)) { mnn=append(mnn,mean(ress[[i]]$cover[1:20])) }
for (i in 1:length(files)) { minn=append(minn,min(ress[[i]]$cover[1:60])) }
mnn
mxx
topcover= max(mxx)
for (i in 1:length(files)) {
#	ress[[i]]$cover[ress[[i]]$cover>topcover]=topcover
	ress[[i]]$cover=rescaler(ress[[i]]$cover)
}

hist(ress[[4]]$cover)

i=5640
library(viridis)
library(gridExtra)

for (i in c(5480:5700)){
	plo=list()
	txtcol="grey90"
	if(i>5500) { txtcol="yellow"}
	for (f in 1:length(files)) {
		sgg=subset(ress[[files[f]]],g==i)
		plo[[f]]=ggplot() +
			geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey60', color='black', lwd = 0.1) +
			geom_point(data=sgg, aes(lon,lat,color=cover,cex=km2), alpha=0.5) +
			scale_colour_viridis(option="B",limits=c(0,topcover))+
			theme_void()+
			annotate("text",
				label=paste("RCP ",rcp[f],"; mu ",mu[f],sep=""),
				col="grey20",cex=5,x=max(sgg$lon)-13,y=max(sgg$lat))+
			theme(legend.position = "none",
				plot.margin=unit(c(-0.05,-0.05,-0.1,-0.1), "cm"))
#			theme(panel.background = element_rect(fill = "grey40",colour = "grey40",size = 0.5, linetype = "solid"))

		
		if (f==1) { 	plo[[f]]=plo[[f]]+annotate("text",
				label=i-5500,
				col=txtcol,cex=7,x=min(sgg$lon)+5,y=max(sgg$lat)-1)
		}
	}
	jpeg(filename=paste("jpg_jj/",i,".jpg",sep=""),width=640,height=750)
	grid.arrange(plo[[1]],plo[[2]],plo[[3]],plo[[4]],ncol=2,heights=c(1,1))
#	grid.arrange(plo[[1]],plo[[2]],ncol=2,heights=c(1))
	dev.off()
}
	

# ress2=ress
# ress[[2]]=ress2[[3]]
# ress[[3]]=ress2[[2]]	
i=5500
plo=list()
txtcol="skyblue"
if(i>5500) { txtcol="grey10"}
for (f in 1:length(files)) {
	sgg=subset(ress[[files[f]]],g==i)
	plo[[f]]=ggplot() +
			    geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey90', color='grey40', lwd = 0.1) +
			geom_point(data=sgg, aes(lon,lat,color=cover,cex=km2), alpha=0.5) +
			scale_colour_viridis(option="B",limits=c(mnn[f],mxx[f]))+
			theme_void()+
			# annotate("text",
				# label=paste("RCP ",rcp[f],"; mu ",mu[f],sep=""),
				# col="grey90",x=max(sgg$lon)-12,y=max(sgg$lat))+
			theme(legend.position = "none",
				plot.margin=unit(c(-0.05,-0.05,-0.1,-0.1), "cm"))
#			theme(panel.background = element_rect(fill = "grey40",colour = "grey40",size = 0.5, linetype = "solid"))

}
load('~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/slim_triangle_results_apr2019/connectivity_map_ggplot.RData')
conmap2=conmap+	theme_void()+theme(legend.position = "none", plot.margin=unit(c(-0.05,-0.05,-0.1,-0.1), "cm"))
pdf(file=paste(i,".pdf",sep=""),width=7.5, height=9)
grid.arrange(conmap2,plo[[1]],plo[[3]],plo[[4]],ncol=2,heights=c(1,1))
dev.off()

grid.arrange(plo[[1]],plo[[2]],plo[[3]],plo[[4]],ncol=2,heights=c(1,1))


load('~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/slim_triangle_results_apr2019/connectivity_map_ggplot.RData')

names(sst)
m45=emptymap+geom_point(data=sst,aes(lon,lat,color=DT_RCP45,cex=km2),alpha=0.5)+scale_color_gradient( low="blue",high="coral")+	theme_void()
m85=emptymap+geom_point(data=sst,aes(lon,lat,color=DT_RCP85,cex=km2),alpha=0.5)+scale_color_gradient( low="blue",high="coral")+	theme_void()
mph=emptymap+geom_point(data=sst,aes(lon,lat,color=lprophots,cex=km2),alpha=0.5)+scale_color_viridis(option="B")+	theme_void()
pdf(file="prophots_RCPs.pdf",width=9.5, height=9)
grid.arrange(m45,m85,plo[[4]]+theme_void(),mph,ncol=2,heights=c(1,1))
dev.off()

# ----- variation in coral cover
str(ress[[1]])

p=1
sp=subset(ress[[3]],reef_id==p)
sp$g=sp$g-5500
lo=loess(cover~g,sp,span=0.25)
sp$Wlo=predict(lo,newdata=sp)
plot((sp$cover-sp$Wlo)/sp$cover~sp$g,type="l",col=rgb(0,0,0,alpha=0.01),ylim=c(-0.25,0.25))
for (p in 2:680){
	sp=subset(ress[[3]],reef_id==p)
	sp$g=sp$g-5500
	lo=loess(cover~g,sp,span=0.25)
	sp$Wlo=predict(lo,newdata=sp)
	lines((sp$cover-sp$Wlo)/sp$cover~sp$g,type="l",col=rgb(0,0,0,alpha=0.01))
}

#---------- SD, CV

	
readcover=function(x) {
		g85=read.table(x,sep="\t",skip=19)
		names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
		tail(g85)
		g85$reef_id=g85$reef_id+1
		g85$reef_id=as.factor(g85$reef_id)
		g85=subset(g85,g>=5300 & g<5700)
		g85=merge(g85,sst[,c(1,4)],by="reef_id",all.x=T)
		gg=g85 %>% mutate(cover=nad/km2)
		gg=gg %>% arrange(g)
		gg$dec=ceiling((gg$g-min(gg$g)+1)/10)
		cov.all=summarySE(gg, measurevar="cover",groupvars=c("reef_id","dec"))
		cov.all=cov.all %>% mutate(CV=sd/cover)
		return(cov.all)
}

library(tidyverse)
library(Rmisc)
ll=load("~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/sst_modified.RData")
c4=readcover("triRE452sinA.clean")
c8=readcover("triRE852sinA.clean")
c8nv=readcover("triRE852nvA.clean")
c8rnd=readcover("triRE852A.clean")

dat=c8
head(dat)
dat$dec=dat$dec-20
ggplot(dat,aes(dec,CV,group=reef_id))+geom_line(alpha=0.01)+theme_bw()+ylim(0,0.1)

run=c("A","B","C","D","E")
for (s in c("852hsin") ){
	respranges=list();prewarms=list();raws=list();stds=list();relresps=list();hs=list()
	
	;index=""
 	
	r="C";sets=s;index=""
	message(sets)
	for (r in run){
		message(paste("   ",r))
		g85=read.table(paste("triRE",sets,r,index,".clean",sep=""),sep="\t",skip=19)
		names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
		tail(g85)
		g85$reef_id=g85$reef_id+1
		g85$reef_id=as.factor(g85$reef_id)
		g85=subset(g85,g>=5300 & g<5700)
		g85=merge(g85,sst[,c(1,4)],by="reef_id",all.x=T)
		gg=g85 %>% mutate(cover=nad/km2)
		gg=gg %>% arrange(g)
		gg$dec=ceiling((gg$g-min(gg$g)+1)/10)
		cov.all=summarySE(gg, measurevar="cover",groupvars=c("reef_id","dec"))
		if (r=="A") { 
			covall=cov.all 
			} else { 
			covall[,4]=covall[,4]+cov.all[,4]
			covall[,5]=covall[,5]+cov.all[,5]
		}
		h.all=summarySE(gg, measurevar="h",groupvars=c("reef_id","dec"))
		age.all=summarySE(gg, measurevar="age",groupvars=c("reef_id","dec"))
		if (r=="A") { 
			ageall=age.all 
			} else {					
			ageall[,4]=ageall[,4]+age.all[,4]
			ageall[,5]=ageall[,5]+age.all[,5]
		}
		ct=cov.all %>% dplyr::select(reef_id,dec,cover) %>% spread(reef_id,cover) %>% dplyr::select(-dec)
		ht=h.all %>% dplyr::select(reef_id,dec,h) %>% spread(reef_id,h) %>% dplyr::select(-dec)
		mid=t(data.frame(cbind(ct,toons)))
		hid=t(data.frame(ht))
		hs[[r]]=hid
		raws[[r]]=mid
		means=apply(mid,1,mean)
		sds=apply(mid,1,sd)
		stds[[r]]=(mid-means)/sds
		prewarms[[r]]=apply(mid[,1:20],1,mean)
		respranges[[r]]=apply(mid,1,function(x) { quantile(x,0.9)-quantile(x,0.1) })
		direction=as.numeric(prewarms[[r]]<apply(mid[,21:40],1,mean))
		direction[direction==0]=(-1)
		table(direction[1:680])
		relresps[[r]]=direction*respranges[[r]]/prewarms[[r]]
	}
	covall$dec=covall$dec-20	
	ageall$dec=ageall$dec-20	
	covall[,4]=covall[,4]/length(run)
	covall[,5]=covall[,5]/length(run)
	ageall[,4]=ageall[,4]/length(run)
	ageall[,5]=ageall[,5]/length(run)
		
	names(raws)=run
	names(hs)=run
	names(stds)=run
	names(prewarms)=run
	names(respranges)=run
	
	rra=do.call(cbind,respranges)
	mrra=apply(rra,1,mean)
	pw=do.call(cbind,prewarms)
	mpw=apply(pw,1,mean)
	save(raws,hs,stds,prewarms, relresps,mpw,mrra,covall,ageall,file=paste("RE_sin_",sets,index,".RData",sep=""))

}


reef=10
ra=raws[["A"]][reef,]
rb=raws[["B"]][reef,]
plot(rb~ra,pch=16,cex=0.5,col=rgb(1,0,0,alpha=1),ylim=c(0,1),xlim=c(0,1))
reef=100
ra=raws[["A"]][reef,]
rb=raws[["B"]][reef,]
points(rb~ra,pch=16,cex=0.5,col=rgb(0,0,1,alpha=1))
reef=200
ra=raws[["A"]][reef,]
rb=raws[["B"]][reef,]
points(rb~ra,pch=16,cex=0.5,col=rgb(1,0,1,alpha=1))
reef=300
ra=raws[["A"]][reef,]
rb=raws[["B"]][reef,]
points(rb~ra,pch=16,cex=0.5,col=rgb(1,1,0,alpha=1))







#load("3runs_85.5.esd.sin_age.RData")

std.mid=stds[["A"]]
for (i in 1:nrow(std.mid)) { std.mid[i, std.mid[i,]>3]=3 }
colnames(std.mid)=seq(10,400,10)
pheatmap(std.mid,cluster_cols=F,clustering_distance_rows="correlation",clustering_method="average")
#pheatmap(t(sdt),cluster_cols=F,clustering_distance_rows="correlation")

load('~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/RE_sin45.RData')
rdas=list()
quartz()
par(mfrow=c(3,8))
run=c("A","B","C","D","E","F")
for (r in run){
	std.mid=stds[[r]]
	mm=rda(std.mid~1,scale=T)
	plot(mm$CA$u[,1:2],pch=16,col=rgb(0,0,0,0.25),main=paste("mean",r))
	str(mm$CA$u)
	text(mm$CA$u[681:688,1:2],row.names(mm$CA$u)[681:688],col="red")
	rdas[[r]]=mm
	std.mid=SDstds[[r]]
	mm=rda(std.mid~1,scale=T)
	plot(mm$CA$u[,1:2],pch=16,col=rgb(0,0,0,0.25),main=paste("sd",r))
	str(mm$CA$u)
	text(mm$CA$u[681:688,1:2],row.names(mm$CA$u)[681:688],col="red")
	plot(procrustes(rdas[[r]],mm))
}

toplot=670
plot(stds[[r]][toplot,]~c(1:40),ylim=c(-3,3))
lines(SDstds[[r]][toplot,]~c(1:40),col="red")

plot(SDstds[[r]][toplot,]~stds[[r]][toplot,])

cors=c()
for (r in run) {
	for (reef in 1:680) { 
		cors=append(cors,cor(SDstds[[r]][reef,],stds[[r]][reef,])) 
	}
}
Cors=data.frame(cbind(cors,"reef"=1:680))
Cors$cors=as.numeric(as.character(Cors$cors))
cc=summarySE(Cors,measurevar="cors",groupvars="reef")
str(cc)
head(cc)
plot(cors~sd,cc)

plot(density(cc$cors))
sst$corResp=cc$cors
names(sst)
ggplot(sst, aes(km2,corResp,colour=income))+geom_point(alpha=0.5)+scale_color_gradientn(colours=c("cyan3","grey","coral"))+theme_bw()+scale_x_log10()
ggplot(sst, aes(income,corResp,colour=logarea))+geom_point(alpha=0.5)+scale_color_gradientn(colours=c("cyan3","grey","coral"))+theme_bw()
summary(lm(corResp~income,sst))
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.64141    0.02320  -27.65   <2e-16 ***
# income       0.36387    0.04048    8.99   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.3296 on 678 degrees of freedom
# Multiple R-squared:  0.1065,	Adjusted R-squared:  0.1052 

summary(lm(corResp~logarea,sst))
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.83744    0.05952   14.07   <2e-16 ***
# logarea     -0.42513    0.01912  -22.24   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2651 on 678 degrees of freedom
# Multiple R-squared:  0.4217,	Adjusted R-squared:  0.4209 


par(mfrow=c(2,3))
for (r in run){ plot(rdas[[r]]) }

par(mfrow=c(2,3))
errs=c()
for (r1 in 1:(length(run)-1)){
	for (r2 in (r1+1):length(run)){
		pp=procrustes(rdas[[run[r1]]],rdas[[run[r2]]])		
		errs=append(errs,summary(pp)$rmse)
		plot(pp,main=paste(run[r1],run[r2]))
	}
}

 pcoas=list()
 par(mfrow=c(2,3))
 for (r in run){
	 std.mid=stds[[r]]
	plot(mm$CA$u[,1:2],pch=16,col=rgb(0,0,0,0.25))
	str(mm$CA$u)
	text(mm$CA$u[681:688,1:2],row.names(mm$CA$u)[681:688],col="red")
	 pcoas[[r]]=mm
 }

# quartz()
# par(mfrow=c(2,3))
# errsd=c()
# for (r1 in 1:(length(run)-1)){
	# for (r2 in (r1+1):length(run)) {
		# pp=procrustes(pcoas[[run[r1]]],pcoas[[run[r2]]])		
		# errsd=append(errsd,summary(pp)$rmse)
		# plot(pp,main=paste(run[r1],run[r2]))
	# }
# }

# plot(errs,errsd)
abline(0,1)

#-------
# mean PCA:

rws=raws[[run[1]]]
str(raws[[r]])
str(rws)
for (r in run[-1]) {
	rws=rws+raws[[r]]
}

rws=rws/length(run)

# sdr=SDstds[[run[1]]]
# for (r in run[-1]) {
	# sdr=sdr+SDstds[[r]]
# }

# sdr=sdt/length(sdr)
# pheatmap(ssds,cluster_cols=F)


means=apply(rws,1,mean)
sds=apply(rws,1,sd)
stdm=(rws-means)/sds

mm=rda(stdm~1,scale=T)
plot(mm$CA$u[,1:2],pch=16,col=rgb(0,0,0,0.25))
cv=1-abs(cc$sd/cc$cors)
cv[cv<0]=0
plot(density(cv))
plot(density(cc$cors),bty="n",mgp=c(2.3,1,0),main="",ylab="",yaxt="n",xlab="cover ~ SD(cover) correlation")
abline(v=0,lty=3)
dats=data.frame(mm$CA$u[1:680,])
dats=dats[order(cv),]
ggplot()+geom_point(data=dats, aes(PC1,PC2,colour=sst$corResp),alpha=cv[order(cv)])+scale_colour_gradientn(colours=c("cyan3","coral"))+theme_void()
ggplot()+geom_point(data=dats, aes(PC1,PC2,colour=sst$corResp),alpha=0.5)+scale_colour_gradientn(colours=c("cyan3","grey","coral"))+theme_void()

str(mm$CA$u)
text(mm$CA$u[681:688,1:2],row.names(mm$CA$u)[681:688],col="red")

#mm=capscale(vegdist(stdm,method="manhattan")~1)
plot(mm$CA$eig)

plot(mm)
plot(mm,dis="sp")
ef=envfit(mm$CA$u[1:680,],sst[,c(5,8,16,17)],perm=3999)
ef

ra45=mm
plot(mm)
plot(ef,cex=0.9,col="cyan3")
ra45=mm
plot(ra45)
plot(procrustes(ra45$CA$u[1:680,],ra85$CA$u[1:680,]))

#------ rotate ordination to match some previous analysis

rotato=function(x,rmat){ 
	return (c(x[1]*rmat[1,1]+x[2]*rmat[1,2],x[1]*rmat[2,1]+x[2]*rmat[2,2]))
}

standard=rdas[[1]]
pp=procrustes(ra85$CA$u,standard)
rotsA=t(apply(ra85$CA$u[,1:2],1,function (x) {rotato(x,pp$rotation)}))
plot(rotsA)
# plot(scores(mm)[["sites"]])
# plot(ra.sin.85)
# plot(ra.sin.45$CA$u[,1:2])

quartz()
plot(rotsA,pch=16,col=rgb(0,0,0,0.25))
text(rotsA[681:688,1:2],row.names(rotsA)[681:688],col="red")
ef=envfit(rotsA[1:680,],sst[,c(5,8,16,17)],perm=3999)
ef
plot(ef)

rotsA=mm$CA$u[,1:2]
rotsA[,2]=(-1)*rotsA[,2]
#----------------
# map:

# color hue = quadrant of PCA plot
# alpha level = intensity of response (sd of raw among-decadal variation)
# dot size - reef area
angle0 <- function(M,N=c(1,0)){
  aa=acos( sum(M*N) / ( sqrt(sum(M*M)) * sqrt(sum(N*N)) ) )
  if (M[2]<0) { aa=2*pi-aa }
  len=sqrt(sum(M*M))
  return (c(aa,len))
}
colnames(rotsA)=c("PC1","PC2")
rr=rotsA
dat=data.frame(rr)
head(dat)
summary(dat)
dat[,1]=dat[,1]+0.03 # no change is at PC1=-0.05
d2=c()
for(i in 1:nrow(dat)) { d2=data.frame(rbind(d2,angle0(dat[i,]))) }
summary(d2)
d2[row.names(rr)=="x5",]
d2[,1]=d2[,1]*4/(2*pi)

ggplot()+geom_point(data=data.frame(rr), aes(PC1,PC2,colour=d2[,1]),alpha=0.5)+scale_colour_gradientn(colours=c("grey40","cyan3","green3","coral","grey40"))+geom_point(data=data.frame(rr[681:688,]), aes(PC1,PC2),pch=1,size=5,col="grey40")+theme_void()
ggplot()+geom_point(data=data.frame(rr[1:680,]), aes(PC1,PC2,colour=d2[1:680,1]),alpha=0.5)+scale_colour_gradientn(colours=c("grey40","cyan3","green3","coral","grey40"))+theme_void()

ssr=sst
ssr$resp.type=d2[1:680,1]
load("RE_sin45.RData")
ssr$magnitude=apply(magnitudes[1:680,],1,mean)

rasin45=rotsA
plot(rasin45[1:680,1]~ra45$CA$u[1:680,1],pch=16,col=rgb(0,0,0,alpha=0.25),cex=0.75)
plot((-1)*rasin45[1:680,2]~ra45$CA$u[1:680,2],pch=16,col=rgb(0,0,0,alpha=0.25),cex=0.75)
plot(ra85$CA$u[1:680,1]~rasin45[1:680,1],pch=16,col=rgb(0,0,0,alpha=0.25),cex=0.75)
plot((-1)*ra85$CA$u[1:680,2]~rasin45[1:680,2],pch=16,col=rgb(0,0,0,alpha=0.25),cex=0.75)
summary(lm((-1)*ra85$CA$u[1:680,2]~rasin45[1:680,2]))

ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=resp.type,cex=magnitude))+scale_colour_gradientn(colours=c("grey40","cyan3","green3","coral","grey40")) + theme_void()+ggtitle("average")

ggplot(ssr,aes(x_cent,y_cent,label=reef_id))+geom_point(aes(colour=resp.type,alpha=magnitude,cex=km2))+scale_colour_gradientn(colours=c("grey40","cyan3","green3","coral","grey40")) + theme_void()+ggtitle("average")
ggsave("map_average_sin45_age.pdf")

#----------
# plotting PCA and envfit results with gglplot

ll=load("sst_mig.RData")
names(sst)
ef=envfit(rotsA[1:680,],sst[,c(5,8,16,17)],perm=3999)
                   # PC1      PC2     r2  Pr(>r)    
# MEAN_meanT     0.93752  0.34792 0.1525 0.00025 ***
# MEAN_DT_RCP45  0.38189 -0.92421 0.0307 0.00025 ***
# income        -0.75161  0.65961 0.0496 0.00025 ***
# logarea        0.74313 -0.66915 0.0206 0.00150 ** 

colnames(rotsA)=c("PC1","PC2")
dat=rotsA[1:680,]
label_offset=1.5 # offset for environment labels
scale.env.arrows=2
scale.datapoints=12
color_arrows="red"
color_sites="red"
color_points=rgb(0,0,0,alpha=0.3)

str(ef)
ef.df<-as.data.frame(ef$vectors$arrows*scale.env.arrows*sqrt(ef$vectors$r))
ef.df2<-as.data.frame(ef$vectors$arrows*scale.env.arrows*sqrt(ef$vectors$r)*label_offset)
ef.df2$species<-rownames(ef.df2)
dats=data.frame(apply(dat,2,function(x) { x*scale.datapoints}))
str(dats)
ggplot()+
	geom_point(data=dats, aes(PC1,PC2,colour=d2[1:680,1]),alpha=1)+
	scale_colour_gradientn(colours=c("grey40","cyan3","green3","coral","grey40"))+
#	geom_text(data=data.frame(scores(mm)[["species"]]),aes(PC1,PC2),colour=color_sites,label=as.character(c(1:80)),cex=2.5)+
	geom_segment(data=ef.df[ef$vectors$pvals<0.05,],aes(x=0,xend=PC1,y=0,yend=PC2),arrow = arrow(length = unit(0.2, "cm")),colour=color_arrows)+
	geom_text(data=ef.df2[ef$vectors$pvals<0.05,],aes(PC1,PC2,label=species),col=color_arrows,size=3)+
	xlim(-2,1.8)+
	theme_void()
#	geom_point(data=data.frame(dats[681:688,]), aes(PC1,PC2),pch=1,size=5,col="grey70")

adonis(stdm[1:680,]~income+logarea+MEAN_meanT+MEAN_DT_RCP45,data=sst,method="manhattan")
               # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# income          1     11820   11820  17.732 0.02225  0.001 ***
# logarea         1      3123    3123   4.686 0.00588  0.016 *  
# MEAN_meanT      1     61107   61107  91.672 0.11505  0.001 ***
# MEAN_DT_RCP45   1      5152    5152   7.730 0.00970  0.001 ***
# Residuals     675    449943     667         0.84712           
# Total         679    531145                 1.00000           


#------- quadrants:

load("RE_45.RData")
#load("rdas.RData")
ll=load("sst_mig.RData")
#rr=ra.sin.45$CA$u[1:680,1:2]
rr=rotsA[1:680,]
dat=data.frame(rr)
dat[,1]=dat[,1]+0.03 # no change is at PC1=-0.05
plot(dat)
abline(h=0,lty=3)
abline(v=0,lty=3)

angle0 <- function(M,N){
  aa=acos( sum(M*N) / ( sqrt(sum(M*M)) * sqrt(sum(N*N)) ) )
  if (M[2]<0) { aa=2*pi-aa }
  len=sqrt(sum(M*M))
  return (c(aa,len))
}

d2=c()
for(i in 1:nrow(dat)) { d2=data.frame(rbind(d2,angle0(dat[i,],c(1,0)))) }
d2[,1]=d2[,1]*4/(2*pi)
ssr=sst
ssr$resp.type=d2[,1]
ssr$magnitude=apply(magnitudes[1:680,],1,mean)

rws=raws[[run[1]]]
str(raws[[r]])
str(rws)
for (r in run[-1]) {
	rws=rws+raws[[r]]
}

rws=rws/length(run)

str=stack(data.frame(t(rws[1:680,])))
str(rws)
levels(str$ind)
str$reef_id=sub("X","",str$ind)
str$dec=seq(1:40)
names(str)[1]="cover"
covsst=merge(str,ssr,by="reef_id",all.x=T)
str(covsst)
#covsst$dec=factor(covsst$dec)

c1=subset(covsst,resp.type>3 | resp.type<0 )
c1$resp.type[c1$resp.type>3]=c1$resp.type[c1$resp.type>3]-4
c1$reef_id=factor(c1$reef_id, levels=ssr$reef_id[order(ssr$resp.type,decreasing=T)])
head(c1)
ggplot(c1,aes(dec,cover,group=reef_id))+geom_line(aes(colour=resp.type,alpha=magnitude))+scale_colour_gradientn(colours=c("coral","grey25"))+theme_bw()
ggsave("quadrant4_re45.pdf",device="pdf")

c1=subset(ssr,resp.type>0 & resp.type<=1 )
hist(c1$resp.type)
c11=subset(covsst,resp.type>0 & resp.type<=1 )
c11$reef_id=factor(c11$reef_id, levels=ssr$reef_id[order(ssr$magnitude,decreasing=T)])
ggplot(c11,aes(dec,cover,group=reef_id))+geom_line(aes(colour=resp.type,alpha=magnitude))+scale_colour_gradientn(colours=c("grey25","cyan3"))+theme_bw()
ggsave("quadrant1_re45.pdf",device="pdf")

c1=subset(ssr,resp.type>1 & resp.type<=2 )
hist(c1$resp.type)
c11=subset(covsst,resp.type>1 & resp.type<=2 )
c11$reef_id=factor(c11$reef_id, levels=ssr$reef_id[order(ssr$magnitude)])
ggplot(c11,aes(dec,cover,group=reef_id))+geom_line(aes(colour=resp.type,alpha=magnitude))+scale_colour_gradientn(colours=c("cyan3","green3"))+theme_bw()
ggsave("quadrant2_re45.pdf",device="pdf")

c1=subset(ssr,resp.type>2 & resp.type<=3 )
hist(c1$resp.type)
c11=subset(covsst,resp.type>2 & resp.type<=3 )
c11$reef_id=factor(c11$reef_id, levels=ssr$reef_id[order(ssr$magnitude)])
ggplot(c11,aes(dec,cover,group=reef_id))+geom_line(aes(colour=resp.type,alpha=magnitude))+scale_colour_gradientn(colours=c("green3","coral"))+theme_bw()
ggsave("quadrant3_re45.pdf",device="pdf")

#----------

# pc1: high - overall decline under warming, low - increase under warming
ggplot(sst,aes(x_cent,y_cent,colour=pc1.85.5.esd,label=reef_id))+geom_point(cex=(sst$logarea^4)/50)+scale_colour_gradient( low="blue",high="coral") + theme_void()
# pc2: high - initial decline then recovery, low - initial persistence (or increase) then decline
ggplot(sst,aes(x_cent,y_cent,colour=pc2.85.5.esd,label=reef_id))+geom_point(cex=(sst$logarea^4)/50)+scale_colour_gradient( low="blue",high="coral") + theme_void()


#--------
# PCA circle heatmaps

setwd("~/Dropbox/Documents/coral_triangle_SLiM/slimQTL_v6/")
ll=load("standardized_cover.RData")
ll=load("sst_mig.RData")


# 85.5.esd:
plot(pc2.85.5.esd~pc1.85.5.esd,sst,pch=16,col=rgb(0,0,0,alpha=0.2),xlim=c(-0.16,0.08),ylim=c(-0.15,0.15),xlab="PC1",ylab="PC2",mgp=c(2.3,1,0))
abline(h=0,lty=3,col="grey70")
abline(v=0,lty=3,col="grey70")
reefs=identify(sst$pc1.85.5.esd,sst$pc2.85.5.esd,labels=sst$reef_id,n=40,col="red",cex=0.7)
reefs=c(32,33,378,38,384,36,379,410,388,497,125,565,671,567,672,159,152,195,246,523,189,184,349,501,194,351,302,210,3,628,610,594,168,577,595,181,21,231,390,393,383,130)
#reefs=c(54,143, 276, 422 ,619)
pheatmap(std.mid[rev(reefs),],cluster_cols=F,cluster_rows=F,clustering_distance_rows="correlation",clustering_method="average",border_color=F)



ss=t(cov85.5.esd)
str(ss)
length(seq(1,800,10))
warm.year=seq(-390,400,10)
length(warm.year)
nrow(ss)
ss=cbind(ss,warm.year)
head(ss)
plot(ss[,1]~warm.year,bty="n",xlab="year of warming",yaxt="n",ylab="",col="white",mgp=c(2.3,1,0))

plot(mm$CA$u[,3]~sst$pc2.85.5.esd,pch=16,col=rgb(0,0,0,alpha=0.2),xlim=c(-0.16,0.15),ylim=c(-0.12,0.17))
abline(h=0,lty=3)
abline(v=0,lty=3)
identify(sst$pc2.85.5.esd,mm$CA$u[,3],labels=sst$reef_id,n=40,col="red",cex=0.7)
reefs=c(607,250,8,349,501,428,426,295,475,282,383,231,131,677,605,404,403,667,646,28)
pheatmap(std.mid[reefs,],cluster_cols=F,cluster_rows=F,clustering_distance_rows="correlation",clustering_method="average",border_color=F)

# 85esd:
plot(pc2.85esd~pc1.85esd,sst,pch=16,col=rgb(0,0,0,alpha=0.2),xlim=c(-0.16,0.08),ylim=c(-0.15,0.15))
abline(h=0,lty=3)
abline(v=0,lty=3)
identify(sst$pc1.85esd,sst$pc2.85esd,labels=sst$reef_id,n=40,col="red",cex=0.7)
reefs=c(593,635,329,190,151,152,172,326,511,599,194,598,510,302,210,3,675,594,21,168,507,390,167,13,181,605,502,231,206,400,408,411,410,386,399,491,476,496,284)
pheatmap(cov85esd[rev(reefs),],cluster_cols=F,cluster_rows=F,clustering_distance_rows="correlation",clustering_method="average",border_color=F)

# 85:
plot(pc2.85~pc1.85,sst,pch=16,col=rgb(0,0,0,alpha=0.2),xlim=c(-0.16,0.08),ylim=c(-0.16,0.13))
abline(h=0,lty=3)
abline(v=0,lty=3)
identify(sst$pc1.85,sst$pc2.85,labels=sst$reef_id,n=40,col="red",cex=0.7)

reefs=c(593,527,95,470,188,419,173,388,331,194,218,214,247,246,511,351,598,512,664,3,13,507,615,605,663,612,119,120,595,572,224,677,21,411,407,385,5,491,409,476,496,37,92,404,74)
pheatmap(std.mid[reefs,],cluster_cols=F,cluster_rows=F,clustering_distance_rows="correlation",clustering_method="average",border_color=F)


cov85=std.mid
save(cov85.5.esd,cov85esd,cov85,file="standardized_cover.RData")

ll=lm(pc1~income+logarea+MEAN_meanT+MEAN_DT_RCP45,sst)
summary(ll)
