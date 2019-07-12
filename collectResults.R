# ------ cartoons:
library(ggplot2)
x=rep(0.2,40)
b=seq(0,0.3,length.out=20)
y=c(1:40)

x1=x2=x3=x4=x5=x6=x7=x8=x
x1[y>20]=0.05
x2[y>20]=x[y>20]-0.15+b/2
x3[y>20]=x[y>20]-0.15+b
x4[y>20]=x[y>20]+b/2
x5[y>20]=0.35
x6[y>20]=x[y>20]+0.15-b/2
x7[y>20]=x[y>20]+0.15-b
x8[y>20]=x[y>20]-b/2

toons=data.frame(cbind(x1,x2,x3,x4,x5,x6,x7,x8,x1))
names(toons)=paste("type",c(1:9),sep="")
stoons=stack(toons)
names(stoons)=c("cover","resp.type")
stoons$dec=c(1:40)
stoons$rt=(as.numeric(stoons$resp.type)-1)/2

ggplot(stoons,aes(dec,cover,group=resp.type))+geom_line(aes(colour=rt),size=2)+scale_colour_gradientn(colours=c("plum","cyan3","green3","coral","plum"))  + theme(panel.background = element_blank(),panel.grid.major = element_blank(),strip.background = element_blank(),strip.text.x = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1))+geom_vline(xintercept=20,lty=3)+facet_wrap(~resp.type)+ylim(0.01,0.4)

toons=data.frame(cbind(x1,x2,x3,x4,x5,x6,x7,x8))

#------------

library(tidyverse)
library(Rmisc)
library(vegan)
library(pheatmap)
ll=load("sst_mig.RData")
getwd()
run=c("A","B","C","D","E")
for (s in c("85","45") ){
	respranges=list();prewarms=list();raws=list();stds=list();relresps=list();hs=list()
	r="C";sets=s;index=""
	message(sets)
	for (r in run){
		message(paste("   ",r))
		g85=read.table(paste("tri",sets,r,index,".clean",sep=""),sep="\t",skip=19)
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

