bruno=read.csv("Global_Cover_1_5_10_bruno_may28_2019.csv")
load("sst_modified.RData")
names(sst)
ll=load('reef_responses.RData')
sst$rr50=rls50[,which(meta[,1]=="8" & meta[,2]=="SIN")]
sst$rr100=rls100[,which(meta[,1]=="8" & meta[,2]=="SIN")]

bruno=subset(bruno, LONGITUDE>min(sst$lon)-1 & LONGITUDE<max(sst$lon)+1 & LATITUDE>min(sst$lat)-1 & LATITUDE<max(sst$lat)+1)

install.packages("gpclib")
library(gpclib)
if (!rgeosStatus()) gpclibPermit()
setwd("~/Dropbox/Documents/coral_triangle_SLiM/")
mapfile="/HiRes_shoreline/gshhs_l.b"
sf1 <- getRgshhsMap(mapfile, xlim = c(min(sst$lon)-1, max(sst$lon)+1), ylim = c(min(sst$lat)-5, max(sst$lat)+3)) 
ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey90', color='black', lwd = 0.1) +
  geom_point(data=sst, aes(x=lon, y=lat,color=lprophots,cex=km2),alpha=0.8)+
  geom_point(data=bruno, aes(x=LONGITUDE, y=LATITUDE),color="cyan3",cex=0.5)+
  scale_color_viridis(option="B")+
  theme_bw()


lat=c(sst$lat,bruno$LATITUDE)
lon=c(sst$lon,bruno$LONGITUDE)
treefs=paste("t",sst$reef_id,sep="")
breefs=paste("b",bruno$REEF_ID_GR,sep="")
coords=data.frame(cbind(lat,lon,"reef"=c(treefs,breefs)))
c2=c()
for (r in unique(coords$reef)) {
	c2=data.frame(rbind(c2,subset(coords,reef==r)[1,]))
}
row.names(c2)=c2$reef
br=grep("b",c2$reef)
tr=grep("t",c2$reef)
c2$set="bruno"
c2$set[tr]="tri"

hc=hclust(dist(c2[,c("lat","lon")]),method="complete")
c2$grp=cutree(hc,h=1)


curves=list();c=1;goods=c();breefs=c();treefs=c();sizes=c()
for (g in unique(c2$grp)) {
		message(g)
		ss=subset(c2,grp==g)
		if("bruno" %in% ss$set &  "tri" %in% ss$set) {
			sizes[g]=nrow(ss)
			goods[g]=1
			sb=subset(ss,set=="bruno")
			st=subset(ss,set=="tri")
			breefs[g]=nrow(sb)
			treefs[g]=nrow(st)
			for (i in 1:nrow(sb)){
				for (j in 1:nrow(st)){
					curves[[cc]]=data.frame(
					x1=sb[i,"lon"],
					y1=sb[i,"lat"],
					x2=st[j,"lon"],
					y2=st[j,"lat"]
					)
					cc=cc+1
				}
			}
		}
}
goods
plot(sizes)

length(curves)
curv=do.call(rbind,curves)
curv[,1]=as.numeric(as.character(curv[,1]))
curv[,2]=as.numeric(as.character(curv[,2]))
curv[,3]=as.numeric(as.character(curv[,3]))
curv[,4]=as.numeric(as.character(curv[,4]))
head(curv)

#-------- which bruno samples are near triangle?
ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey90', color='black', lwd = 0.1) +
  geom_point(data=sst, aes(x=lon, y=lat),color="blue",cex=0.3)+
  geom_point(data=bruno, aes(x=LONGITUDE, y=LATITUDE),color="red",cex=2,alpha=0.2)+
  geom_segment(data=curv,aes(x=x1,y=y1,xend=x2,yend=y2),size=1,color="cyan3")+
  theme_bw()


#------------

c2$lat=as.numeric(as.character(c2$lat))
c2$lon=as.numeric(as.character(c2$lon))
cc=1
g=475; dec=list(); cc=1
for (g in unique(c2$grp)) {
	message(g)
	ss=subset(c2,grp==g)
	if("bruno" %in% ss$set &  "tri" %in% ss$set) {
			lon=mean(ss$lon)
			lat=mean(ss$lat)
			sb=sub("b","",as.character(subset(ss,set=="bruno")$reef))
			ye=c();pe=c()
			for (b in sb){
				brd=subset(bruno,REEF_ID_GR==b)
				if (range(brd$YEAR)[2]-range(brd$YEAR)[1]>10) {
					ye=append(ye,summary(lm(HARD_COR_P~YEAR,brd))$coefficients[2,1])
					pe=append(pe,summary(lm(HARD_COR_P~YEAR,brd))$coefficients[2,4])
				}
			}
			if(length(ye)>0) {
#					plot(HARD_COR_P~YEAR,brd)
				st=sub("t","",as.character(subset(ss,set=="tri")$reef))
				d50=mean(subset(sst,reef_id %in% st)$rr50)
				d100=mean(subset(sst,reef_id %in% st)$rr100)
				lph=mean(subset(sst,reef_id %in% st)$lprophots)
				dec[[cc]]=data.frame(cbind(lon=lon,lat=lat,ye=mean(ye),pe=mean(pe),d50=d50,d100=d100,lph=lph))
				cc=cc+1
			}
	}
}
dec=do.call(rbind,dec)
str(dec)


ggplot(dec,aes(lph,ye))+geom_point()
ggplot(dec,aes(d50,ye))+geom_point()
ggplot(dec,aes(d100,ye))+geom_point()
ggplot(dec,aes(lph,d50))+geom_point()
summary(lm(change~lprophots,sst2))


ggplot() + 
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = 'grey90', color='black', lwd = 0.1) +
#  geom_point(data=dec, aes(x=lon, y=lat),color="blue",cex=0.3)+
  geom_point(data=dec, aes(x=lon, y=lat,color=ye,size=-log(pe,10)))+
 # geom_segment(data=curv,aes(x=x1,y=y1,xend=x2,yend=y2),size=1,color="cyan3")+
  scale_color_viridis(option="B")+
  theme_bw()

