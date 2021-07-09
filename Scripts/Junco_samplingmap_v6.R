require("sp")
require("rgdal")
require("xlsx")
require("maps")
require("mapdata")
require("raster")
require("rgeos")
require("gplots")
require("rvertnet")
require("dismo")
require("maptools")
require("RColorBrewer")
require("dplyr")
require("colorspace")
require("png")

### Read in pngs ###
junco_male<-readPNG("~/Desktop/Manuscripts/Juncos/Plates/JuncoMale.png")
junco_female<-readPNG("~/Desktop/Manuscripts/Juncos/Plates/JuncoFemale.png")
male_symbol<-readPNG("~/Desktop/Manuscripts/Juncos/Plates/male.png")
female_symbol<-readPNG("~/Desktop/Manuscripts/Juncos/Plates/female.png")

### Read in our data set from csv ###
orju_specimens<-read.csv("~/Desktop/Manuscripts/Juncos/Data/OregonJuncoColoration_AppendixA_v6.csv",stringsAsFactors=F)
orju_specimens$Subspecies<-factor(orju_specimens$Subspecies,levels=rev(c("oreganus","montanus","shufeldti","thurberi","pinosus")))

### Read in GADM Shapefiles ###
usa0<-readOGR("~/GIS/GADM/USA_adm/USA_adm0.shp")
can0<-readOGR("~/GIS/GADM/CAN_adm/CAN_adm0.shp")
mex0<-readOGR("~/GIS/GADM/MEX_adm/MEX_adm0.shp")
usa0.simp<-gSimplify(usa0,0.05)
can0.simp<-gSimplify(can0,0.05)
mex0.simp<-gSimplify(mex0,0.05)
na0.simp<-rbind(usa0.simp,can0.simp,mex0.simp, makeUniqueIDs = TRUE)

usa1<-readOGR("~/GIS/GADM/USA_adm/USA_adm1.shp")
can1<-readOGR("~/GIS/GADM/CAN_adm/CAN_adm1.shp")
mex1<-readOGR("~/GIS/GADM/MEX_adm/MEX_adm1.shp")
usa1.simp<-gSimplify(usa1,0.05)
can1.simp<-gSimplify(can1,0.05)
mex1.simp<-gSimplify(mex1,0.05)
na1.simp<-rbind(usa1.simp,can1.simp,mex1.simp, makeUniqueIDs = TRUE)

usa2<-readOGR("~/GIS/GADM/USA_adm/USA_adm2.shp")
can2<-readOGR("~/GIS/GADM/CAN_adm/CAN_adm2.shp")
mex2<-readOGR("~/GIS/GADM/MEX_adm/MEX_adm2.shp")
usa2.simp<-gSimplify(usa2,0.05)
can2.simp<-gSimplify(can2,0.05)
mex2.simp<-gSimplify(mex2,0.05)
na2.simp<-rbind(usa2.simp,can2.simp,mex2.simp, makeUniqueIDs = TRUE)

### Read in elevational data ###
na_alt<-raster("/Users/NickMason/Google Drive/StellersJay/GIS/na_alt.grd")

### Read in Junco Range Map ###
junco<-readOGR("~/Desktop/Manuscripts/Juncos/GIS/junco_breeding.shp")

### Read in modified subspecies ranges ###
ssp_shps<-list.files("~/Desktop/Manuscripts/Juncos/GIS/breedingrangeshapefiles",pattern=".shp",full.names=T)
ssp_shps_short<-list.files("~/Desktop/Manuscripts/Juncos/GIS/breedingrangeshapefiles",pattern=".shp",full.names=F)
ssp_list<-list()
for(i in 1:length(ssp_shps)){
	ssp_list[[i]]<-readOGR(ssp_shps[i])
}

names(ssp_list)<-gsub("[.]shp","",ssp_shps_short)
ssp_list<-ssp_list[tolower(levels(orju_specimens$Subspecies))] #Reorder list based on factor in dataframe

### Determine winter / breeding status ### Breeding defined as inclusive April through August
orju_specimens$breeding<-c("winter","breeding")[as.numeric(orju_specimens$Month >=4 & orju_specimens$Month <=8)+1]

### Sampling map with separate subspecies colors and shapes for breeding vs winter ###
junco_cols<-qualitative_hcl(length(table(orju_specimens$Subspecies)), "Dark 3")

### Individual pngs ### 

for(i in 1:length(levels(orju_specimens$Subspecies))){
png(file=paste0("~/Desktop/Manuscripts/Juncos/Figures/JuncoSamplingMap_",names(ssp_list)[i],"_v2.png"),width=3.25,height=4.25,units="in",res=400)
#quartz(width=3.25,height=4.25)
par(mar=c(0,0,0,0))

plot(na2.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray90")
plot(na1.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray70",add=T)
plot(na0.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray40",add=T)

plot(ssp_list[[i]],col= paste0(junco_cols[i],"80"),border=NA,add=T)

junco_spec_spp<-orju_specimens[orju_specimens$Subspecies==levels(orju_specimens$Subspecies)[i],]
points(junco_spec_spp$Longitude, junco_spec_spp$Latitude,pch=21,bg=paste0(junco_cols[i],50))
text(label=paste0(junco_spec_spp$Museum,junco_spec_spp$Catalog.Number),x=as.numeric(junco_spec_spp$Longitude), y=as.numeric(junco_spec_spp$Latitude)+0.1,pch=21,cex=0.4)

	# }else{
		# points(junco_spec_spp$Longitude[junco_spec_spp$breeding=="breeding"], junco_spec_spp$Latitude[junco_spec_spp$breeding=="breeding"],pch=21,bg=paste0(junco_cols[i],70))
		# points(junco_spec_spp$Longitude[junco_spec_spp$breeding=="winter"], junco_spec_spp$Latitude[junco_spec_spp$breeding=="winter"],pch=22,bg=paste0(junco_cols[i],70))
	# }	
dev.off()
}

png(file="~/Desktop/Manuscripts/Juncos/Figures/JuncoSamplingMap_V3.png",width=3.25,height=4.25,units="in",res=400)
#quartz(width=3.25,height=4.25)
par(mar=c(0,0,0,0))

plot(na2.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray90")
plot(na1.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray70",add=T)
plot(na0.simp,xlim=c(-137.8534,-107.7421),ylim=c(31.17948,59.16771),border="gray40",add=T)

for(i in 1:length(levels(orju_specimens$Subspecies))){
	plot(ssp_list[[i]],col= paste0(junco_cols[i],"80"),border=NA,add=T)
}

for(i in 1:length(levels(orju_specimens$Subspecies))){
	junco_spec_spp<-orju_specimens[orju_specimens$Subspecies==levels(orju_specimens$Subspecies)[i],]
	# if(i==1){
		points(junco_spec_spp$Longitude, junco_spec_spp$Latitude,pch=21,bg=paste0(junco_cols[i],50))
	# }else{
		# points(junco_spec_spp$Longitude[junco_spec_spp$breeding=="breeding"], junco_spec_spp$Latitude[junco_spec_spp$breeding=="breeding"],pch=21,bg=paste0(junco_cols[i],70))
		# points(junco_spec_spp$Longitude[junco_spec_spp$breeding=="winter"], junco_spec_spp$Latitude[junco_spec_spp$breeding=="winter"],pch=22,bg=paste0(junco_cols[i],70))
	# }	
}

### plot PNGs ###
rasterImage(junco_male,xleft=-138.3311,ybottom= 47.46083,xright=-129.7868,ytop= 51.54730)
rasterImage(male_symbol,xleft=-133.5,ybottom= 50.45453,xright=-132.4365,ytop= 51.2)

rasterImage(junco_female,xleft=-138.3311,ybottom= 42,xright=-129.7868,ytop= 46.08)
rasterImage(female_symbol,xleft=-133.2976,ybottom= 45,xright=-132.4365,ytop= 46)

### key of colors ###
points(x=rep(-137.6,5),y=seq(39.5,33.5,length.out=5),pch=22,col=NA,bg=rev(junco_cols))
text(labels= rev(levels(orju_specimens$Subspecies)),x=rep(-137,5),y=seq(39.6,33.6,length.out=5),adj=c(0,0.5),cex=0.7,font=3)

#scalebar(100,xy=c(-125.6,31.1),label="",type="bar",lwd=0.5)
map.scale(x=-126,y=32,cex=0.5,ratio=F)
dev.off()
