require("rvertnet")
require("colorspace")
require("pavo")
require("scatterplot3d")
require("rgl")
require("agricolae")#Provides HSD.test()
require("plotrix")
require("MuMIn")
require("mclust")
require("RColorBrewer")

### DFA ###
require(MASS)
citation("MASS")

setwd("~/oregonjuncocolor")

### Read in Junco csv file ### 
junco_data<-read.csv("Data/OregonJuncoColoration_AppendixA_v6.csv",stringsAsFactors=F)

junco_data$Subspecies<-factor(junco_data$Subspecies,levels=c("oreganus","montanus","shufeldti","thurberi","pinosus"))

junco_cols<-qualitative_hcl(length(table(junco_data$Subspecies)), "Dark 3")

### 206 total specimens ###
nrow(junco_data)

### What if we restrict to just breeding specimens? ###
junco_data$Month<-as.numeric(junco_data$Month)
junco_data_breeding<-junco_data[4 <= junco_data$Month & junco_data$Month <=8,]
nrow(junco_data_breeding)

### As of v 9 [02 april 2020] assigning unknown age classes as adults ###
junco_data$Age[junco_data$Age=="U"]<-"A"

### include factor that correponds to breeding vs non-breeding ###
junco_data$breeding<-factor(c("non-breeding","breeding")[as.numeric(4 <= junco_data$Month & junco_data$Month <=8)+1])
head(junco_data)

### Subspecies counts ###
table(junco_data$Subspecies)
table(junco_data_breeding$Subspecies)

### Figure out the days since molt ###
fallmolt<-strptime("09.01.2000",format="%m.%d.%Y")
springmolt<-strptime("03.15.2000",format="%m.%d.%Y")
collectiondays<-strptime(paste(as.character(junco_data$Month),as.character(junco_data$Day),"2000",sep="."),format="%m.%d.%Y")

dayssincemolt<-rep(NA,nrow(junco_data))
for(i in 1:nrow(junco_data)){
	foo<-round(difftime(collectiondays[i],fallmolt,units="days"),0)
	if(foo<0){
		foo<-round(difftime(collectiondays[i],springmolt,units="days"),0)
	}
	if(foo<0){
		foo<-round(foo+difftime(strptime("12.31.2000",format="%m.%d.%Y"),strptime("09.01.2000",format="%m.%d.%Y"))+difftime(strptime("03.15.2000",format="%m.%d.%Y"),strptime("01.01.2000",format="%m.%d.%Y")),0)
	}
	dayssincemolt[i]<-foo
}
junco_data$dayssincemolt<-dayssincemolt

### Calculate the age of each specimen ###
junco_data$specimenage<-2018-as.numeric(junco_data$Year)

### Create a factor for junco subspecies, ordered from north to south ###
junco_data$Subspecies<-factor(junco_data$Subspecies,levels=c("oreganus","montanus","shufeldti","thurberi","pinosus"))
junco_data$Age <-factor(junco_data$Age,levels=c("A","I"))
junco_data$Sex <-factor(junco_data$Sex,levels=c("M","F"))

### Simplify GLM methods to just look at additive effects of subspecies, sex, specimen age, and days since molt 
glm_output<-list()
for(i in 1:6){
	junco_data[,(15:20)[i]]<-as.numeric(junco_data[,(15:20)[i]])
	glm_output[[i]]<-glm(formula(paste0(colnames(junco_data)[15:20][i],"~ Sex+Age+Subspecies+dayssincemolt+specimenage")),data=junco_data)	
}
names(glm_output)<-colnames(junco_data)[15:20]

### Create table of GLM output for SI / main ###
png(file="Figures/GLM_residuals_v4.png",width=6.5,height=6.5*(3/2),units="in",res=500)
par(mfrow=c(3,2))
foo_tab<-list()

for(i in 1:length(glm_output)){
	foo_tab[[i]]<-round(summary(glm_output[[i]])$coefficients,3)
	foo_tab[[i]][,4]<-round(summary(glm_output[[i]])$coefficients,4)[,4]
	write.csv(foo_tab[[i]],file=paste0("Tables/GLMoutput_",colnames(junco_data)[15:20][i],"_v5.csv"))
	hist(residuals(glm_output[[i]]),main=colnames(junco_data)[15:20][i],xlab="Residuals")
}
dev.off()

### Correct values by specimen age for i = 3, 5 , 6###
for(i in 1:length(glm_output)){
	if(summary(glm_output[[i]])$coefficients["specimenage",4]<0.05){
		junco_data[,paste0(names(glm_output)[i],"_","agecorrected")]<-junco_data[,(15:20)[i]]+ junco_data$specimenage * summary(glm_output[[i]])$coefficients["specimenage",1]
	}
}

### Get effect sizes of specimen ages out of models to report ###
names(glm_output)[3]
summary(glm_output[[3]])$coefficients['specimenage',1:2]
names(glm_output)[5]
summary(glm_output[[5]])$coefficients['specimenage',1:2]
names(glm_output)[6]
summary(glm_output[[6]])$coefficients['specimenage',1:2]

### BOX PLOTS BY SEX ###
### Set up colors for plotting ###
junco_sexes_col<-c("gray95","gray60")

### Set up sex as a factor ###
junco_data$Sex<-factor(junco_data$Sex)

#quartz(width=3.25,height=(3.25*2/3))

png(file="~/Desktop/Manuscripts/Juncos/Figures/LABboxplots_sex_v6.png",width=3.25,height=(3.25*2/3),units="in",res=400)
par(mfrow=c(2,3))
par(mar=c(0.5,2.25,0.5,0.5))
par(oma=c(2.5,0,0,0))
par(lwd=0.5)

### Make sure to use corrected values for Hood.b, Back. a, and Back.b!
ncol(junco_data)

### Box plots and Tukeys for sex and ssp ###
for(i in 1:6){
	char_glm<-glm(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Subspecies+Sex+Age+dayssincemolt+specimenage")),data=junco_data)	
	char_aov<-aov(char_glm)
	
	char_HSD<-HSD.test(char_aov,"Sex",group=T,console=T)
	char_HSD$groups<-char_HSD$groups[levels(junco_data$Sex),]
	
	### Set up y limits ###
	char_ylim<-c(range(junco_data[,(c(15,16,24,18,25,26))[i]])[1],diff(range(junco_data[,(c(15,16,24,18,25,26))[i]]))*0.2 + range(junco_data[,(c(15,16,24,18,25,26))[i]])[2])
	
	### Generate plot ###
	char_bp<-boxplot(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Sex")),data=junco_data,col= junco_sexes_col,ylim=char_ylim,axes=F,xlab="",ylab="",lwd=0.5)
	box()
	axis(2,mgp=c(0,0.25,0),cex.axis=0.7,tck=-0.05,lwd=0.5)
	
	### Add Y labels ###
	mtext(text=gsub("[.]"," ",colnames(junco_data[c(15,16,24,18,25,26)])[i]),side=2,line=1.2,cex=0.6)
	
	### Add x labels ###
	par(xpd=NA)
	if(i %in% 4:6){
		text(c("female","male"),x=1:2,y=par("usr")[3]-diff(par("usr")[3:4])*0.075,srt=90,adj=c(1,0.5),cex=0.8,font=1)
		}
	par(xpd=T)

	### Add Tukey's Labels ###
	text(label=toupper(char_HSD$groups$groups),y= diff(range(junco_data[,(c(15,16,24,18,25,26))[i]]))*0.15 + range(junco_data[,(c(15,16,24,18,25,26))[i]])[2],x=1:2,cex=0.8)
	
	### Add figure labels ###
	text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.08,y=par("usr")[3]+diff(c(par("usr")[3],par("usr")[4]))*0.08,label=LETTERS[i],font=2,cex=0.8)
	
}
dev.off()

### BOX PLOTS BY Age ###
### Set up colors for plotting ###
junco_ages_col<-c("gray95","gray60")

### Set up sex as a factor ###
junco_data$Age<-factor(junco_data$Age,levels=c("I","A"))

#quartz(width=3.25,height=(3.25*2/3))

png(file="~/Desktop/Manuscripts/Juncos/Figures/LABboxplots_age_v5.png",width=3.25,height=(3.25*2/3),units="in",res=400)
par(mfrow=c(2,3))
par(mar=c(0.75,2.25,0.5,0.5))
par(oma=c(2.5,0,0,0))
par(lwd=0.5)

### Box plots and Tukeys for sex and ssp ###
for(i in 1:6){
	junco_data[,(c(15,16,24,18,25,26))[i]]<-as.numeric(junco_data[,(c(15,16,24,18,25,26))[i]])
	char_glm<-glm(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Subspecies+Sex+Age+dayssincemolt+specimenage")),data=junco_data)	
	char_aov<-aov(char_glm)
	
	char_HSD<-HSD.test(char_aov,"Age",group=T,console=T)
	char_HSD$groups<-char_HSD$groups[levels(junco_data$Age),]
	
	### Set up y limits ###
	char_ylim<-c(range(junco_data[,(c(15,16,24,18,25,26))[i]])[1],diff(range(junco_data[,(c(15,16,24,18,25,26))[i]]))*0.2 + range(junco_data[,(c(15,16,24,18,25,26))[i]])[2])
	
	### Generate plot ###
	char_bp<-boxplot(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Age")),data=junco_data,col= junco_sexes_col,ylim=char_ylim,axes=F,xlab="",ylab="",lwd=0.5)
	box()
	axis(2,mgp=c(0,0.25,0),cex.axis=0.7,tck=-0.05,lwd=0.5)
	
	### Add Y labels ###
	mtext(text=gsub("[.]"," ",colnames(junco_data[c(15,16,24,18,25,26)])[i]),side=2,line=1.2,cex=0.6)
	
	### Add x labels ###
	par(xpd=NA)
	if(i %in% 4:6){
		text(c("immature","adult"),x=1:2,y=par("usr")[3]-diff(par("usr")[3:4])*0.075,srt=90,adj=c(1,0.5),cex=0.8,font=1)
		}
	par(xpd=T)

	### Add Tukey's Labels ###
	text(label=toupper(char_HSD$groups$groups),y= diff(range(junco_data[,(c(15,16,24,18,25,26))[i]]))*0.15 + range(junco_data[,(c(15,16,24,18,25,26))[i]])[2],x=1:2,cex=0.8)
	
	### Add figure labels ###
	text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.08,y=par("usr")[3]+diff(c(par("usr")[3],par("usr")[4]))*0.08,label=LETTERS[i],font=2,cex=0.8)
	
}
dev.off()

### BOX PLOTS BY SUBSPECIES ###
### Set up colors for plotting ###
junco_cols<-qualitative_hcl(length(table(junco_data$Subspecies)), "Dark 3")

#quartz(width=6.5,height=(6.5*2/3))
png(file="~/Desktop/Manuscripts/Juncos/Figures/LAB_subspecies_boxplots_v6.png",width=6.5,height=(6.5*2/3),units="in",res=300)
par(mfrow=c(2,3))
par(mar=c(0.5,3,0.5,1))
par(oma=c(5,0,0,0))
par(lwd=1)
### Box plots and Tukeys for sex and ssp ###
i<-1

for(i in 1:6){
	char_glm<-glm(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Subspecies+Sex+Age+dayssincemolt+specimenage")),data=junco_data)	
	char_aov<-aov(char_glm)
	
	char_HSD<-HSD.test(char_aov,"Subspecies",group=T,console=T)
	char_HSD$groups<-char_HSD$groups[levels(junco_data$Subspecies),]
	
	### Set up y limits ###
	char_ylim<-c(range(junco_data[,(c(15,16,24,18,25,26))[i]])[1],diff(range(junco_data[,(c(15,16,24,18,25,26))[i]]))*0.1 + range(junco_data[,(c(15,16,24,18,25,26))[i]])[2])
	
	### Generate plot ###
	char_bp<-boxplot(formula(paste0(colnames(junco_data)[c(15,16,24,18,25,26)][i],"~Subspecies")),data=junco_data,col= rev(junco_cols),ylim=char_ylim,axes=F,xlab="",ylab="")
	box()
	axis(2,mgp=c(0,0.75,0))
	mtext(text=gsub("[.]"," ",colnames(junco_data[c(15,16,24,18,25,26)])[i]),side=2,line=1.75)
	
	par(xpd=NA)
	if(i %in% 4:6){
		text(tolower(levels(junco_data$Subspecies)),x=1:5,y=par("usr")[3]-diff(par("usr")[3:4])*0.05,srt=90,adj=c(1,0.5),cex=1.25,font=3)
		}
	par(xpd=T)

	### Add Tukey's Labels ###
	text(label=toupper(char_HSD$groups$groups),y=char_ylim[2],x=1:5)
	
	### Add figure labels ###
	text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[3]+diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[i],font=2,cex=1.5)
	
}
dev.off()

#############################################################
### Discriminant function analysis on each sex separately ###
#############################################################
junco_data_f<-junco_data[junco_data$Sex=="F" & junco_data$Age=="A",]
junco_data_m<-junco_data[junco_data$Sex=="M" & junco_data$Age=="A",]

nrow(junco_data_m)
nrow(junco_data_f)

### Males ###
junco_lda_nocv_m<-lda(Subspecies~.,data=junco_data_m[c(6,c(15,16,24,18,25,26))],CV=F)
junco_lda_nocv_m$scaling

junco_lda_m<-lda(Subspecies ~.,data= junco_data_m[c(6,c(15,16,24,18,25,26))],CV=T)

CV_tab_m <- table(junco_data_m $Subspecies, junco_lda_m $class)
CV_tab_m <- as.data.frame(matrix(CV_tab_m,nrow=5))
rownames(CV_tab_m)<-paste("actual",levels(junco_data_m$Subspecies),sep="_")
colnames(CV_tab_m)<-paste("predicted", levels(junco_data_m$Subspecies),sep="_")

diag(as.matrix(CV_tab_m))
round(sum(diag(as.matrix(CV_tab_m))) / sum(CV_tab_m) * 100,2)

junco_predict_m<-predict(junco_lda_nocv_m)

### Females ###
junco_lda_nocv_f<-lda(Subspecies~.,data=junco_data_f[c(6,c(15,16,24,18,25,26))],CV=F)
junco_lda_nocv_f$scaling

junco_lda_f<-lda(Subspecies ~.,data= junco_data_f[c(6,c(15,16,24,18,25,26))],CV=T)

CV_tab_f <- table(junco_data_f $Subspecies, junco_lda_f $class)
CV_tab_f <- as.data.frame(matrix(CV_tab_f,nrow=5))
rownames(CV_tab_f)<-paste("actual",levels(junco_data_f$Subspecies),sep="_")
colnames(CV_tab_f)<-paste("predicted", levels(junco_data_f$Subspecies),sep="_")

diag(as.matrix(CV_tab_f))
round(sum(diag(as.matrix(CV_tab_f))) / sum(CV_tab_f) * 100,2)

junco_predict_f<-predict(junco_lda_nocv_f)

### Create figure for DFA analysis ###
png(file="~/Desktop/Manuscripts/Juncos/Figures/DFA_output_v7.png",width=6.5,height=6.5/2,units="in",res=500)
#quartz(width=6.25,height=6.5/2)

par(mfrow=c(1,2))
par(mar=c(2.5,3,2,1))
plot(-junco_predict_m$x[,1],-junco_predict_m$x[,2]*-1,axes=F,pch=21,bg=paste0(rev(junco_cols)[as.numeric(junco_data_m$Subspecies)],"80"),col= rev(junco_cols)[as.numeric(junco_lda_m$class)],lwd=2)
title(main=paste0("Male DFA (",round(sum(diag(as.matrix(CV_tab_m))) / sum(CV_tab_m) * 100,2),"%)"),cex.main=1)
box()
axis(1,mgp=c(0,0.5,0),tck=-0.035,cex.axis=0.75)
mtext(text="LD1",side=1,line=1.5)
axis(2,mgp=c(0,0.5,0),tck=-0.035,cex.axis=0.75)
mtext(text="LD2",side=2,line=1.5)
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[1],font=2,cex=1)
legend("bottomright",legend=paste0(tolower(levels(junco_data_m$Subspecies))," (",diag(as.matrix(CV_tab_m)),"/", apply(as.matrix(CV_tab_m),1,sum)," = ",round(diag(as.matrix(CV_tab_m))/apply(as.matrix(CV_tab_m),1,sum)*100,2),"%)"),pt.cex=0.75,cex=0.3,pch=21,pt.bg=paste0(rev(junco_cols),"80"), col=junco_cols,lwd=0,x.intersp=0,box.lwd=1,bg="transparent")

plot(junco_predict_f$x[,1]*-1,junco_predict_f$x[,2],axes=F,pch=21,bg=paste0(rev(junco_cols)[as.numeric(junco_data_f$Subspecies)],"80"),col= rev(junco_cols)[as.numeric(junco_lda_f$class)],lwd=2)
title(main=paste0("Female DFA (",round(sum(diag(as.matrix(CV_tab_f))) / sum(CV_tab_f) * 100,2),"%)"),cex.main=1)
box()
axis(1,mgp=c(0,0.5,0),tck=-0.035,cex.axis=0.75)
mtext(text="LD1",side=1,line=1.5)
axis(2,mgp=c(0,0.5,0),tck=-0.035,cex.axis=0.75)
mtext(text="LD2",side=2,line=1.5)
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[2],font=2,cex=1)
legend("bottomright",legend=paste0(tolower(levels(junco_data_f$Subspecies))," (",diag(as.matrix(CV_tab_f)),"/", apply(as.matrix(CV_tab_f),1,sum)," = ",round(diag(as.matrix(CV_tab_f))/apply(as.matrix(CV_tab_f),1,sum)*100,2),"%)"),pt.cex=0.75,cex=0.3,pch=21,pt.bg= paste0(rev(junco_cols),"80"),lwd=0,x.intersp=0,box.lwd=1,bg="transparent")

dev.off()

### Patten and Unitt 75% rule ###
patten_unitt<-function(x1,x2){
	mu1<-mean(x1)
	mu2<-mean(x2)

	s1<-sd(x1)
	s2<-sd(x2)

	if(mu1<mu2){
		t1<-mu1+s1*abs(qt(0.01, length(x1)-1))
		t2<-mu2-s2*abs(qt(0.25, length(x2)-1))
		d21<-t2-t1
		
		q1<-mu1+s1*abs(qt(0.25, length(x1)-1))
		q2<-mu2-s2*abs(qt(0.01, length(x2)-1))
		d12<-q2-q1
	}else{
		t1<-mu1-s1*abs(qt(0.25, length(x1)-1))
		t2<-mu2+s2*abs(qt(0.01, length(x2)-1))
		d12<-t1-t2
	
		q1<-mu1-s1*abs(qt(0.01, length(x1)-1))
		q2<-mu2+s2*abs(qt(0.25, length(x2)-1))
		d21<-q1-q2
	}
	
	if(d12 > 0){
		cat(paste0("D12 ≥ 0; Species 1 is diagnosable from Species 2 under the 75% rule"))
	}else{
		cat(paste0("D12 < 0; Species 1 is NOT diagnosable from Species 2 under the 75% rule"))
	}

	if(d21 > 0){
		cat(paste0("D21 ≥ 0; Species 2 is diagnosable from Species 1 under the 75% rule"))
	}else{
		cat(paste0("D21 < 0; Species 2 is NOT diagnosable from Species 1 under the 75% rule"))
	}
	
	return(list(d12=d12,d21=d21))
}

### Set up five pair-wise tests that we wish to conduct ###
test1<-c("oreganus","montanus")
test2<-c("oreganus","shufeldti")
test3<-c("montanus","shufeldti")
test4<-c("shufeldti","thurberi")
test5<-c("thurberi","pinosus")

all_tests<-list(test1,test2,test3,test4,test5)
test_output<-list()
for(i in 1:length(all_tests)){
	junco_m_trim<-junco_data_m[as.character(junco_data_m$Subspecies) %in% all_tests[[i]],]
	junco_m_trim_lda<-lda(Subspecies~.,data= junco_m_trim[c(6,c(15,16,24,18,25,26))],CV=F)
	junco_m_trim_predict<-predict(junco_m_trim_lda)
	
	x1<-junco_m_trim_predict$x[,1][as.character(junco_m_trim$Subspecies) == all_tests[[i]][1]]
	x2<-junco_m_trim_predict$x[,1][!as.character(junco_m_trim$Subspecies) == all_tests[[i]][1]]
	
	D_m<-patten_unitt(x1,x2)
	
	junco_f_trim<-junco_data_f[as.character(junco_data_f$Subspecies) %in% all_tests[[i]],]
	junco_f_trim_lda<-lda(Subspecies~.,data= junco_f_trim[c(6,c(15,16,24,18,25,26))],CV=F)
	junco_f_trim_predict<-predict(junco_f_trim_lda)
	
	x1<-junco_f_trim_predict$x[,1][as.character(junco_f_trim$Subspecies) == all_tests[[i]][1]]
	x2<-junco_f_trim_predict$x[,1][!as.character(junco_f_trim$Subspecies) == all_tests[[i]][1]]
	
	D_f<-patten_unitt(x1,x2)
	test_output[[i]]<-list(D_m,D_f)
}

sum_df<-cbind(t(data.frame(all_tests)),round(matrix(unlist(test_output),ncol=4,byrow=T),2))
rownames(sum_df)<-NULL
write.xlsx(sum_df,file="/Users/NickMason/Desktop/Manuscripts/Juncos/Tables/PattenUnittTable_v3.xlsx")