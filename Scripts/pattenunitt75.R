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