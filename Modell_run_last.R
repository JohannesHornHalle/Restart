setwd("C:/temp/Restart19") # setting of work directory; 
# Important Dataset "Start.Rdata" must be in this folder, if not adjust line 
# In this folder also the genererated Data is saved
library("foreach")
library("doSNOW")
library("dplyr")
#loading starting Parameter; the file was created by Modell_restart_preset.R
load("Start.Rdata")

# ll is the counter of the run; one run includes all scenarios each with the 4 different mass events settings 
for(ll in 1:10){
registerDoSNOW(makeCluster(4)) # depending on the number of your computers cores and RAM   

scens<-length(events)*length(cterm)*length(facinf)*length(aesn)-1
gc()
#loop over parameter set 
Erg<-foreach(rr = 0:scens) %dopar%{
fi<-rr%/%(length(cterm)*length(events)*length(aesn))+1
ev<-(rr%/%(length(cterm)*length(aesn)))%%(length(events))+1
ct<-(rr%/%length(aesn))%%(length(cterm))+1
ae<-rr%%(length(aesn))+1

# Inz3 is the summary Outcome matrix for each run; cols are the szenarios   
Inz3<-array(0,dim=c(4,14))
Inz3[,10]<-rep(0:3)
dimnames(Inz3)[[2]]<-c("ASY","LI","SC","H","ICU","Tot","Q","PI","PNI","Szen","event","mask","Inc","aesn")
Inz3[,11]<-as.integer(events[ev])
Inz3[,12]<-facinf[fi]
Inz3[,13]<-cterm[ct]
Inz3[,14]<-aesn[ae]

A<-c(rep(F_list[[1]][,"P1"],reploop[[1]]),rep(c(F_list[[2]][,"P1"],F_list[[2]][,"P2"]),rep(reploop[[2]],2)),rep(c(F_list[[3]][,"P1"],F_list[[3]][,"P2"],F_list[[3]][,"P3"]),rep(reploop[[3]],3)),rep(c(F_list[[4]][,"P1"],F_list[[4]][,"P2"],F_list[[4]][,"P3"],F_list[[4]][,"P4"]),rep(reploop[[4]],4)),rep(c(F_list[[5]][,"P1"],F_list[[5]][,"P2"],F_list[[5]][,"P3"],F_list[[5]][,"P4"],F_list[[5]][,"P5"]),rep(reploop[[5]],5)))
Npop<-length(A)

# Defining schools and workplaces based on age
Sch<-rep(0,Npop)
Arbeiter<-sum(A>4 & A<15)
Arbeiter[2:5]<-rmultinom(1,Arbeiter,as.numeric(BT[2,2:5]))
Klassen<-ceiling(sum(A<5)/classsize)
temp2<-c(Klassen,0,0,0,0)
temp2[2:5]<-ceiling(c(Arbeiter[2]/2,Arbeiter[3:5]/as.numeric(BT[3,3:5])))
temp<-c(rep((1:temp2[2])+temp2[1],2),rep((1:temp2[3])+sum(temp2[1:2]),ceiling(BT[3,3])),rep((1:temp2[4])+sum(temp2[1:3]),ceiling(BT[3,4])),rep((1:temp2[5])+sum(temp2[1:4]),ceiling(BT[3,5])))
temp<-temp[order(runif(length(temp)))][1:Arbeiter[1]]
#rm(temp2)
Sch[A>4 & A<15]<-temp
#rm(temp)
Lehrer<-((1:Npop)[A>4 & A<15])[sample(sum(A>4 & A<15),Klassen)]
Sch[A<5]<-rep(1:Klassen,each=classsize)
Sch[Lehrer] <- (1:Klassen)
Sch[Sch==0]<-NA
Sch2<-Sch
Sch3<-Sch
for(i in 1:max(Sch,na.rm = T)){
  wbool<-which(Sch==i) 
  Sch2[wbool]<-length(wbool)
  Sch3[wbool]<-1:length(wbool)
}
Sch[is.na(Sch)==F]<-as.numeric(factor(Sch[is.na(Sch)==F]))

Dat<-array(0,dim=c(Npop,8))
Dat[,1]<-A
Dat[,2]<-fam2
Dat[,3]<-fam
Dat[,4]<-c(rep(1,sum(fam2==1)),rep(1:2,sum(fam2==2)/2),rep(1:3,sum(fam2==3)/3),rep(1:4,sum(fam2==4)/4),rep(1:5,sum(fam2==5)/5))
Dat[,5]<-Sch2
Dat[,6]<-Sch
Dat[,7]<-Sch3
Dat[,8]<-1:Npop
dimnames(Dat)[[2]]<-c("A","HS","Fam","IF","HW","Work","IW","I")
#rm(Sch)
#rm(Sch2)
#rm(Sch3)
#rm(fam)
#rm(fam2)

#Description of Dat (Contains all information about contact settings) 
# 1.	A    Age
# 2.	HS   Houshold size
# 3.	Fam  Indexnumber of Houshold/family
# 4.	IF   Indexnumber within Houshold/family from 1 to HS
# 5.	HW   Number of persons in workplace/school
# 6.	Work Indexnumber of workplace/school 
# 7.	IW   Indexnumber within workplace/school 1 bis HW
# 8.	I    general Indexnumber

N_fam<-max(Dat[,3])
N_work<-max(Dat[,6],na.rm = T)
N_workm<-max(Dat[,5],na.rm = T)

# Lists describing the Population decreases model runtime 
# Flist Family list
Flist<-array(as.numeric(NA),dim=c(N_fam,7))
Flist[,c(1,6,7)]<-Dat[Dat[,4]==1,c(8,2,3)]
Flist<-Flist[order(Flist[,7]),]  # eigentlich unnötig denn das sollte schon geordnet sein 
for(i in 2:5){
  ID<-Dat[Dat[,4]==i,3]
  Flist[ID,i]<-Dat[Dat[,4]==i,8]
}

# Wlist workplace list
Wlist<-array(as.numeric(NA),dim=c(N_work,N_workm+2))
Wlist[,c(1,1:2+N_workm)]<-Dat[Dat[,7]==1 & is.na(Dat[,7])==F,c(8,5,6)]
Wlist<-Wlist[order(Wlist[,2+N_workm]),]
for(i in 2:N_workm){
  ID<-Dat[Dat[,7]==i & is.na(Dat[,7])==F,6]
  Wlist[ID,i]<-Dat[Dat[,7]==i & is.na(Dat[,7])==F,8]
}

# Rlist Random list/ age list
Rlist<-list()
for(i in 1:18){
  Rlist[[i]]<-which(Dat[,1]==i)}

#a=asymptomatic course of disease, l=mild course of disease, s=severe course of disease
boola<-runif(Npop)<per_s1[A]
booll<- rep(F,Npop)
booll[!boola]<-runif(sum(!boola))<per_s2[A[!boola]]/(1-per_s1[A[!boola]])
bools<- (!(boola | booll))

#Hospitalisation, ICU und death, precisly 1 -
boolh<-runif(Npop)<per_h[A]
booli<-runif(Npop)<per_i[A]  
boolm<-runif(Npop)<per_m[A]

asp<-(runif(Npop)<asy) +1 # asp defines if a person when infected will later get symptoms(value 1) or remain completely asymptotic (value 2)
#setting starting condition
S<-rep(0,Npop)     # Status 0 susceptible; 1 latent ; 2 pre-symptomatic; 3 symptomatic; 7 resistant only current status is saved; 4 hospital; 5 ICU; 6 Tot 
S[sample(1:Npop,round((FL[2]-FL[1])*DZ))]<-7 # resistant persons
S[sample(1:Npop,round(cterm[ct]/12*FL[1]*DZ*(ZZ1[1]*ZZ2[1]+1)/sum(ZZ1[1:3]*ZZ2[1:3]+1)))]<-1     # asymptomatic persons
S[sample(1:Npop,round(cterm[ct]/12*FL[1]*DZ*(ZZ1[2]*ZZ2[2]+1)/sum(ZZ1[1:3]*ZZ2[1:3]+1)))]<-2     # presymptomatic
S[sample(1:Npop,round(cterm[ct]/12*FL[1]*DZ*(ZZ1[3]*ZZ2[3]+1)/sum(ZZ1[1:3]*ZZ2[1:3]+1)))]<-3     # symptomatic
S[sample(1:Npop,FL[5])]<-4                   # hospital
Z<-rep(-1,Npop)                              # Day counter S, R
bool<-S>0 & S<7
Z[bool]<-rbinom(sum(bool),ZZ1[S[bool]],ZZ2[S[bool]])+1 # Day counter other stages

# disease stages
# 0 susceptibel
# 1 latent
# 2 Presymptomatic infectious
# 3 Infectious symptomatic
# 4 Hospital
# 5 ICU 
# 6 Dead
# 7 Resistent 

ZI<-rep(-100,Npop)    # Day counter when a Person becomes infected
ZS<-rep(-100,Npop)    # Day Counter when a Person becomes symptomatic 

#tps<-seq(0.8,1,length=Days) # probability to be detected changes/improvement in time symptomatic
#tpa<-seq(0.3,1,length=Days) # probability to be detected changes/improvement in time asymptomatic  
#apppro<-seq(0,0.7,length=Days) # proportion with app
#dla<-c(0.9,3) # effectiveness app, delay detection possible 
con<-c(2.5,2.5,2.5)
INFEK<-rep(0,Npop)       # Day of infection
Symp<-rep(-10,Npop)      # Day of getting Symptoms
Symp[S==3]<-  -1
Detect<-rep(-1,Npop)     # Day of getting detected
Qu<-rep(-10,Npop)         # Days person will be still in quarantine (only temporal since person can get quarantained twice)  

#defining Matrices
Qul<-rep(0,Days+burnindays)
Det<-array(0,dim=c(6,Days+burnindays))
Out<-array(0,dim=c(Days+burnindays,9))
dimnames(Out)[[2]]<-c("S","E","PS","IS","H","ICU","D","R","A")
Out1<-rep(0,Days+burnindays)

#Contact lists used for tracking subdivided in H household,W workplace, R random and also A asymptomatic S symptomatic  
CHA<-list()
CHS<-list()
CWA<-list()
CWS<-list()
CRA<-list()
CRS<-list()
for(i in 1:(Days+burnindays)){
  CHA[[i]]<-array(numeric(),dim=c(0,7))
  CHS[[i]]<-array(numeric(),dim=c(0,7))
  CWA[[i]]<-array(numeric(),dim=c(0,7))
  CWS[[i]]<-array(numeric(),dim=c(0,7))
  CRA[[i]]<-array(numeric(),dim=c(0,5))
  CRS[[i]]<-array(numeric(),dim=c(0,5)) 
  
  dimnames(CHA[[i]])[[2]]<-c("P","HS","Fam","IF","IC","C","INF")
  dimnames(CHS[[i]])[[2]]<-c("P","HS","Fam","IF","IC","C","INF")
  dimnames(CWA[[i]])[[2]]<-c("P","WS","WG","IW","IC","C","INF")
  dimnames(CWS[[i]])[[2]]<-c("P","WS","WG","IW","IC","C","INF")
  dimnames(CRA[[i]])[[2]]<-c("P","AP","AC","C","INF")
  dimnames(CRS[[i]])[[2]]<-c("P","AP","AC","C","INF")
}
evlist<-list()
for(i in 1:(Days+burnindays)){
evlist[[i]]<-array(numeric(),dim=c(0,3))
dimnames(evlist[[1]])[[2]]<-c("P","C","INF")}

for(Szen in 0:4){
Inz<-array(0,dim=c(Days+burnindays,9))
#dimnames(Inz)[[1]]<-c("Burnin","Baseline","Szen1","Szen2","Szen3")
dimnames(Inz)[[2]]<-c("ASY","LI","SC","H","ICU","Tot","Q","PI","PNI") 

if(Szen==0){Daysr<-c(1:burnindays)}else{Daysr<-c((1+burnindays):(Days+burnindays))}
if(Szen>0){
  Z<-Zsave
  S<-Ssave
  Qu<-Qusave
  Qul[(1+burnindays):(Days+burnindays)]<-0
  Symp<-Sympsave
  Detect[Detect>burnindays]<- -1
  INFEK[INFEK>burnindays]<- 0
  Det[,(1+burnindays):(Days+burnindays)]<-0
  for(i in (burnindays+1):(Days+burnindays)){
    CHA[[i]]<-array(numeric(),dim=c(0,7))
    CHS[[i]]<-array(numeric(),dim=c(0,7))
    CWA[[i]]<-array(numeric(),dim=c(0,7))
    CWS[[i]]<-array(numeric(),dim=c(0,7))
    CRA[[i]]<-array(numeric(),dim=c(0,5))
    CRS[[i]]<-array(numeric(),dim=c(0,5)) 
    
    dimnames(CHA[[i]])[[2]]<-c("P","HS","Fam","IF","IC","C","INF")
    dimnames(CHS[[i]])[[2]]<-c("P","HS","Fam","IF","IC","C","INF")
    dimnames(CWA[[i]])[[2]]<-c("P","WS","WG","IW","IC","C","INF")
    dimnames(CWS[[i]])[[2]]<-c("P","WS","WG","IW","IC","C","INF")
    dimnames(CRA[[i]])[[2]]<-c("P","AP","AC","C","INF")
    dimnames(CRS[[i]])[[2]]<-c("P","AP","AC","C","INF")
  }
}

#Schleife
for(dd in Daysr){ 
  # Persons under Quarantaine
  # 1. Symptomatic detction 
  # severe course then light course than asymptomatic 
  temp<-which(Symp == dd- detects[1,1] & bools)
  temp<-temp[runif(length(temp))<detects[1,2]]
  temp<-temp[Detect[temp]== -1]
  Qu[temp]<-detects[1,3]
  Detect[temp]<-dd
  Det[1,dd]<-length(temp)
  temp<-which(Symp == dd - detects[2,1] & booll)
  temp<-temp[runif(length(temp))<detects[2,2]]
  temp<-temp[Detect[temp]== -1]
  Qu[temp]<-detects[2,3]
  Detect[temp]<-dd
  Det[2,dd]<-length(temp)
  temp<-which(Symp == dd - detects[3,1] & boola)
  temp<-temp[runif(length(temp))<detects[3,2]]
  temp<-temp[Detect[temp]== -1]
  Qu[temp]<-detects[3,3]
  Detect[temp]<-dd
  Det[3,dd]<-length(temp)
  Inz[dd,7]<-sum(Det[1:3,dd])
  
  # 2.Detection via contact tracing
  # household contact, work contact, random contact  
  if(dd>detectc[1,1]+ 1){
  temp<-which(Detect+detectc[1,1]==dd)
  if(length(temp)>0){
  tempx<-unname(Dat[temp,2])
  tempy<-unname(Dat[temp,3])
  temp2<-c(unname(Flist[tempy,1]),unname(Flist[tempy,2]),unname(Flist[tempy,3]),unname(Flist[tempy,4]),unname(Flist[tempy,5]))
  temp2<-temp2[!(temp2 %in% temp) & is.na(temp2)==F]
  temp2<-temp2[Detect[temp2]== -1]
  Qu[temp2]<-detectc[1,3]-detectc[1,1]
  Inz[dd,7]<-Inz[dd,7]+length(temp2)
  temp2<-temp2[S[temp2] %in% c(2,3)]
  Detect[temp2]<-dd
  Det[4,dd]<-length(temp2)
  }}
  if(dd>detectc[2,1]+ 1){
  temp<-which(Detect+detectc[2,1]==dd)
  if(length(temp)>0){
  vec<-max(dd-detectc[2,1]-6,1):(dd-detectc[2,1])
  temp2<-c(unname(CWA[[vec[1]]][CWA[[vec[1]]][,1] %in% temp,6]),unname(CWS[[vec[1]]][CWS[[vec[1]]][,1] %in% temp,6]))
  if(length(vec)>1){
    for(i in 2:length(vec)){
      temp2<-c(temp2,unname(CWA[[vec[i]]][CWA[[vec[i]]][,1] %in% temp,6]),unname(CWS[[vec[i]]][CWS[[vec[i]]][,1] %in% temp,6]))
    }}
  temp2<-unique(temp2[runif(length(temp2))<detectc[2,2]])
  temp2<-temp2[Detect[temp2]== -1]
  Qu[temp2]<-detectc[2,3]-detectc[2,1]
  Inz[dd,7]<-Inz[dd,7]+length(temp2)
  temp2<-temp2[S[temp2] %in% c(2,3)]
  Detect[temp2]<-dd
  Det[5,dd]<-length(temp2)
  }}
  if(dd>detectc[3,1]+ 1){
  temp<-which(Detect+detectc[3,1]==dd)
  if(length(temp)>0){
  vec<-max(dd-detectc[3,1]-6,1):(dd-detectc[3,1])
  temp2<-c(unname(CRA[[vec[1]]][CRA[[vec[1]]][,1] %in% temp,4]),unname(CRS[[vec[1]]][CRS[[vec[1]]][,1] %in% temp,4]))
  if(length(vec)>1){
    for(i in 2:length(vec)){
      temp2<-c(temp2,unname(CRA[[vec[i]]][CRA[[vec[i]]][,1] %in% temp,4]),unname(CRS[[vec[i]]][CRS[[vec[i]]][,1] %in% temp,4]))
    }}
  temp2<-unique(temp2[runif(length(temp2))<detectc[3,2]])
  temp2<-temp2[Detect[temp2]== -1]
  Qu[temp2]<-detectc[3,3]-detectc[3,1]
  Inz[dd,7]<-Inz[dd,7]+length(temp2)
  temp2<-temp2[S[temp2] %in% c(2,3)]
  Detect[temp2]<-dd
  Det[6,dd]<-length(temp2)
  }}
  Qul[dd]<-sum(Qu> -1)
  
  # contacts  
  c1<-which(S==2)
  c2<-which(S==3)
  c1q<-c1[Qu[c1]== -10] # removing contacts which are in quarantine
  c2q<-c2[Qu[c2]== -10] # removing contacts which are in quarantine
  #ss<-(S==0)
  
  if(length(c1)>1){
    #presymptomatic Household
    c1x<-c1[Dat[c1,2]>1]
    CHAP<-rep(c1x,rpois(length(c1x),rco[A[c1x]+(Dat[c1x,2]-1)*acl,1]))
    m<-length(CHAP)
    if(m>0){
    CHAPI<-array(Dat[CHAP,c(2:4)],dim=c(length(CHAP),3),dimnames = list(NULL,dimnames(Dat)[[2]][2:4]))
    CHAIC<-ceiling(runif(m)*(Dat[CHAP,2]-1))
    if(m>1){
    bool<-CHAIC==CHAPI[,"IF"]
    CHAIC[bool]<-CHAPI[bool,1]}else{
    bool<-CHAIC==CHAPI["IF"] & is.na(CHAIC)==F
    CHAIC[bool]<-CHAPI[bool & c(T,F,F)]}
    CHAC<-rep(0,m)
    for(i in 1:m){if(CHAIC[i]>0){
    CHAC[i]<-Flist[CHAPI[i,2],][CHAIC[i]]}}
    CHAINF<-runif(m)<as.numeric(pinf[1,2-boola[CHAP]])
    CHA[[dd]]<-cbind(CHAP,CHAPI,CHAIC,CHAC,CHAINF)}
    
    #presymptomatic Workplace
    CWAP<-rep(c1q,rpois(length(c1q),rco[A[c1q]+(Dat[c1q,2]-1)*acl,2]))
    m<-length(CWAP)
    if(m>0){
    CWAPI<-array(Dat[CWAP,c(5:7)],dim=c(length(CWAP),3),dimnames = list(NULL,dimnames(Dat)[[2]][5:7]))
    CWAIC<-ceiling(runif(m)*(Dat[CWAP,5]-1))
    if(m>1){
    bool<-CWAIC==CWAPI[,"IW"]
    CWAIC[bool]<-CWAPI[bool,1]}else{
    bool<-CWAIC==CWAPI["IW"] & is.na(CWAIC)==F
    CWAIC[bool]<-CWAPI[bool & c(T,F,F)]}
    CWAC<-rep(0,m)
    for(i in 1:m){if(CWAIC[i]>0){
      CWAC[i]<-Wlist[CWAPI[i,2],][CWAIC[i]]}}
    CWAINF<-runif(m)<as.numeric(pinf[1,2-boola[CWAP]])
    CWA[[dd]]<-cbind(CWAP,CWAPI,CWAIC,CWAC,CWAINF)}
    
    #presymptotic random contacts
    CRA5<-rep(c1q,rpois(length(c1q),rco[A[c1q]+(Dat[c1q,2]-1)*acl,3]))
    m<-length(CRA5)
    if(m>0){
    CRA1<-Dat[CRA5,1]
    CRA2<-CRA1
    CRA3<-CRA5
    for(i in 1:acl){
      bool<-CRA1==i
      CRA2[bool]<-apply(rmultinom(sum(bool),1,prob = prop[i,])==1,2,which)}
    for(i in 1:acl){
      bool<-CRA2==i
      CRA3[bool]<-sample(Rlist[[i]],sum(bool),replace = T)}
    CRA4<-runif(m)<as.numeric(pinf[1,2-boola[CRA5]])
    CRA[[dd]]<-cbind(CRA5,CRA1,CRA2,CRA3,CRA4)}
  }
  
  if(length(c2)>1){
    #symptomatic Household
    c2x<-c2[Dat[c2,2]>1]
    CHSP<-rep(c2x,rpois(length(c2x),rco[A[c2x]+(Dat[c2x,2]-1)*acl,1]))
    m<-length(CHSP)
    if(m>0){
    CHSPI<-array(Dat[CHSP,c(2:4)],dim=c(length(CHSP),3),dimnames = list(NULL,dimnames(Dat)[[2]][2:4]))
    CHSIC<-ceiling(runif(m)*(Dat[CHSP,2]-1))
    if(m>1){
    bool<-CHSIC==CHSPI[,"IF"]
    CHSIC[bool]<-CHSPI[bool,1]}else{
    bool<-CHSIC==CHSPI["IF"] & is.na(CHSIC)==F
    CHSIC[bool]<-CHSPI[bool & c(T,F,F)]}
    CHSC<-rep(0,m)
    for(i in 1:m){if(CHSIC[i]>0){
      CHSC[i]<-Flist[CHSPI[i,2],][CHSIC[i]]}}
    CHSINF<-runif(m)<as.numeric(pinf[1,2-boola[CHSP]])
    CHS[[dd]]<-cbind(CHSP,CHSPI,CHSIC,CHSC,CHSINF)}
    
    #symptomatic Workplace
    CWSP<-rep(c2q,rpois(length(c2q),rco[A[c2q]+(Dat[c2q,2]-1)*acl,2]))
    m<-length(CWSP)
    if(m>0){
    CWSPI<-array(Dat[CWSP,c(5:7)],dim=c(length(CWSP),3),dimnames = list(NULL,dimnames(Dat)[[2]][5:7]))
    CWSIC<-ceiling(runif(m)*(Dat[CWSP,5]-1))
    if(m>1){
    bool<-CWSIC==CWSPI[,"IW"] & is.na(CWSIC)==F
    CWSIC[bool]<-CWSPI[bool,"HW"]}else{
    bool<-CWSIC==CWSPI["IW"] & is.na(CWSIC)==F
    CWSIC[bool]<-CWSPI[bool & c(T,F,F)]}
    CWSC<-rep(0,m)
    for(i in 1:m){if(CWSIC[i]>0){
    CWSC[i]<-Wlist[CWSPI[i,2],][CWSIC[i]]}}
    CWSINF<-runif(m)<as.numeric(pinf[1,2-boola[CWSP]])
    CWS[[dd]]<-cbind(CWSP,CWSPI,CWSIC,CWSC,CWSINF)}
    
    #symptotic random contacts
    CRS5<-rep(c2q,rpois(length(c2q),rco[A[c2q]+(Dat[c2q,2]-1)*acl,3]))
    m<-length(CRS5)
    if(m>0){
    CRS1<-Dat[CRS5,1]
    CRS2<-CRS1
    CRS3<-CRS5
    for(i in 1:acl){
      bool<-CRS1==i
      CRS2[bool]<-apply(rmultinom(sum(bool),1,prob = prop[i,])==1,2,which)}
    for(i in 1:acl){
      bool<-CRS2==i
      CRS3[bool]<-sample(Rlist[[i]],sum(bool),replace = T)}
    CRS4<-runif(m)<as.numeric(pinf[1,2-boola[CRS5]])
    CRS[[dd]]<-cbind(CRS5,CRS1,CRS2,CRS3,CRS4)}
  }
  
  if(Szen>1){
  # exclusion of persons younger than 15 years and persons in quarantaine 
  temp<-sample(Dat[A>3 & Qu== -10,8],round(events[ev]/30))
  tempx<-temp[S[temp]==2]
  tempy<-temp[S[temp]==3]
  Inz[dd,8]<-length(tempx)+length(tempy)
  m1<-rep(tempx,rpois(length(tempx),eventsc[Szen-1]))
  m2<-rep(tempy,rpois(length(tempy),eventsc[Szen-1]))
  m3<-rep(tempx,rpois(length(tempx),ac[Szen-1,aesn[ae]]))
  m4<-rep(tempy,rpois(length(tempy),ac[Szen-1,aesn[ae]]))
  tempz<-boola[c(m1,m3,m2,m4)]
  m<-rep(0,3)
  m[1]<-length(m1)+length(m2)+length(m3)+length(m4)
  m[2]<-length(m1)+length(m3)
  m[3]<-length(m2)+length(m4)
  NInfe<-numeric(0)
  if(sum(m)>0){
  evlist[[dd]]<-matrix(c(m1,m3,m2,m4,rep(0,2*m[1])),nrow=m[1],ncol=3)
  evlist[[dd]][,2]<-sample(temp[1:(m[1]-1)],m[1],replace = T)
  evlist[[dd]][evlist[[dd]][,1]==evlist[[dd]][,2],2]<-temp[m[1]]
  bool<-c(rep(T,length(m1)),rep(T,length(m3)),rep(F,length(m2)),rep(F,length(m4)))
  if(m[2]>0){evlist[[dd]][bool,3]<-as.numeric(runif(sum(bool))<facinf[fi]*unname(pinf[1,2-tempz[1:m[2]]]))}
  if(m[3]>0){evlist[[dd]][!bool,3]<-as.numeric(runif(sum(!bool))<facinf[fi]*unname(pinf[2,2-tempz[(m[2]+1):m[1]]]))}
  NInfe<-unname(evlist[[dd]][evlist[[dd]][,3]==1,2])
  NInfe<-NInfe[S[NInfe]==0]
  Inz[dd,9]<-length(NInfe)
  }}
    
  if(length(c1)>0 & length(c2)>0){
    NInf<-sort(unique(c(CHA[[dd]][CHA[[dd]][,7]==1,6],CHS[[dd]][CHS[[dd]][,7]==1,6],CWA[[dd]][CWA[[dd]][,7]==1,6],CWS[[dd]][CWS[[dd]][,7]==1,6],CRA[[dd]][CRA[[dd]][,5]==1,4],CRS[[dd]][CRS[[dd]][,5]==1,4])))
  }else{
    if(length(c1)==0){NInf<-sort(unique(c(CHS[[dd]][CHS[[dd]][,7]==1,6],CWS[[dd]][CWS[[dd]][,7]==1,6],CRS[[dd]][CRS[[dd]][,5]==1,4])))}
    if(length(c2)==0){NInf<-sort(unique(c(CHA[[dd]][CHA[[dd]][,7]==1,6],CWA[[dd]][CWA[[dd]][,7]==1,6],CRA[[dd]][CRA[[dd]][,5]==1,4])))}  
  }
  
  if(Szen>1){
  if(length(NInfe)>0){NInf<-unique(c(NInf,NInfe))}  
  }
  
  if(length(c1)>0 | length(c2)>0){
    NInf<-NInf[S[NInf]==0]
    NInf<-NInf[NInf>0 & is.na(NInf)==F]
    INFEK[NInf]<-dd
    #New Infections
    S[NInf]<-1
    Z[NInf]<-rbinom(length(NInf),ZZ1[1],ZZ2[1])+1
    Inz[dd,1]<-sum(boola[NInf])
    Inz[dd,2]<-sum(booll[NInf])
    Inz[dd,3]<-sum(bools[NInf])
  }
  
  bool<-(Z==0)
  if(sum(bool)>0){
  Symp[bool & S==2]<-dd
  Inz[dd,4]<-sum(bool & S==3 & bools & (boolh))
  Inz[dd,5]<-sum(bool & S==4 & (booli))
  Inz[dd,6]<-sum(bool & S==5 & (boolm))
  bool2<- bool & ((S==3 & !(bools)) | (S==3 & bools & (!boolh))   | (S==4 & (!booli))   |  (S==5 & (!boolm)))
  bool <- bool & (!bool2)
  Z[bool]<-rbinom(sum(bool),ZZ1[S[bool]+1],ZZ2[S[bool]+1])+1
  S[bool]<-S[bool]+1
  S[bool2]<-  7
  bool3<-bool & (S==3 & !(bools & boolh)) | (S==4 & (!booli)) | (S==5  & (!boolm))
  Z[bool3]<-rbinom(sum(bool3),ZZ3[S[bool3]-2],ZZ4[S[bool3]-2])+1
  Z[S > 5]<- -1}
  bool<-(Qu> -1)
  Qu[bool]<-Qu[bool]-1
  Qu[!bool]<- -10
  Z <- Z-1
  Z[Z < -1]<- -1
for(kk in c(0,4,5,6,7)){Out[dd,kk+1]<-sum(S == kk)}
for(kk in c(1,2,3)){
Out[dd,kk+1]<-sum(S == kk & (!boola))}
Out[dd,9]<-sum(S %in% c(1,2,3) & (boola))
Out1[dd]<-length(NInf)
}
if(Szen==0){
  Zsave<-Z
  Ssave<-S
  Qusave<-Qu
  Sympsave<-Symp}else{Inz3[Szen,1:9]<-apply(Inz,2,sum)}
}
return(Inz3)}

for(i in 1:length(Erg)){names(Erg[[i]])<-dimnames(Erg[[i]])[2]}
Erg2<-bind_rows(Erg)
write.csv(Erg2,paste("Frun",ll,".csv",sep=""))
}
