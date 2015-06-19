# Algorithm to calculate number of cosine values above 0.6 for a set of genes:

#a simulates similarity matrix
#x= matrix(runif(9000000,0,1),3000,3000);x=(x+t(x))/2;dim(x)
#write(t(x), file='B:/Academic/R/GCAT package/simMatrix.txt', sep='\t', ncolumns=3000)

x=t(as.matrix(read.table('B:/Academic/R/GCAT package/simMatrix.txt')));dim(x)

x[1:1000,1:1000]=0.61
#rownames(x)=c(1:200)
#colnames(x)=c(1:200)
#g=c(1:200)
  
#######################################
#counting the number of cosines greater than 0.6
#g: is the index of input gene list
#x: is the similarity matrix
countEdge=function(g,x){
  tmp=x[g,g]
  diag(tmp)=0
  #tmp1=tmp>0.6
  #n=(  length(tmp[tmp>0.6])-length( diag(tmp) [diag(tmp)>0.6] ) ) /2
  n=length(tmp[tmp>0.6])/2
  return(n)
}
g=sample(1:3000,100,replace=FALSE)
countEdge(g,x)
########################################
#re-sampling, we randomly choose r genes (length of vector g: number of the original sample genes)
#and re-sample from the cosine matrix m times
#and find # of cosine >0.6 and take average
nedge_boot=function(g,x){
  m=200
  l=length(g)
  tmp=matrix(0,m)
  for (i in 1:m){
    x1=sample(1:length(x[1,]),l,replace=TRUE)
    tmp[i]=countEdge(x1,x)
  }
  tmp1=floor(sum(tmp)/m)
  return(tmp1)
}
nedge_boot(g,x)

################################################################
#calculating Fisher exact test: literature p-value LPv
#m is number of re-sampling performed to calculate simulated LPV
Lpv=function(g,x){
  n1=countEdge(g,x)
  n2=nedge_boot(g,x)
  l=length(g)
  m=matrix(c(n1,n2,l^2/2-n1-l/2,l^2/2-n2-l/2),2,2,dimnames=list(c('TRUE','simulates'),c('> 0.6','< 0.6')))
  lpv=fisher.test(m,alternative='greater')
  return(lpv$p.value)
}
Lpv(g,x)
system.time(Lpv(g,x))

####################################################################################################
#function to obtain differentially expressed genes:
  # TO BE ADDED : how to do correction for multiple testing?
#####################################################################

# literature based functional significance:
# input: list of differentially expressed genes
#output: number of Lpvs <0.05 in the 1000 simulation of 50 genes from DEGs!!

N_lpvSimul=function(g,x){
  r=50
  lpvs=matrix(0,1000)
  for (i in 1:1000){
    g1=sample(g,r,replace=FALSE)
    lpvs[i]=Lpv(g1,x)
    }
  n1=length(lpvs[lpvs<= 0.05])
  return(n1)
  }
N_lpvSimul(g,x)
system.time(N_lpvSimul(g,x))
      #system time elapsed:   39.8 sec
##############################
#everything works up to this point.
#########################################################################
# d is the original gene expression data before test of significance: columns: different samples
#ind: index of group 1 data
#ind2: index of group 2 data

#alpha:level of significance to be used in t test!
#output "LCI *1000" which is number of Lpv's less than 'alpha' in 1000 re-sampling
#assume we perform student t test for now
#performs t.test for the total number of genes and returns the index of differentially expressed genes
        #it can be adjusted to return the name of the differentially expressed genes
#######################################################
#we need a function that performs at least 3 different test of significance in addition to t test
###################################
  #how to calculate many t.tests in parallel?? no need! it is fast enough!!!!
# find R packages that does diff expressed genes!!
      #we need to make it faster
##############################################################################
LCIOL=function(d,ind,alpha,x){
  g1=d[,ind];dim(g1)
  g2=d[,-ind];dim(g2)
  diffExpr=0
  for(j in 1:length(d[,1])){
    t1=t.test(g1[j,],g2[j,],altrnative='two.sided')
    if( t1$p.value < alpha ){
      diffExpr=c(diffExpr,j)
    }
  }
  n=N_lpvSimul(diffExpr[-1],x)
  return(n)
}

alpha=0.05
l1=LCIOL(d,ind,alpha,x);l1
system.time(LCIOL(d,ind,alpha,x))

#system time elapsed:   41.39 sec

##########################################################
#s: number of samples in group1
#permutes the sample labels 10 times 
#calculates p-value for each gene in the data set
#finds differentially expressed genes in each re-labelling
#calculates # of Lpv < alpha in each iteration
#and calculates avg # of Lpvs less than alpha in 10 relabelling (total/10)
s=27
LCIPL=function(d,s,alpha,x){
  l=0
  for (i in 1:10){
    ind=sample(length(d[1,]),s,replace=FALSE)
    g1=d[,ind];dim(g1)
    g2=d[,-ind];dim(g2)
    diffExpr=0
    for(j in 1:length(d[,1])){
      t1=t.test(g1[j,],g2[j,],altrnative='two.sided')
      if(t1$p.value <alpha){
        diffExpr=c(diffExpr,j)
       }
    }
    l=l+N_lpvSimul(diffExpr[-1],x)
    
  }
  
  avgLpv=floor(l/10)
  return(avgLpv)
}

LCIPL(d,s,alpha,x)
system.time(LCIPL(d,s,alpha,x))
#system time elapsed    406.09
################################################
#some simulated data
# d1=matrix(rnorm(30000,mean=3,sd=1),nrow=3000);dim(d1)
# d2=matrix(rnorm(30000,mean=0,sd=2),nrow=3000);dim(d2)
# d=cbind(d1,d2);dim(d);d[1,]
# s=10
# LCIPL(d,s,0.05,x)
# ind=c(1:10)
##################################################


###################################################

#calculating LBFS
LBFS=function(d,ind,alpha,x){
  r=50
  n1=LCIOL(d,ind,alpha,x)
  n2=LCIPL(d,length(ind),alpha,x)

  m=matrix(c(n1,n2,r^2/2-n1-r/2,r^2/2-n2-r/2),2,2,dimnames=list(c('TRUE','simulates'),c('> 0.6','< 0.6')))
  lpv=fisher.test(m,alternative='greater')
  pval=lpv$p.value
  return(pval)
}
alpha=0.04
LBFS(d,ind,alpha,x)
system.time(LBFS(d,ind,alpha,x))
#system time elapsed : 467.4 for 3000 genes in 38 samples!!!

#######################################################
#Literature aided statistical significance threshold: LAAST

  #1- specify an increasing sequence of EPv thresholds alpha=(a1,...,am)
  #2- for each DEGs generated at each alpha level, calculate LCI   (L1,...,Lm)

#LCI is calculated as the fraction of the sampled subsets that have LPv<0.05
        #output of N_lPVSimul divided by 1000 gives LCI
  #3-choose an integer m0=3 and perform two piece linear fit:
      #for k=m0,m0+1,...m-m0
          #fit straight line to points (aj,Lj), j=1,...,k and record beta0L, beta1L
          #fit a straight line to points(aj,Lj), j=k+1,...,m record beta0R , beta1R
          #vk=(beta0L-beta0R)+(beta1R-beta1L)
          #let K_star be the first local maxima of vk, k=m0,..,m-m0
          # ie: k_star=min{j: vj>vj+1}
          # take the k_star th entry on the alpha sequence specified as EPv threshold.

###############################################################################


LAASTLCI=function(d,ind,alpha,x){
  
  lci=matrix(0,length(alpha))
  for (i in 1:length(alpha)){
    lci[i]=LCIOL(d,ind,alpha[i],x)/1000
    
  }
  lci1=cbind(alpha,lci)
  return (lci1)
}
###########################
#now the final function:
LAAST=function(d,ind,x){
  alpha=seq(0.02,1, by=0.02);alpha
  lci=LAASTLCI(d,ind,alpha,x)
  p=length(lci[,1])-3
  tmp1=matrix(0,p,2)
  tmp2=matrix(0,p,2)
  
  for(k in 3:p){
    mod1=lci[1:k,]
    tmp1[k,]=lm(mod1[,2]~mod1[,1])$coefficients
    ll=k+1
    ul=length(lci[,1])
    mod2=lci[ll:ul,]
    tmp2[k,]=lm(mod2[,2]~mod2[,1])$coefficients
    
  }
  tmp=tmp1-tmp2;tmp=tmp[-c(1:2),]
  v=tmp[,1]-tmp[,2]
  vk=0
  i=2
  while(vk==0 && i<length(v)){
    
    if (v[i] >v[i-1]){
        vk=i
    }
    i=i+1
  }
  alph=alpha[vk]
  return(alph)
}
#####################################################################

EPthreshold=LAAST(d,ind,x)
system.time(LAAST(d,ind,x))
#system time elapsed: 2106.78


############################
#using Glob data set and artificial sililarity matrix:

d=as.matrix(read.table('B:/Academic/R/GCAT package/golubOrig-Normalized.txt'))[1:3000,];dim(d)
ind=1:27
mm=LAAST(d,ind,x)
system.time(LAAST(d,ind,x))

# Ls=lm(mm[,2]~mm[,1])
# m1=Ls$coefficients
# m1[2]


#g=matrix(1:12,3,4)
#g[,-2]