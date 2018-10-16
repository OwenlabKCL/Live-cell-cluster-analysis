source("internal_live.R")
foldernames=c("Movie 1")
MPR<<-77

###############################################
limit_high_th<<-500
limit_high_r<<-300
trust_in_radii<<-0.75 
noise_over_time<<-c()
radii_over_time<<-c()
number_of_clusters_over_time<<-c()
sapply(foldernames, function(foldername){
r = readLines(con=file.path(paste(foldername, "/config.txt", sep="")))
get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
model=get("model")
{if (model=="Gaussian(prec)"){
  xlim = as.v(get("xlim"))
  ylim = as.v(get("ylim"))
  histbins = as.v(get("histbins")) 
  histvalues = as.v(get("histvalues"))
  histbins<<-histbins
  histvalues_over_time<<-c(rep(0, length(histbins)))
  histvalues<<-histvalues
  if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
    useplabel=FALSE; pb=NULL; alpha=NULL
  }
  else {
    useplabel=TRUE; 
    pb=as.numeric(get("pbackground"))
    pb<<-pb
    alpha=as.numeric(get("alpha"))
  }
  if (length(grep("bestonly",r))==0) bestonly=FALSE
  else bestonly=as.numeric(get("bestonly"))>0
  if (length(grep("rseq",r))==0) rseq=seq(10, 200, by=5)
  else {
      rparams=as.v(get("rseq"))
      rseq=seq(rparams[1], rparams[2], by=rparams[3])
  }
  if (length(grep("thseq",r))==0) thseq=seq(5, 500, by=5)
  else {
      thparams=as.v(get("thseq"))
      thseq=seq(thparams[1], thparams[2], by=thparams[3])
  }
  if (length(grep("clustermethod",r))==0) clustermethod="K"
  else {
      method=as.numeric(get("clustermethod"))
      if (method==1) clustermethod="K"
      else clustermethod="DBSCAN"
  }
}
else {stop("Haven't implemented anything else!")}}

o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
#print("o")
#print(o)
f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value
psd <- function(sd){
  log(f(sd))-log(cst) 
}
minsd = histbins[1]; maxsd = histbins[length(histbins)]
minsd<<-minsd
maxsd<<-maxsd
#print("psd")
#print(psd)

ld=list.dirs(foldername, recursive=FALSE)
ld=ld[ld!=foldername]

nb_folder=length(ld)
u<<-0
cut_in_analysis<<-3
R_init<<-rseq[1]
TH_init<<-thseq[1]
R_by<<-rseq[2]-rseq[1]
TH_by<<-thseq[2]-thseq[1]
Trajectory_picked_proposalr<<-c()
Trajectory_picked_proposalth<<-c()


for (h in 1:nb_folder)
{ 
  data= read.csv(file.path(paste(foldernames,"/",h, "/data.txt", sep="")))
  S=dim(data)
  s=S[1]
  if (s>MPR){

    
    print(h)
    pts = data[,1:2]; sds = data[,3];
    
    
    if (u==1){

      res=Kclust(pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb, score=TRUE, rlabel=TRUE, report=TRUE, rseq=rseq, thseq=thseq, clustermethod=clustermethod)
      writeRes(res, file.path(paste(foldernames,"/",h,"/r_vs_thresh.txt", sep="")), file.path(paste(foldernames,"/",h, "/labels", sep="")), bestonly=bestonly)

      thsize=length(thseq)
      scores=res[["scores"]]
      a=length(thseq)
      b=length(rseq)
      score_mat=matrix(0,a,b)

      for (i in 1:b){
        for (j in 1:a){
          score_mat[j,i]=scores[1+(i-1)*a+(j-1)]
        }
      }

      MAX_index=which.max(scores)
      MAX_value=max(scores)
      #print("MAX_value")
      #print(MAX_value)
      th_with_max=c()
      r_with_max=c()
      for (i in 1:b){
        for (j in 1:a){
          if(score_mat[j,i]==MAX_value){
            th_with_max=c(th_with_max,j)
            r_with_max=c(r_with_max,i)
          }
        }
      }
      
      TH<<-0
      R<<-0
      
     R=R_init+(ceiling(mean(r_with_max))-1)*R_by #should be incrementation but I ll add it up later, I think it is done now
     R<<-R
     print("best R")
     print(R)
     R_best<<-R
     TH=TH_init+(ceiling(mean(th_with_max))-1)*R_by
     TH<<-TH
     print("best th")
     print(TH)
     TH_best<<-TH


       if (h<3){
          
           R=floor((Trajectory_picked_proposalr[1]+Trajectory_picked_proposalr[1]+R)/3)
           R<<-R
           TH=floor((Trajectory_picked_proposalth[1]+Trajectory_picked_proposalth[1]+TH)/3)
           TH<<-TH

         
       }else{
         #pb used over time is an average
         if (Trajectory_picked_proposalr[h-1]!=-1 && Trajectory_picked_proposalr[h-2]!=-1){
           R=floor((Trajectory_picked_proposalr[h-2]+Trajectory_picked_proposalr[h-1]+R)/3)
           R<<-R
           TH=floor((Trajectory_picked_proposalth[h-2]+Trajectory_picked_proposalth[h-1]+TH)/3)
           TH<<-TH
         }else{
           R<<-R
           TH<<-TH
         }

       }
       
       R<<-R
       TH<<-TH
       print("chosen TH")
       print(TH)
       
       print("chosen r")
       print(R)

 
      Trajectory_picked_proposalr<<-c(Trajectory_picked_proposalr, R) #kipping the average one
      Trajectory_picked_proposalth<<-c(Trajectory_picked_proposalth, TH)

      
      rseq=c()
      thseq=c()
      initr=0
      initth=0
      
      for (y in -3:3){
        
       
        e=R+y*R_by 
        if (e>0){
          rseq=c(rseq, e)

        }
        e=TH+y*TH_by
        if (e>0){
          thseq=c(thseq, e)   

        }
        
      }
      rseq<<-rseq
      thseq<<-thseq
      R_init<<-rseq[1]
      TH_init<<-thseq[1]
      v<<-0
     
      
      
    }

    if(u==0){
      res=Kclust(pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb, score=TRUE, rlabel=TRUE, report=TRUE, rseq=rseq, thseq=thseq, clustermethod=clustermethod)
      writeRes(res, file.path(paste(foldernames, "/",h,"/r_vs_thresh.txt", sep="")), file.path(paste(foldernames, "/",h,"/labels", sep="")), bestonly=bestonly)
      #best label for next
      thsize=length(thseq)
      scores=res[["scores"]]
      MAX_index=which.max(scores)
      a=length(thseq)
      b=length(rseq)
      score_mat=matrix(0,a,b)
     
      for (i in 1:b){
        for (j in 1:a){
          score_mat[j,i]=scores[1+(i-1)*a+(j-1)]
        }
      }

      MAX_index=which.max(scores)
      MAX_value=max(scores)
      th_with_max=c()
      r_with_max=c()
      for (i in 1:b){
        for (j in 1:a){
          if(score_mat[j,i]==MAX_value){
            th_with_max=c(th_with_max,j)
            r_with_max=c(r_with_max,i)
          }
        }
      }
      
      TH<<-0
      R<<-0
      R=R_init+(ceiling(mean(r_with_max))-1)*R_by #should be incrementation but I ll add it up later, I think it is done now
      R<<-R
      print("best R")
      print(R)
      TH=TH_init+(ceiling(mean(th_with_max))-1)*R_by
      TH<<-TH
      print("best TH")
      print(TH)
      
      R_best<<-R
      TH_best<<-TH
      Trajectory_picked_proposalr<<-c(Trajectory_picked_proposalr, R)
      Trajectory_picked_proposalth<<-c(Trajectory_picked_proposalth, TH)
      rseq=c()
      thseq=c()
      initr=0
      initth=0
      R_by<<-10
      TH_by<<-10
    
      for (y in -3:3){
        
        
        e=R+y*R_by
        if (e>0){
          rseq=c(rseq, e)

        }
        e=TH+y*TH_by 
        if (e>0){
          thseq=c(thseq, e)   

        }
        
      }
      rseq<<-rseq
      thseq<<-thseq
      R_init<<-rseq[1]
      TH_init<<-thseq[1]
      
      u<<-1
    }
    
   
    bfile=file.path(paste(foldername, "/",h,"/labels/clusterscale", R_best, " thresh", TH_best, "labels.txt", sep=""))
    nbfile=bfile
    labelsbest=c()
    labelsbest=strsplit(readLines(nbfile),",")[[1]]
    pb=1-as.numeric(percentageInCluster(labelsbest))/100
    noise_over_time<<-c(noise_over_time,pb)
    print("folder")
    print(h)
    print("cut_in_analysis")
    print(cut_in_analysis)
    print("detected noise")
    print(pb)
    pb<<-pb
    
    
    n_cluster=nClustersi(labelsbest)
    histvalues=c(rep(0, length(histbins)))
    number_of_clusters_over_time<<-c(number_of_clusters_over_time,length(n_cluster))
    radii_av=0
    
  
    if (h<3 || cut_in_analysis<3){
      if (h==1 || cut_in_analysis==1){
        pb=(0.5+0.5+pb)/3
        pb<<-pb
      }else{
        pb=(0.5+noise_over_time[length(noise_over_time)]+pb)/3
        pb<<-pb
      }
      
    }else{
    
      pb=(noise_over_time[length(noise_over_time)-1]+noise_over_time[length(noise_over_time)]+pb)/3
      pb<<-pb
    }
    
    
    if (pb>0.75){
      pb=0.75
      pb<<-pb
    }
    
    print("prior noise")
    print(pb)

    
    for (k in 1:length(n_cluster)){
      st=0
      l=n_cluster[k]
      H=0
      H=clusterRadiii(pts,labelsbest,l)
      radii_av=radii_av+H
      #print(H)
      for (j in 1:length(histbins)){
        #you have to go backwards
        if (st==0){
          if (H>histbins[length(histbins)-j+1]){
            st=1
            if (j==1){
              histvalues[length(histbins)-j+1]=histvalues[length(histbins)-j+1]+1
            }else{
              histvalues[length(histbins)-j+2]=histvalues[length(histbins)-j+2]+1 
            }
            
          }
        }
      }
      if (st==0){
        histvalues[1]=histvalues[1]+1
      }
    }
    
    radii_av=radii_av/length(n_cluster)
    radii_over_time<<-c(radii_over_time,radii_av)
    
    binsatzeros=0
    for (j in 1:length(histbins)){
      if (histvalues[j]==0){
        binsatzeros=binsatzeros+1
      }
    }
    tot=length(n_cluster)/trust_in_radii
    if (binsatzeros!=0){
      forbinsatzeros=((1-trust_in_radii)*tot)/binsatzeros
      for (j in 1:length(histbins)){
        if (histvalues[j]==0){
          histvalues[j]=forbinsatzeros
        }
      }
      
    }
    
    histvalues<<-histvalues
    
    print("detected histvalues")
    print(histvalues)
    
    if (h==1 || cut_in_analysis==1){
      histvalues_over_time=histvalues
      histvalues_over_time<<-histvalues_over_time
    }
    if (h==2 || cut_in_analysis==2){
      histvalues_over_time=rbind(histvalues_over_time,histvalues)
      histvalues_over_time<<-histvalues_over_time
      for (s in 1:dim(histvalues_over_time)[2]){
        histvalues[s]=mean(histvalues_over_time[,s])
      }
      histvalues<<-histvalues
      
    }
    if (h>2 || cut_in_analysis>2){
      histvalues_over_time=rbind(histvalues_over_time,histvalues)
      histvalues_over_time<<-histvalues_over_time
      for (s in 1:dim(histvalues_over_time)[2]){
        histvalues[s]=mean(histvalues_over_time[(dim(histvalues_over_time)[1]-2):dim(histvalues_over_time)[1],s])
      }
      histvalues<<-histvalues
    }
    
    print("prior histvalues")
    print(histvalues)
    
    o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
    f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
    cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value

    psd <- function(sd){
      log(f(sd))-log(cst) 
    }
    minsd = histbins[1]; maxsd = histbins[length(histbins)]
    minsd<<-minsd
    maxsd<<-maxsd
    cut_in_analysis=cut_in_analysis+1
    cut_in_analysis<<-cut_in_analysis

    

  }else{
    

    Trajectory_picked_proposalr<<-c(Trajectory_picked_proposalr, -1)
    Trajectory_picked_proposalth<<-c(Trajectory_picked_proposalth, -1)
    
    cut_in_analysis<<-1
    if (u==1){
      
      print("increasing scan")
      v<<-v+1
      r=-5-v
      rr=10+v

      rseq=c()
      thseq=c()
      initth=0
      initr=0
      for (y in r:rr){

        e=R+y*R_by 
        if (e>0 && e<limit_high_r){
          rseq=c(rseq, e)
          if(initr==0){
            R_init<<-e
            initr=1
          }
        }
        e=TH+y*TH_by 
        if (e>0  && e<limit_high_th){
          thseq=c(thseq, e) 
          if(initth==0){
            TH_init<<-e
            initth=1
          }
        }
        
      }
      rseq<<-rseq
      thseq<<-thseq
    
      
      r = readLines(con=file.path(paste(foldername, "/config.txt", sep="")))
      get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
      as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
      model=get("model")
      {if (model=="Gaussian(prec)"){
        histbins = as.v(get("histbins")) 
        histvalues = as.v(get("histvalues"))
        histbins<<-histbins
        histvalues<<-histvalues
        if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
          useplabel=FALSE; pb=NULL; alpha=NULL
        }
        else {
          useplabel=TRUE; 
          pb=as.numeric(get("pbackground"))
          pb<<-pb
         }
        }}
      o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
      f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
      cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value
      psd <- function(sd){
        log(f(sd))-log(cst) 
      }
      minsd = histbins[1]; maxsd = histbins[length(histbins)]
      minsd<<-minsd
      maxsd<<-maxsd

    }

  }
  
}
write.csv(Trajectory_picked_proposalr,file = paste(foldername, "/Radii_picked_for_analysis.csv", sep="/"), sep=",")
write.csv(Trajectory_picked_proposalth,file = paste(foldername, "/Threshold_picked_for_analysis.csv", sep="/"), sep=",")
write.csv(noise_over_time,file = paste(foldername, "/p_unclustered_loc_over_time.csv", sep="/"), sep=",")
write.csv(radii_over_time,file = paste(foldername, "/radii_over_time.csv", sep="/"), sep=",")
write.csv(number_of_clusters_over_time,file = paste(foldername, "/number_of_clusters_over_time.csv", sep="/"), sep=",")

})
