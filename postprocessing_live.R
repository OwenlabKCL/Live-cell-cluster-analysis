source("internal_live.R")
foldernames=c("Movie 1")
Starting_frame_id<<-2304 
nb_folder<<-20
MPR<<-77
######################################

makeplot=TRUE
superplot=FALSE
skeleton=FALSE
start_analysis<<-0

contribution=rep(0,11)


faith=rep(0,11)



merging=c(0,0)



splitting=c(0,0)


contribution<<-contribution
faith<<-faith
merging<<-merging
splitting<<-splitting

###########################################################################

first=0


colors=sample(rainbow(1000))

sapply(foldernames, function(expname){ #In the case you put files to post processin serie
  
  #################################################################
  #########################Settings################################
    nexpname=expname
    if (skeleton){
        nexpname=paste("R_", expname, sep="")
        dir.create(file.path(nexpname))
        file.copy(file.path(paste(expname, "/config.txt", sep="")), paste(nexpname, "/config.txt", sep=""))
        file.copy(file.path(paste(expname, "/sim_params.txt", sep="")), paste(nexpname, "/sim_params.txt", sep=""))
    }
    r = readLines(con=file.path(paste(nexpname, "/config.txt", sep="")))
    get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
    as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
    model=get("model")
    {if (model=="Gaussian(prec)"){
        xlim = as.v(get("xlim"))
        ylim = as.v(get("ylim"))
        histbins = as.v(get("histbins"))
        histvalues = as.v(get("histvalues"))
        if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
            useplabel=FALSE; pb=NULL; alpha=NULL
        }
        else {
            useplabel=TRUE; 
            pb=as.numeric(get("pbackground"))
            alpha=as.numeric(get("alpha"))
        }
    }
     else {stop("Haven't implemented anything else!")}}
    
    all=list.files(expname)
    dirnames=all[file.info(file.path(paste(expname, "/", all, sep="")))$isdir]
    axes=FALSE;
    cex=1/(3*sqrt(length(dirnames)))
    if (makeplot & superplot) {
        ##These settings control image size and resolution
        png(file.path(paste(nexpname, "/together.png",sep="")), width=10, height=10, units="cm", res=1200)
        nrow=ceiling(sqrt(length(dirnames)))
        par(mfrow=c(nrow, nrow))
    }
    
    beenthereonce<<-0 
    coming_from_break<<-1 

    ###############################################################################
    ###############################################################################
    ##############Analysis per point resolution frame starts here##################
    
      for (dirname in 1:nb_folder){
        print(dirname)
        ###loading the files required for the analysis
        foldername=file.path(paste(expname, "/", dirname, sep=""))
        nfoldername=file.path(paste(nexpname, "/", dirname, sep=""))
        if (skeleton){
            dir.create(nfoldername)  
            file.copy(file.path(paste(foldername, "/data.txt", sep="")), file.path(paste(nfoldername, "/data.txt", sep="")))
        }
        data=read.csv(file.path(paste(nfoldername, "/data.txt", sep="")))
        pts = data[,1:2]; sds = data[,3]; frame=data[,4];
        ###
        
        ###The analysis continues only for files that have been analysed ie >point_resolution 
        S=dim(data)
        s=S[1]
        s<<-s
        if (s>MPR){
        beenthereonce<<-1
        
        ###Loading the best scores and localiations labels
        if (skeleton){
            file.copy(file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(paste(nfoldername, "/r_vs_thresh.txt", sep="")))
        }
        r = read.csv(file.path(paste(nfoldername, "/r_vs_thresh.txt",sep="")), header=FALSE, sep="\t")
        m = as.matrix(r)
        cs=(m[1,])[-1]
        thr=(m[,1])[-1]
        m = as.matrix(m[2:length(m[,1]),2:length(m[1,])])
        which.maxm <- function(mat){
            indcol <- rep(1:ncol(mat), each=nrow(mat))[which.max(mat)] 
            indrow <- rep(1:nrow(mat), ncol(mat))[which.max(mat)]
            c(indrow, indcol)
        }
        best=which.maxm(m)
        bestcs=cs[best[2]]
        bestcs<<-bestcs 
        bestthr=thr[best[1]]
        bestthr<<-bestthr
        bfile=file.path(paste(foldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
        nbfile=bfile
        if (skeleton){
            dir.create(paste(nfoldername, "/labels", sep=""))    
            nbfile=file.path(paste(nfoldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
            file.copy(bfile, nbfile)
        }
        ###
        
        
        ################################################################################
        #######################Initialisation on the first data#########################
        if (start_analysis==0){
          Starting_frame_id=frame[1]
          Starting_frame_id<<-Starting_frame_id
          start_analysis<<-1
          labelsbest=c()
          labelsbest=strsplit(readLines(nbfile),",")[[1]]
          labelsbest<<-labelsbest
          cluster_id_previous=c()
          n_cluster=nClustersi(labelsbest)
          new_C_ID<<-length(n_cluster)+1
          for (k in 1:length(n_cluster)){
            l=n_cluster[k] #which label has been attributed to it during the run
            cluster_id_previous=c(cluster_id_previous, l)
          }
          new_C_ID<<-length(n_cluster)+1
          cluster_id_previous<<-cluster_id_previous
          cluster_id=c()
          cluster_id=cluster_id_previous
          cluster_id<<-cluster_id

        }else{

        label_previous=c()
        label_previous=labelsbest 
        last_over<<-Starting_frame_id+dirname-2
        count<<-0
        
        for (g in 1:s_prev){
          if (frame_prev[g]==last_over){
            count<<-count+1
           }
          }
        label_previous_cropped=c()
        label_previous_cropped=label_previous[(count+1):s_prev]
        s_previous_cropped=length(label_previous_cropped)
        nb_new=s-s_previous_cropped 
        
        labelsbest=c()
        labelsbest=strsplit(readLines(nbfile),",")[[1]]
        labelsbest<<-labelsbest
        labelsbest_with_new_id=c()
        labelsbest_with_new_id=labelsbest
        labelsbest_with_new_id<<-labelsbest_with_new_id
        
        
        if (coming_from_break==0){

        }else{
          
          #################################################################################
          ###############################General case######################################
        cluster_id=c()
        frame_i_am_in=Starting_frame_id+dirname-1 
        frame_i_am_in<<-frame_i_am_in 
        n_cluster=nClustersi(labelsbest) 
        init_noise=max(cluster_id_previous) 
        print("max id from previous frame")
        print(init_noise)
        
      
        for (p in 1:s){ 
          noiseornot=0
          for (k in 1:length(n_cluster)){
            if (as.numeric(labelsbest[p])==n_cluster[k]){
              noiseornot=1
            }
          }
          if (noiseornot==0){
            labelsbest_with_new_id[p]=as.numeric(labelsbest[p])+init_noise 
          }
        }
        ##
        
        
 
        nb_clusters_previous=length(cluster_id_previous)
        already_contributed=c()
        Keep_track_of_id=c()
        for (k in 1:length(n_cluster)){ 
          l=n_cluster[k] 
          from=c()
          from<<-from
          at_least_from_one_cluster=0
          noise=0
          New=0
          New<<-New
          noise<<-noise
          MERGING=0
          SPLITTING=0
          TOT_NB_POINTS=0
          at_least_from_one_cluster=0
          TOT_NB_POINTS<<-TOT_NB_POINTS
          ##
          for (ll in 1:s){ 
            at_least_from_one_cluster1=0
            if (as.numeric(labelsbest[ll])==l){ 
              TOT_NB_POINTS=TOT_NB_POINTS+1
              TOT_NB_POINTS<<-TOT_NB_POINTS
              if (ll<=s_previous_cropped){
                for (lll in 1:nb_clusters_previous){ 
                  if (label_previous_cropped[ll]==cluster_id_previous[lll]){
                    
                    at_least_from_one_cluster=1
                    at_least_from_one_cluster1=1
              
                    s_from=length(from)/2 
                    already_in=0
                    
                    if (s_from!=0){
                      for (x in 1:s_from){
                        if (from[x,1]==label_previous_cropped[ll]){
                          already_in=1
                          from[x,2]=from[x,2]+1
                          
                          from<<-from
                        }
                        
                      }
                    }
                 
                    if (already_in==0){ # It is a new track involved
                      if(s_from!=0){
                        print("merging")
                        MERGING=1
                        M=c(dirname,cluster_id_previous[lll])
                        merging=rbind(merging,M)
                        merging<<-merging
                        
                      }
                      v=c(cluster_id_previous[lll],1)
                      from=rbind(from, v) 
                    }
                   
                  }
                } 
                
                if (at_least_from_one_cluster1==0){
                  noise=noise+1
                  noise<<-noise
                }
              }else{ # for new points!
                 New=New+1
                 New<<-New
              }
            }
          }
          
          
          
          from<<-from

          
          gh=length(Keep_track_of_id) 
          ghh=length(from)/2
          if (gh!=0 && ghh!=0){
            for (t in 1:gh){
              for (tt in 1:ghh){
                if(Keep_track_of_id[t]==from[tt,1]){ 
                  SPLITTING=1
                  
                  sp=c(dirname,from[tt,1]) 
                  splitting=rbind(splitting, sp)
                  splitting<<-splitting
                  
                }
              }
            }
          }
          
          Keep_track_of_id=c(Keep_track_of_id, from[,1])

           if (at_least_from_one_cluster==0){
            
            new_ID=new_C_ID
            new_C_ID<<-new_C_ID+1
            Coord=c()
            pp=0
            
            for (ll in 1:s){
              if (labelsbest[ll]==l){
                labelsbest_with_new_id[ll]=new_ID
                if (pp==0){
                  Coord=pts[ll,]
                  pp=1
                }else{
                  Coord=rbind(Coord,pts[ll,])
                }
                
              }
            }
            av_x=mean(Coord[,1])
            av_y=mean(Coord[,2])

    
            H=0
            H=molsPerClusteri(labelsbest,l)
            J=0
            J=clusterRadiii(pts, labelsbest, l)
            if (J==-1){
              print("size")
              print(J)
            }
            
            if (first==0){
              s_contribution1=1
              s_contribution2=length(contribution)
              first=1
              inn=0
            }else{
              s_contribution1=dim(contribution)[1]
              s_contribution2=dim(contribution)[2]
              inn=1
            }
            cluster_id=c(cluster_id,new_ID)
            cluster_id<<-cluster_id

            contribution_tempo=c(rep(0, 11))
            contribution_tempo[1]=frame_i_am_in
            contribution_tempo[2]=new_ID
            contribution_tempo[7]=New/TOT_NB_POINTS
            contribution_tempo[8]=noise/TOT_NB_POINTS
            contribution_tempo[5]=H 
            contribution_tempo[6]=J 
            contribution_tempo[3]=av_x
            contribution_tempo[4]=av_y

            
            init_noise=init_noise+1
            init_noise<<-init_noise
            
     
          }else{

            if (MERGING==0){
              if (first==0){
                s_contribution1=1
                s_contribution2=length(contribution)
                first=1
                inn=0
              }else{
                s_contribution1=dim(contribution)[1]
                s_contribution2=dim(contribution)[2]
                inn=1
              }

              if (SPLITTING==0){

                new_ID=from[1,1]
                Coord=c()
                pp=0
                for (ll in 1:s){
                  if (labelsbest[ll]==l){
                    labelsbest_with_new_id[ll]=new_ID
                    if (pp==0){
                      Coord=pts[ll,]
                      pp=1
                    }else{
                      Coord=rbind(Coord,pts[ll,])
                    }
                  }
                }
                av_x=mean(Coord[,1])
                av_y=mean(Coord[,2])
                
      
                H=0
                H=molsPerClusteri(labelsbest,l)
                J=0
                J=clusterRadiii(pts, labelsbest, l)
                if (J==-1){
                  print("size")
                  print(J)
                }

                
                contribution_tempo=c(rep(0, 11))
                contribution_tempo[1]=frame_i_am_in
                contribution_tempo[2]=new_ID
                contribution_tempo[7]=New/TOT_NB_POINTS
                contribution_tempo[8]=noise/TOT_NB_POINTS
                contribution_tempo[5]=H
                contribution_tempo[6]=J
                contribution_tempo[3]=av_x
                contribution_tempo[4]=av_y
                contribution_tempo[9]=from[2]/TOT_NB_POINTS
                cluster_id=c(cluster_id,new_ID)
                cluster_id<<-cluster_id
                
              
           
                
              }else{
                

                new_ID=new_C_ID
                new_C_ID<<-new_C_ID+1
                cluster_id=c(cluster_id,new_ID)
                cluster_id<<-cluster_id
            
                init_noise=init_noise+1
                init_noise<<-init_noise
                Coord=c()
                pp=0
                for (ll in 1:s){
                  if (labelsbest[ll]==l){
                    labelsbest_with_new_id[ll]=new_ID
                    if (pp==0){
                      Coord=pts[ll,]
                      pp=1
                    }else{
                      Coord=rbind(Coord,pts[ll,])
                    }
                  }
                }
                av_x=mean(Coord[,1])
                av_y=mean(Coord[,2])
                
                H=0
                H=molsPerClusteri(labelsbest,l)
                J=0
                J=clusterRadiii(pts, labelsbest, l)
                if (J==-1){
                  print("size")
                  print(J)
                  print(l)
                  print(labelsbest)
                }
                contribution_tempo=c(rep(0, 11))
                contribution_tempo[1]=frame_i_am_in
                contribution_tempo[2]=new_ID
                contribution_tempo[7]=New/TOT_NB_POINTS
                contribution_tempo[8]=noise/TOT_NB_POINTS
                contribution_tempo[5]=H
                contribution_tempo[6]=J
                contribution_tempo[3]=av_x
                contribution_tempo[4]=av_y
                contribution_tempo[10]=from[2]/TOT_NB_POINTS
                

                
                s_contribution1=dim(contribution)[1]
               
                just_last_occurence=1
             
                new_ID=new_C_ID
                new_C_ID<<-new_C_ID+1
             
                cluster_id=c(cluster_id,new_ID)
                cluster_id<<-cluster_id
                init_noise=init_noise+1
                init_noise<<-init_noise
                Coord=c()
                pp=0
                for (ll in 1:s){
                  if (labelsbest_with_new_id[ll]==from[1,1]){
                    labelsbest_with_new_id[ll]=new_ID
                  }
                }
                for (d in 1:(k-1)){
                  if (contribution[(s_contribution1-d),2]==from[1,1]){
                      contribution[(s_contribution1-d),2]=new_ID
                      
                      print(from[1,1])
                    
                      print(new_ID)
                  }
                }
              }
             
             }else{
              

              new_ID=new_C_ID
              new_C_ID<<-new_C_ID+1
              init_noise=init_noise+1
              init_noise<<-init_noise
              Coord=c()
              pp=0
              for (ll in 1:s){
                if (labelsbest[ll]==l){
                  labelsbest_with_new_id[ll]=new_ID
                  if (pp==0){
                    Coord=pts[ll,]
                    pp=1
                  }else{
                    Coord=rbind(Coord,pts[ll,])
                  }
                }
              }
              av_x=mean(Coord[,1])
              av_y=mean(Coord[,2])
              
              
              H=0
              H=molsPerClusteri(labelsbest,l)
              J=0
              J=clusterRadiii(pts, labelsbest, l)
              if (J==-1){
                print("size")
                print(J)
                print(l)
                print(labelsbest)
              }
              
              
              cluster_id=c(cluster_id,new_ID)
              cluster_id<<-cluster_id
              

              contribution_tempo=c(rep(0, 11))
              contribution_tempo[1]=frame_i_am_in
              contribution_tempo[2]=new_ID
              contribution_tempo[7]=New/TOT_NB_POINTS
              contribution_tempo[8]=noise/TOT_NB_POINTS
              contribution_tempo[5]=H
              contribution_tempo[6]=J
              contribution_tempo[3]=av_x
              contribution_tempo[4]=av_y
              for (x in 1:dim(from)[1]){
                contribution_tempo[11]=contribution_tempo[11]+from[x,2]/TOT_NB_POINTS
              }
            }
          }
          contribution=rbind(contribution, contribution_tempo)
          contribution<<-contribution
        
        }

        labelsbest=c()
        labelsbest=labelsbest_with_new_id
        labelsbest<<-labelsbest
        contribution<<-contribution
        
        
        cluster_id_previous=c()
        cluster_id_previous=cluster_id #first rota you ll have to increment it!
        cluster_id_previous<<-cluster_id_previous

        }
        }

        

        wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
        cat("The best: clusterscale", bestcs, " thresh", bestthr, "labels.txt\nNumber of clusters:", nClusters(labelsbest), "\nPercentage in clusters: ", percentageInCluster(labelsbest), "%\nMean number of molecules per cluster: ", nMolsPerCluster(labelsbest), "\nMean radius: ", mean(clusterRadii(pts, labelsbest)), sep="", file=wfile)
   
        if (makeplot){
            if (!superplot){
                    tiff(file.path(paste(nfoldername, "/plot.tiff", sep="")))
                    axes=TRUE
                    cex=1
                }           
            if (dim(data)[2]==5){
                labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})

                par(pty="s")
                par(mfrow=c(1,2))
                par(mar=c(4,4,.5, .5))
                par(oma=c(1,1,1,1))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
                                        #X11()
                par(mar=c(4,4,.5, .5))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcolsi(labelsbest,colors,cluster_id), sub="Estimated",xlab="",ylab="")
                if (!superplot) dev.off()
            }
            else {
                par(pty="s")
                par(mar=c(0,0,0,0))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcolsi(labelsbest,colors,cluster_id), sub="Clustering",xlab="",ylab="", pch=16, cex=cex, axes=axes)
                box()
                if (!superplot) dev.off()
            }
        }
        
        list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), pclustered=percentageInCluster(labelsbest), totalmols=length(labelsbest), reldensity=reldensity(pts, labelsbest, xlim, ylim))

        #Not sure what that is for
        hit=0
        hit<<-hit
        frame_prev=c()
        frame_prev=frame
        frame_prev<<-frame_prev
        s_prev=s
        s_prev<<-s_prev
        
        } 
        

        else{if (beenthereonce==0){
          Starting_frame_id=frame[1]
          Starting_frame_id<<-Starting_frame_id
          print("in before analysis start")
          labelsbest=c()
          #all in the noise
          for (p in 1:s){
            labelsbest=c(labelsbest, p+50) # Just to be sure not to be attributed tp a cluster when the analtysis does start
          }
          labelsbest<<-labelsbest
          ##Some summaries
          wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
          cat(("No analysis for lack of points\nNumber of clusters: NA\nPercentage in clusters: NA \nMean number of molecules per cluster: NA\nMean radius: NA"), sep="", file=wfile)
          ##browser()
          if (makeplot){
            if (!superplot){
              tiff(file.path(paste(nfoldername, "/plot.tiff", sep="")))
              axes=TRUE
              cex=1
            }           
            if (dim(data)[2]==5){
              labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})
              
              par(pty="s")
              par(mfrow=c(1,2))
              par(mar=c(4,4,.5, .5))
              par(oma=c(1,1,1,1))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
              #X11()
              par(mar=c(4,4,.5, .5))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Estimated",xlab="",ylab="")
              if (!superplot) dev.off()
            }
            else {
              par(pty="s")
              par(mar=c(0,0,0,0))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcolsi(labelsbest,colors,cluster_id), sub="Clustering",xlab="",ylab="", pch=16, cex=cex, axes=axes)
              box()
              if (!superplot) dev.off()
            }
          }
          frame_prev=c()
          frame_prev=frame
          frame_prev<<-frame_prev
          s_prev=s
          s_prev<<-s_prev
          
        }else{ 
          print("in break")
          #coming_from_break<<-0
          if (hit==0){last_over<<-Starting_frame_id+dirname-2
     
          }

          count<<-0
          print("prev frame id")
          q=last_over+hit
          print(q)
          
          for (g in 1:s_prev){
            if (frame_prev[g]==last_over+hit){
              #counting the old points in there
              count<<-count+1
            
            }}
          print("count")
          print(count)

          hit=hit+1
          hit<<-hit
          label_croped=c()
          S_prev_label=length(labelsbest)
          MAX=0
          MAX<<-MAX

          for (d in (count+1):S_prev_label){
            label_croped=c(label_croped, labelsbest[d])
            #You copy all the old ones you are just gonna add the one from the new frame (if any)
             if (as.numeric(labelsbest[d])-MAX>=0){
               #keep the highest if id for the noise
               MAX=as.numeric(labelsbest[d])
               MAX<<-MAX
               
             }
          }

          label_max=MAX+1
          label_max<<-label_max
          S_label=length(label_croped)
          
          print("points to add")
          q=s-S_label
          print(q)

          if (s-S_label>0){
            l=s-S_label
            for (p in 1:l) {

              label_croped=c(label_croped, label_max)
              label_max=label_max+1
              label_max<<-label_max
              
            }
            
          }
          #now should be the size of data ie s


          labelsbest=c()
          labelsbest=label_croped
          labelsbest<<-labelsbest
          #print("I am out")
          ##Some summaries
          wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
          cat(("No analysis for lack of points\nNumber of clusters: NA\nPercentage in clusters: NA \nMean number of molecules per cluster: NA\nMean radius: NA"), sep="", file=wfile)
##browser()
          if (makeplot){
            if (!superplot){
              tiff(file.path(paste(nfoldername, "/plot.tiff", sep="")))
              axes=TRUE
              cex=1
            }           
            if (dim(data)[2]==5){
              labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})
              
              par(pty="s")
              par(mfrow=c(1,2))
              par(mar=c(4,4,.5, .5))
              par(oma=c(1,1,1,1))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
              #X11()
              par(mar=c(4,4,.5, .5))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Estimated",xlab="",ylab="")
              if (!superplot) dev.off()
            }
            else {
              par(pty="s")
              par(mar=c(0,0,0,0))
              plot(pts, xlim=xlim, ylim=ylim, col=mkcolsi(labelsbest,colors,cluster_id), sub="Clustering",xlab="",ylab="", pch=16, cex=cex, axes=axes)
              box()
              if (!superplot) dev.off()
            }
          }
          ##browser()
          list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), pclustered=percentageInCluster(labelsbest), totalmols=length(labelsbest), reldensity=reldensity(pts, labelsbest, xlim, ylim))
          
        frame_prev=c()
        frame_prev=frame
        frame_prev<<-frame_prev
        s_prev=s
        s_prev<<-s_prev
        }
        }
        iop=length(splitting)
        if (splitting[iop]!=0){
          splitting[iop]=splitting[iop]+1
        }
        
}
        #})
    
    
    ###############################################################################################################
    #print("I am here")
    merging_to_print=c()
    merging_to_print=merging
    splitting_to_print=c()
    splitting_to_print=splitting
    write.csv(contribution,file = paste(nexpname, "/Contribution.csv", sep="/"), sep=",")
    write.csv(merging_to_print,file = paste(nexpname, "/Merging_over_time.csv", sep="/"), sep=",")
    write.csv(splitting_to_print,file = paste(nexpname, "/Splitting_over_time.csv", sep="/"), sep=",")
   

})
