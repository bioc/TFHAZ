# This function finds transcription factor high accumulation DNA zones (TFHAZ).
# Each TFHAZ is formed by contiguous bases with accumulation higher than the threshold (TH).
# input:
#       accumulation: list of four elements containing: a sparse vector with accumulation values (e.g.,obtained with the accumulation function), the accumulation type, a chromosome name, and the half-width of the window used for the accumulation count.
# output: A list of eight elements:
#        n_zones: an integer containing the number of dense zones obtained.
#        n_bases: an integer containing the total number of bases belonging to the dense zones obtained.
#        lengths: a data frame containing the considered threshold value and min, max, mean, median and standard deviation of the dense zone lengths obtained.
#        distances: a dataframe containing the considered threshold value and  min, max, mean, median and standard deviation of the distances between adjacent dense zones obtained with each of the threshold values.
#        TH: a number with the threshold value found.
#        acctype: a string with the accumulation type used.
#        chr: a string with the chromosome name associated with the accumulation vector used.
#        w: an integer with half-width of the window used to calculate the accumulation vector.
# A ".bed" file with the chromosome and genomic coordinates of the dense zones found is created.
# A ".png" file with the plot of the TFHAZ found along the cromosome is created.

high_accumulation_zones= function(accumulation) {

    # input accumulation vector
    acc_tf=accumulation$accvector

    # finding the threshold
    TH=mean(acc_tf[acc_tf>0])+2*sd(acc_tf[acc_tf>0])

    #d: 0-1 sparse vector. d=1 if the number of accumulated TFs is equal or larger than threshold
    d=sparseVector(0,0,length(acc_tf)) #initialize d: zeros sparse vector with length equal to the length of the accumulation vector
    d[acc_tf>TH]=1          #set d=1 where acc_tf>=threshold
    n_bases=sum(d==1)
    n_zones=sum((sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1))==-1)
    #finding dense zones starting points
    start_zone=which(sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1)==-1)
    #finding dense zones ending points
    end_zone=which(sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1)==1)-1

    # plot
    png(filename=paste("high_accumulation_zones_TH_",round(TH,digits = 1),"_",accumulation$acctype,"_acc_w_",accumulation$w,"_",accumulation$chr,".png",sep=""))
    plot(accumulation$accvector,type="l",xlab="base",ylab = paste("# of ",accumulation$acctype,"s",sep=""))
    points(which(d==1),rep(0,length(which(d==1))),col="red",pch=15)
    abline(h=TH,col="red")
    legend("topleft",legend=c("threshold","high accumulation zones"),col="red",pch = c(NA,15),lty=c(1,NA))
    dev.off()

    # writing on files chromosome and positions of starting and ending points
    chr=rep(accumulation$chr,length(start_zone))
    f=data.frame(chr,start_zone,end_zone)
    write.table(f,file = paste("high_accumulation_zones_TH_",round(TH,digits = 1),"_",accumulation$acctype,"_acc_w_",accumulation$w,"_",accumulation$chr,".bed",sep=""),row.names = FALSE, col.names = FALSE, quote = FALSE)
    # finding elements of length dataframe
    length_zone=end_zone-start_zone+1
    length_zone_min=min(length_zone)
    length_zone_max=max(length_zone)
    length_zone_mean=mean(length_zone)
    length_zone_median=median(length_zone)
    length_zone_sd=sd(length_zone)


    # calculating distances of dense zones
    if (length(start_zone)>1) {

        dist_zone=start_zone[-1]-end_zone[-length(end_zone)]
        dist_zone_min=min(dist_zone)
        dist_zone_max=max(dist_zone)
        dist_zone_mean=mean(dist_zone)
        dist_zone_median=median(dist_zone)
        dist_zone_sd=sd(dist_zone)
    }
    else {
        dist_zone=NA
        dist_zone_min=NA
        dist_zone_max=NA
        dist_zone_mean=NA
        dist_zone_median=NA
        dist_zone_sd=NA
    }

    # datasets creation
    lengths=data.frame(TH,n_zones,length_zone_min,length_zone_max,length_zone_mean,length_zone_median,length_zone_sd)
    distances=data.frame(TH,n_zones,dist_zone_min,dist_zone_max,dist_zone_mean,dist_zone_median,dist_zone_sd)
    return(list("n_zones"=n_zones,"n_bases"=n_bases,"lengths"=lengths,"distances"=distances,"TH"=TH,"chr"=accumulation$chr,"w"=accumulation$w,"acctype"=accumulation$acctype))

}
