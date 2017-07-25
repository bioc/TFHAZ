# Finds transcription factor dense DNA zones for different accumulation threshold values.
# input:
#       accumulation: a list of four elements, as the output of the accumulation function.
#       threshold_step: an integer, the step used to calculate the threshold values.
# output: a list of seven elements, and, for each threshold value, a ".bed" file with the chromosome and genomic coordinates of the dense zones found:
#       zones_count: a data frame containing the considered threshold values and the number of dense zones obtained with each of the threshold values.
#       bases_count: a data frame containing the considered threshold values and the total number of bases belonging to the dense zones obtained with each of the threshold values.
#       lengths: a data frame containing the considered threshold values and min, max, mean, median and standard deviation of the dense zone lengths obtained with each of the considered threshold values.
#       distances: a dataframe containing the considered threshold values and  min, max, mean, median and standard deviation of the distances between adjacent dense zones obtained with each of the threshold values.
#       acctype: a string with the accumulation type used.
#       chr: a string with the chromosome name associated with the accumulation vector used.
#       w: an integer with half-width of the window used to calculate the accumulation vector.




dense_zones=function(accumulation,threshold_step){


  #initializing vectors
  n_zones=vector()
  n_bases=vector()
  length_zone_min=vector()
  length_zone_max=vector()
  length_zone_mean=vector()
  length_zone_median=vector()
  length_zone_sd=vector()
  dist_zone_min=vector()
  dist_zone_max=vector()
  dist_zone_mean=vector()
  dist_zone_median=vector()
  dist_zone_sd=vector()
  acc_tf=accumulation$accvector
  threshold=seq(1,max(acc_tf),threshold_step)


  for(i in 1:length(threshold)) {

    #d: 0-1 sparse vector. d=1 if the number of accumulated TFs is equal or larger than threshold
    d=sparseVector(0,0,length(acc_tf)) #initialize d: zeros sparse vector with length equal to the length of the accumulation vector
    d[acc_tf>=threshold[i]]=1          #set d=1 where acc_tf>=threshold



    #counting the number of dense zones that are formed by contiguos bases given a threshold.
    n_bases[i]=sum(d==1)
    n_zones[i]=sum((sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1))==-1)
    #finding dense zones starting points
    start_zone=which(sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1)==-1)
    #finding dense zones ending points
    end_zone=which(sparseVector(d@x,d@i+1,d@length+1)-sparseVector(d@x,d@i,d@length+1)==1)-1

    # writing on files chromosome and positions of starting and ending points
     chr=rep(accumulation$chr,length(start_zone))
     f=data.frame(chr,start_zone,end_zone)
     write.table(f,file = paste(accumulation$acctype,"_acc_w_",accumulation$w,"_",accumulation$chr,"_dense_zones_th_",threshold[i],".bed",sep=""),row.names = FALSE, col.names = FALSE, quote = FALSE)
    # finding elements of length dataframe
    length_zone=end_zone-start_zone+1
    length_zone_min[i]=min(length_zone)
    length_zone_max[i]=max(length_zone)
    length_zone_mean[i]=mean(length_zone)
    length_zone_median[i]=median(length_zone)
    length_zone_sd[i]=sd(length_zone)


    # calculating distances of dense zones
    if (length(start_zone)>1) {

      dist_zone=start_zone[-1]-end_zone[-length(end_zone)]
      dist_zone_min[i]=min(dist_zone)
      dist_zone_max[i]=max(dist_zone)
      dist_zone_mean[i]=mean(dist_zone)
      dist_zone_median[i]=median(dist_zone)
      dist_zone_sd[i]=sd(dist_zone)
    }
    else {
      dist_zone=NA
      dist_zone_min[i]=NA
      dist_zone_max[i]=NA
      dist_zone_mean[i]=NA
      dist_zone_median[i]=NA
      dist_zone_sd[i]=NA


    }
  }

  # result dataframes creation

  zones_count=data.frame(threshold,n_zones)
  bases_count=data.frame(threshold,n_bases)
  lengths=data.frame(threshold,n_zones,length_zone_min,length_zone_max,length_zone_mean,length_zone_median,length_zone_sd)
  distances=data.frame(threshold,n_zones,dist_zone_min,dist_zone_max,dist_zone_mean,dist_zone_median,dist_zone_sd)
  return(list("zones_count"=zones_count,"bases_count"=bases_count,"lengths"=lengths,"distances"=distances,"chr"=accumulation$chr,"w"=accumulation$w,"acctype"=accumulation$acctype))

}
