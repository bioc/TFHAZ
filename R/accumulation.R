# Creates a vector with accumulation counts of transcription factors for each chromosome base.
# input:
#        data: a GRanges object containing coordinates of TF binding regions and their TF name.
#        acctype: a string with the name of the accumulation type: "TF", "region", "base".
#        chr: a string with the name of the chromosome (e.g., "chr1").
#        w: an integer, half-width of the window that defines the neighborhood of each base.
# output: A list of four elements.
#        accvector: a sparse vector containing the accumulation for each base of the selected chromosome.
#        acctype: a string with the accumulation type used.
#        chr: a string with the chromosome name associated with the accumulation vector.
#        w: an integer with the half-width of the window used to calculate the accumulation vector.
#
    accumulation=function(data,acctype,chr,w){

    if (acctype=="TF") {
    # TF accumulation

    #dataset with only regions of the selected chromosome
    chrom=data[seqnames(data)==chr]

    #find the number of different TF present in the chromosome
    tf=unique(elementMetadata(data))[,1]

    # acc: for each base, it is the sum of the number of TF accumulated in the neighborhood of the considered base
    # the computation of the TF accumulation vector is performed by considering the contribution of all the single regions containing TF peaks on each base
    
    #sparse vector initializing
    
    acc=sparseVector(0,0,max(end(chrom))+w)
                        for (k in 1:length(tf)) {
      
                        #dataframe for a single TF and the selected chromosome 
                        chrom_tf=data[seqnames(data)==chr & elementMetadata(data)[,1]==tf[k]]
      
      #v_tf: accumulation sparse vector for a certain TF. Each postion of the vector corresponds to a base of the genome,while the value corresponds to the number of different TFs present in that base. In order to sum all the TF sparse vectors (v_tf) all of them are initialized with the same length, corrisponding to the maximum ending position of the input regions plus the half window width selected
      
      #v_tf initializing
      v_tf=sparseVector(0,0,max(end(chrom))+w)
      
      #for each input region, the value related to the bases covered by the window (region start position-w : region end position+w) is set to 1.
      for (i in 1:length((start(chrom_tf)))) {
                if (start(chrom_tf[i])-w <=0) {
                    v_tf[1:(end(chrom_tf[i])+w)]=v_tf[1:(end(chrom_tf[i])+w)]+1
                } 
                else {
                    
                    v_tf[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]=v_tf[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]+1
                }
                
            }
            #For each TF the max possible acccumlation value is 1 
            v_tf[v_tf>=1]=1
            acc=acc+v_tf
      
    }
    return(list("accvector"=acc,"chr"=chr,"w"=w,"acctype"=acctype))
    }

    else if (acctype=="region") {
    # region accumulation           
    #dataset with only regions of the selected chromosome
    chrom=data[seqnames(data)==chr]

    #find the number of different TF present in the chromosome
    tf=unique(elementMetadata(data))[,1]

    # acc: for each base, it is the sum of the number of regions containing TF in the neighborhood of the considered base
    # the computation of the region accumulation vector is performed by considering the contribution of all the single regions containing TF peaks on each base
    
    #sparse vector initializing
    acc=sparseVector(0,0,max(end(chrom))+w)
                        for (k in 1:length(tf)) {
      
                        #dataframe for a single TF and the selected chromosome 
                        chrom_tf=data[seqnames(data)==chr & elementMetadata(data)[,1]==tf[k]]
      
      #v_tf: accumulation sparse vector for a certain TF. Each postion of the vector corresponds to a base of the genome,while the value corresponds to the number of different TFs present in that base. In order to sum all the TF sparse vectors (v_tf) all of them are initialized with the same length, corrisponding to the maximum ending position of the input regions plus the half window with selected
      
      #v_tf initializing
      v_tf=sparseVector(0,0,max(end(chrom))+w)
      
      #for each input region, the value related to the bases covered by the window (region start position-w+region end position+w) is set to 1.
      for (i in 1:length(start(chrom_tf))) {
                if (start(chrom_tf[i])-w <=0) {
                    v_tf[1:(end(chrom_tf[i])+w)]=v_tf[1:(end(chrom_tf[i])+w)]+1
                } 
                else {
                    
                    v_tf[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]=v_tf[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]+1
                }
                
            }
           
            acc=acc+v_tf
      
    }
    return(list("accvector"=acc,"chr"=chr,"w"=w,"acctype"=acctype))
    }


    else if (acctype=="base"){


    #dataset with only regions of the selected chromosome
    chrom=data[seqnames(data)==chr]

    #find the number of different TF present in the chromosome
    tf=unique(elementMetadata(data))[,1]

    # acc: for each base, it is the sum of the number of bases containing TF in the neighborhood of the considered base
    # the computation of the base accumulation vector is performed by considering the contribution of all the single regions containing TF peaks on each base.
    # This contribution depends on the number of bases and the position of the considered region
    
    #sparse vector initializing
    acc=sparseVector(0,0,max(end(chrom))+w)

                for (k in 1:length(tf)) {
     chrom_tf=data[seqnames(data)==chr & elementMetadata(data)[,1]==tf[k]]
       #v_b: accumulation sparse vector for a certain TF. Each postion of the vector corresponds to a base of the genome, while the value corresponds to the number of bases containing TFs in the neighborhood of the considered base. In order to sum all the TF sparse vectors (v_b) all of them are initialized with the same length, corrisponding to the maximum ending position of the input regions plus the half window with selected
     
      #v_tf initializing
      v_b=sparseVector(0,0,max(end(chrom))+w)
    
    for (i in 1:length(start(chrom_tf))) {
      d=width(chrom_tf[i])
      
      if (start(chrom_tf[i])-w<1) {
        if (d==1) {
          
          v_b[1:(end(chrom_tf[i])+w)]=v_b[1:(end(chrom_tf[i])+w)]+d 
          
        }  else if (d>1 & d<(2*w+1)) {
          if (end(chrom_tf[i])-w<1) {
            v_b[1:(start(chrom_tf[i])+w)]=v_b[1:(start(chrom_tf[i])+w)]+d
            v_b=v_b+sparseVector(-c((start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w))+start(chrom_tf[i])+w+d,(start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w),length(v_b))
          }     else {
            v_b[(end(chrom_tf[i])-w):(start(chrom_tf[i])+w)]=v_b[(end(chrom_tf[i])-w):(start(chrom_tf[i])+w)]+d
            v_b=v_b+sparseVector(c(1:(end(chrom_tf[i])-w-1))-start(chrom_tf[i])+w+1,1:(end(chrom_tf[i])-w-1),length(v_b))  
            v_b=v_b+sparseVector(-c((start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w))+end(chrom_tf[i])+w+1,(start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w),length(v_b))
          } 
        } else {
          v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]=v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]+(2*w+1)
          v_b=v_b+sparseVector((c(1:(start(chrom_tf[i])+w-1))-start(chrom_tf[i])+w+1),1:(start(chrom_tf[i])+w-1),length(v_b))
          v_b=v_b+sparseVector(c(-((end(chrom_tf[i])-w+1):(end(chrom_tf[i])+w))+w+end(chrom_tf[i])+1),(end(chrom_tf[i])-w+1):(end(chrom_tf[i])+w),length(v_b))
          
        }
        
      } else {
        if(d==1){
          v_b[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]=v_b[(start(chrom_tf[i])-w):(end(chrom_tf[i])+w)]+d
          
        }        else if (d>1 & d<(2*w+1)) {
          
          v_b[(end(chrom_tf[i])-w):(start(chrom_tf[i])+w)]=v_b[(end(chrom_tf[i])-w):(start(chrom_tf[i])+w)]+d
          v_b=v_b+sparseVector(c((start(chrom_tf[i])-w):(end(chrom_tf[i])-w-1))-start(chrom_tf[i])+w+1,(start(chrom_tf[i])-w):(end(chrom_tf[i])-w-1),length(v_b))  
          v_b=v_b+sparseVector(-c((start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w))+end(chrom_tf[i])+w+1,(start(chrom_tf[i])+w+1):(end(chrom_tf[i])+w),length(v_b))
          
        }           else {
                          if (w==0) {
                          v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]=v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]+(2*w+1)
                                      }
                    else {
                        v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]=v_b[(start(chrom_tf[i])+w):(end(chrom_tf[i])-w)]+(2*w+1)
                        v_b=v_b+sparseVector((c((start(chrom_tf[i])-w):(start(chrom_tf[i])+w-1))-start(chrom_tf[i])+w+1),(start(chrom_tf[i])-w):(start(chrom_tf[i])+w-1),length(v_b))
                        v_b=v_b+sparseVector(c(-((end(chrom_tf[i])-w+1):(end(chrom_tf[i])+w))+w+end(chrom_tf[i])+1),(end(chrom_tf[i])-w+1):(end(chrom_tf[i])+w),length(v_b))
                           }
        }
      }
    }
      acc=acc+v_b
    }

    return(list("accvector"=acc,"chr"=chr,"w"=w,"acctype"=acctype))
    }

    else {
    print("Wrong accumulation type in input. Accumulation type input available: TF, region, base")
    }
    }
