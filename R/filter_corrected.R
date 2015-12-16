library(reshape2)
library(plyr)
library(ggplot2)

mass.error <- 0.005
max.rt.drift_in <- 0.05

m1 <- list(min=(1.0034 - mass.error), max=(1.0034 + mass.error))
m2 <- list(min=(2.0043 - mass.error), max=(2.0043 + mass.error))
m3 <- list(min=(3.0077 - mass.error), max=(3.0077 + mass.error))

frag1 <- list(min=(60.0212 - mass.error), max=(60.0212 + mass.error))

to.elim.2 <- list(m1=list(m1,"+"), 
                m2=list(m2,"+"), 
                m3=list(m3,"+")) 

to.elim.1 <- list(frag1=list(frag1, "-"))

testCount1 <- 0
testCount2 <- 0

filter <- function(peakFrameIn, massShiftList1=to.elim.1, massShiftList2=to.elim.2, max.rt.drift=max.rt.drift_in){

    # peakFrameIn: this is the peak information (like mass and RT). Peaks are in the rows.
    # massShiftList1: this is a list of mass shift ranges to eliminate. Each entry in this list
    #       is also a list that has "min" and "max" entries for the range. Peaks that are found as shifts from this list
    #       are removed, but shifts from massShiftList2 are still calculated off of these removed peaks. This is to enable
    #       fragment peaks (subject to removal) to also have isotope peaks (that are also removed).
    # massShiftList2: See massShiftList1.
    # max.rt.drift: this is an RT window to select peaks having the "same" retention time.

    ref.dir <- vector(mode="character", length=dim(peakFrameIn)[[1]])
    ref.type <- vector(mode="character", length=dim(peakFrameIn)[[1]])
    ref.parent <- vector(mode="character", length=dim(peakFrameIn)[[1]])

    #####
    for(i in 1:(dim(peakFrameIn)[[1]]-1)) {

        peak1 <- peakFrameIn[i,]
        peak1.mz <- peak1$mz
        peak1.rt <- peak1$rt

        #print(paste("Working on peak",i,". Counts are",testCount1,"and",testCount2))
        print(paste0("Working on peak ",i,". Num discarded: ", testCount1, " -- ", testCount2 ))
        
        if(!(ref.type[[i]] %in% names(to.elim.1)) && !(ref.type[[i]] %in% names(to.elim.2)) ) {

            for(j in (i+1):(dim(peakFrameIn)[[1]])) {

                if( !(ref.type[[j]] %in% names(to.elim.1)) && !(ref.type[[j]] %in% names(to.elim.2)) ) {
                    # +++
                    # stuff
                    peak2 <- peakFrameIn[j,]
                    peak2.mz <- peak2$mz
                    peak2.rt <- peak2$rt


                    if(abs(peak1$rt - peak2$rt) <= max.rt.drift) {

                        # This peak is a candidate for removal based on RT.
                        # Now check its mass differences.

                        for(elim in names(massShiftList1)){

                            dm.low <- massShiftList1[[elim]][[1]][["min"]]
                            dm.high <- massShiftList1[[elim]][[1]][["max"]]

                            if(massShiftList1[[elim]][[2]] == "+"){
                                # Peak is mass shifted higher
                                if((peak2.mz >= (peak1.mz + dm.low)) && (peak2.mz <= (peak1.mz + dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "higher"
                                    ref.parent[[j]] <- as.character(i)
                                    # Why are you freezing all the time?
                                }
                            }

                            if(massShiftList1[[elim]][[2]] == "-"){
                                # Peak is mass shifted lower
                                if((peak2.mz <= (peak1.mz - dm.low)) && (peak2.mz >= (peak1.mz - dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "lower"
                                    ref.parent[[j]] <- as.character(i)
                                }
                            }

                        }

                        for(elim in names(massShiftList2)){

                            dm.low <- massShiftList2[[elim]][[1]][["min"]]
                            dm.high <- massShiftList2[[elim]][[1]][["max"]]

                            if(massShiftList2[[elim]][[2]] == "+"){
                                # Peak is mass shifted higher
                                if((peak2.mz >= (peak1.mz + dm.low)) && (peak2.mz <= (peak1.mz + dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "higher"
                                    ref.parent[[j]] <- as.character(i)
                                    # Why are you freezing all the time?
                                }
                            }

                            if(massShiftList2[[elim]][[2]] == "-"){
                                # Peak is mass shifted lower
                                if((peak2.mz <= (peak1.mz - dm.low)) && (peak2.mz >= (peak1.mz - dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "lower"
                                    ref.parent[[j]] <- as.character(i)
                                }
                            }
   
                        }
                    }
                    #---
                }
            }
        }
        if((ref.type[[i]] %in% names(to.elim.1)) && !(ref.type[[i]] %in% names(to.elim.2)) ) {
            for(j in (i+1):(dim(peakFrameIn)[[1]])){
                if( !(ref.type[[j]] %in% names(to.elim.1)) && !(ref.type[[j]] %in% names(to.elim.2)) ){
                    # +++
                    # stuff
                    peak2 <- peakFrameIn[j,]
                    peak2.mz <- peak2$mz
                    peak2.rt <- peak2$rt


                    if(abs(peak1$rt - peak2$rt) <= max.rt.drift) {

                        # This peak is a candidate for removal based on RT.
                        # Now check its mass differences.

                        for(elim in names(massShiftList2)){

                            dm.low <- massShiftList2[[elim]][[1]][["min"]]
                            dm.high <- massShiftList2[[elim]][[1]][["max"]]

                            if(massShiftList2[[elim]][[2]] == "+"){
                                # Peak is mass shifted higher
                                if((peak2.mz >= (peak1.mz + dm.low)) && (peak2.mz <= (peak1.mz + dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "higher"
                                    ref.parent[[j]] <- as.character(i)
                                    # Why are you freezing all the time?
                                }
                            }

                            if(massShiftList2[[elim]][[2]] == "-"){
                                # Peak is mass shifted lower
                                if((peak2.mz <= (peak1.mz - dm.low)) && (peak2.mz >= (peak1.mz - dm.high)) ){
                                    testCount1 <<- testCount1 + 1
                                    ref.type[[j]] <- elim
                                    ref.dir[[j]] <- "lower"
                                    ref.parent[[j]] <- as.character(i)
                                }
                            }
   
                        }
                    }
                    #---

                }
            }
        }
    }
    return(list(ref.dir=ref.dir,
                ref.type=ref.type,
                ref.parent=ref.parent,
                keep=which(ref.type==""),
                discard=which(ref.type!="")))
}

test_posMode <- function(){
    d <- read.csv("testData.csv", header=T, stringsAsFactors=F)
    d$mz <- as.numeric(d$mz)
    d$rt <- as.numeric(d$rt)
    print("Starting filtering...")
    t1 <- Sys.time()
    filter(d)
    t2 <- Sys.time()
    print(paste0("Finished. Time taken was ", t2-t1))
    return(list(t1=t1, t2=t2))
}

test_negMode <- function(){
    d <- read.csv("150319_NegIon-AllSamples-hhpoolref_1.csv", header=T, stringsAsFactors=F)
    d2 <- d[,1:3]
    names(d2) <- c("id","mz","rt")
    d2$mz <- as.numeric(d2$mz)
    d2$rt <- as.numeric(d2$rt)
    d3 <- d2[1:1000,]
    wtf <<- d3
    print("Starting filtering...")
    ret <- filter(d3)
    print("Finished filtering.")
    wtf2 <<- ret
    out <- data.frame(d3, ret$ref.parent, ret$ref.type, ret$ref.dir)

    # write the output to csv for Rose.
    out.rose <- out[, 1:5]
    names(out.rose) <- c("id", "mz", "rt", "parent", "type")
    write.csv(out.rose, file="negIonFilterTest_top1000.csv", row.names=F, quote=F)

    return(list(all=out, rem.vec=ret$rem.vec, orig=d3, ret=ret))
}

run_full_negMode <- function(){
    d <- read.csv("150319_NegIon-AllSamples-hhpoolref_1.csv", header=T, stringsAsFactors=F)
    d2 <- d[,1:3]
    names(d2) <- c("id","mz","rt")
    d2$mz <- as.numeric(d2$mz)
    d2$rt <- as.numeric(d2$rt)
    print("Starting filtering...")
    t1 <- Sys.time()
    ret <- filter(d2)
    t2 <- Sys.time()
    print("Finished filtering.")
    print(paste0("Time taken was ", t2-t1))
    out <- data.frame(d2, ret$ref.parent, ret$ref.type, ret$ref.dir)

    return(list(t1=t1, t2=t2, all=out))
}

testLoop <- function(){

    test <- c(1,2,3,4,5,6,7,8,9)

    for(i in 1:(length(test)-1)) {
        for(j in (i+1):length(test)) {
            print(paste(i,"and",j))
        }
    }
}
