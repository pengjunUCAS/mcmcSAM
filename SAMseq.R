
 
    ###
    library(rstan)
    ###
    ###------------------------------------------------------------------------------------------------------------
    ### Function mcmcSAMseq() is used for obtaining the posterior distributions of ages of a sequence of 
    ### samples S1, S2, S3, ..., Sn according to satistical age models (SAM) with ordering constraints 
    ### Age(S1) < Age(S2) < Age(S3) < ... < Age(Sn).
    ###
    ### SAM that can be implemented included: (1) Central Age Model, 
    ### (2) Minimum Age Model(3-parameter), (3) Maximum Age Model (3-parameter).
    ### 
    ### Arguments in function mcmcSAMseq():
    ### [EDdata]   a two-column matrix containing equivalent dose (ED) and associated errors for different samples.
    ###            For example, if the numbers of aliquots (grains) for samples S1, S2 and S3 are 30, 40, and 50  
    ###            respectively, and Age(S1)<Age(S2)<Age(S3), then ED data for samples S1, S2, S3 should be 
    ###            presented in the 1-30, 31-70, and 71-120 rows of [EDdata]. In this case, [EDdata] should be 
    ###            a 120-row 2-column matrix. i.e., EDdata=rbind(ED_data_of_S1, ED_data_of_S2, ED_data_of_S3).  
    ###       
    ### [DRdata]   a two-column matrix containing dose rate (DR) and associated errors for different samples.
    ###            For example, if dose rates and associated errors for samples S1, S2, and S3 are 3.1¡À0.23,
    ###            3.3¡À0.27, and 2.9¡À0.19 Gy/ka, respectively,and Age(S1)<Age(S2)<Age(S3), then DR data of 
    ###            sample S1, S2, and S3 should be presented in the 1, 2, and 3 row of [DRdata].
    ###            In this case, [DRdata] should be a 3-row 2-column matrix.
    ###            i.e., DRdata=rbind(c(3.1, 0.23), c(3.3, 0.27), c(2.9,0.19)).
    ###
    ### [iN]       a vector containing number of ED values for each sample in argument [EDdata]. 
    ###            For example, if the numbers of alquots (grains) for samples S1, S2, and S3 are 30, 40, 
    ###            and 50 respectively, and Age(S1)<Age(S2)<Age(S3), then [iN] should be a 3-element vector.    
    ###            i.e., iN=c(30, 40, 50). NOTE THAT the sum of elements in argument [iN] should be equal 
    ###            to the row number of argument [EDdata].
    ###
    ### [priorAge] a two-element vector containing prior ranges (i.e., the lower and upper bounds of a uniform distribution)
    ###            of ages of different samples. For example, if the priors of ages of samples S1, S2, and S3 are between 
    ###            1 and 100 ka, then priorAge=c(1, 100).
    ###
    ### [Sdata]    a vector specifying the additional uncertainties that want be added in quadrature to the relative 
    ###            standard errors of ED values of different samples. For example, if the extra uncertainty for samples 
    ###            S1, S2, and S3 are set equal to 0.1, 0.15, and 0.2, respectively, and Age(S1)<Age(S2)<Age(S3), 
    ###            then Sdata=c(0.1, 0.15, 0.2).
    ###
    ### [model]    a vector containing character strings indicating the fitting model for each sample.
    ###            For example, if the Minimum, Central, and Maximum age models are used to fit De data for
    ###            samples  S1, S2, and S3, then [model] should be a vector containing three character strings.
    ###            i.e., model=c("minimum", "central", "maximum").
    ###
    ### [sDRc]     a real number specifying the systematic error related to dose rate measurement.
    ###
    ### [iflog]    a logical value indicating if the model should be fitted using log-scaled EDdata or not.
    ###
    ### [ordered]  a logical value indicating if the ages will be constrained to be ordered or not.
    ### 
    ### [iter]     a positive integer specifying the number of iterations for each chain.
    ###
    ### [warmup]   a positive integer specifying the number of warmup (i.e.,burn-in) iterations for each chain.
    ###            For example, if warmup=1e4 then samples simulated from the first 1e4 iterations for each 
    ###            chain will be disgarded. 
    ###
    ### [thin]     a positive integer specifying the period for saving samples. For example, if thin=5
    ###            then only samples 1 out of 5 iterations will be reserved from the remaining samples.
    ###
    ### [chains]   a positive integer specifying the number of Markov chains. The number of chain is 
    ###            allowed to be chosen from 1 to 7.
    ###
    ### [cores]    a positive integer specifying the number of cores to use when executing the MCMC sampling.
    ###
    ### [outpdf]   a character string specifying the name of a file in which MCMC sampling results are saved.
    ###            The saved results will be written to a PDF file in the current working directory. 
    ###
    ### [outfile]  a character string specifying the name of a file in which posterior samples are saved.
    ###            The saved samples will be written to a CSV file in the current working directory. 
    ###
    ### Output of function mcmcSAM():
    ### a list containing posterior draws for each parameter.
    ### 
    ### Dependence: the R software; the R package "rstan".
    ###
    ### Author: Peng Jun, 2019.01.07.
    ###------------------------------------------------------------------------------------------------------------

    mcmcSAMseq <- function(EDdata, DRdata, iN, priorAge, Sdata, model, sDRc,   
                           iflog=TRUE, ordered=TRUE, iter=5e3, warmup=2e3, thin=1, 
                           chains=3, cores=1, outpdf="mcmcSAMseq", outfile="mcmcSAMseq") {
        ###
        if (ncol(EDdata)!=2L) 
            stop("Argument [EDdata] should be a two-column matrix!")
        ###
        if (any(as.numeric(EDdata[,2L])<=0)) 
            stop("The second column of argument [EDdata] should not contain values below zero!")
        ###
        if (ncol(DRdata)!=2L) 
            stop("Argument [DRdata] should be a two-column matrix!")
        ###
        if (any(as.numeric(DRdata[,2])<=0)) 
            stop("The second column of argument [DRdata] should not contain values below zero!")
        ###
        if (!is.numeric(iN)) 
            stop("Argument [iN] should be of type numeric!")
        ###
        if (sum(iN)!=nrow(EDdata)) 
            stop("Inconsistent number of ED values in arguments [EDdata] and [iN]!")
        ###
        if (length(iN)!=nrow(DRdata)) 
            stop("Incorrect length of argument [iN]!")
        ###
        if (length(priorAge)!=2L) 
            stop("Argument [priorAge] should be a two-element vector!")
        ###
        if (!is.vector(Sdata, mode="numeric")) 
            stop("Argument [Sdata] should be a vector!")
        ##
        if (length(Sdata)!=nrow(DRdata))
            stop("Incorrect length for argument [Sdata]!")
        ###
        if (!all(model %in% c("minimum","central","maximum"))) 
            stop("Argument [model] should be chosen from 'minimum', 'central', and 'maximum'!")
        ###
        if (length(model)!=nrow(DRdata)) 
            stop("Incorrect length for argument [model]!")
        ###
        if (!is.numeric(sDRc)) 
            stop("Argument [sDRc] should be of type numeric!")
        ###
        if (!is.logical(iflog)) 
            stop("Argument [iflog] should be of type logical!")
        ###
        if (!is.logical(ordered)) 
            stop("Argument [ordered] should be of type logical!")
        ###
        if (!is.numeric(iter)) 
            stop("Argument [iter] should be of type numeric!")
        ###
        if (!is.numeric(warmup)) 
            stop("Argument [warmup] should be of type numeric!")
        ###
        if (!is.numeric(thin)) 
            stop("Argument [thin] should be of type numeric!")
        ###
        if (thin<1L) 
            stop("Argument [thin] should not below 1!")
        ###
        if (warmup>=floor(iter/thin)) 
            stop("Argument [warmup] should not exceed the ratio of [iter] to [thin]!")
        ###
        if (!is.numeric(chains)) 
            stop("Argument [chains] should be of type numeric!")
        ###
        if (chains<1L || chains>7L) 
            stop("Argument [chains] should lie between 1 and 7!")
        ###
        if (!is.numeric(cores))
            stop("Argument [cores] should be of type numeric!")
        ###
        if (!is.character(outpdf) && !is.null(outpdf)) 
            stop("Argument [outpdf] should be of type character or NULL!")
        ###
        if (!is.character(outfile) && !is.null(outfile)) 
            stop("Argument [outfile] should be of type character or NULL!")
        ###
        ###
        n <- nrow(EDdata)
        m <- nrow(DRdata)
        ###
        if (m<2L) 
            stop("At least two sets of De distribution should be provided!")
        ###
        ###
        fSAM <- function(chainID) {
            p <- runif(n=m, min=0.001, max=0.999)
            ###
            age <- runif(n=m, min=priorAge[1L], max=priorAge[2L])
            if (ordered==TRUE) age <- sort(age)
            ###  
            sigma <- runif(n=m, min=0.0, max=5.0)
            ###
            DR <- rnorm(n=m, mean=DRdata[,1L], sd=DRdata[,2L])  
            ###
            if (is.vector(Sdata, mode="numeric")) {   
                output <- list("p"=p, "age"=age, "sigma"=sigma, "DR"=DR)
            } else if (is.matrix(Sdata)) {
                S <- rnorm(m, mean=Sdata[,1L], sd=Sdata[,2L]) 
                output <- list("p"=p, "age"=age, "sigma"=sigma, "DR"=DR, "S"=S)
            } # end if.
            ###
            return(output)
            ###
        } # end function fSAM.
        ###
        ###
        idx <- cbind(cumsum(iN)-iN+1L, cumsum(iN))
        ###
        startPAR <- lapply(seq(chains), function(id) fSAM(chainID=id))
        ###
        mdl <- vector(length=m) 
        ###
        for (i in seq(m)) {
            if (model[i]=="minimum") {
                mdl[i] <- -1L
            } else if (model[i]=="central") {
                mdl[i] <- 1L
            } else if (model[i]=="maximum") {
                mdl[i] <- -3L
            } # end if.
        } # end for.
        ###
        iflog <- ifelse(iflog==TRUE, 1L, 0L)
        ###
        ###
        if (ordered==FALSE) {
            ###
            fit <- stan(file="SAMseq0.stan", 
                        data=list(n=n, m=m, EDdata=EDdata, DRdata=DRdata, idx=idx,
                        priorAge=priorAge, S=Sdata, mdl=mdl, sDRc=sDRc, iflog=iflog), 
                        iter=iter, warmup=warmup, thin=thin, chains=chains, init=startPAR,
                        cores=cores, control=list(adapt_delta=0.85))
            ###
        } else if (ordered==TRUE) {
            ###
            fit <- stan(file="SAMseq1.stan", 
                        data=list(n=n, m=m, EDdata=EDdata, DRdata=DRdata, idx=idx,
                        priorAge=priorAge, S=Sdata, mdl=mdl, sDRc=sDRc, iflog=iflog), 
                        iter=iter, warmup=warmup, thin=thin, chains=chains, init=startPAR,
                        cores=cores, control=list(adapt_delta=0.85))
            ###
        } # end if.
        ###
        ###
        extractPAR <- c()
        ###
        for (i in seq(m)) {
            if (model[i]=="minimum" || model[i]=="maximum") {
                addPAR <- c(paste("p[", i, "]", sep=""), paste("cDose[", i, "]", sep=""), 
                            paste("sigma[", i, "]", sep=""), paste("age[", i, "]", sep=""))
                extractPAR <- c(extractPAR, addPAR) 
            } else if (model[i]=="central") {
                addPAR <- c(paste("cDose[", i, "]", sep=""), paste("sigma[", i, "]", sep=""), 
                            paste("age[", i, "]", sep=""))
                extractPAR <- c(extractPAR, addPAR) 
            } # end if.
        } # end for.
        ###
        mcmcSamples <- attr(fit,"sim")$samples
        ###
        samples <- vector(mode="list", length=length(extractPAR))
        ###
        N1 <- attr(fit,"sim")$warmup2[1L] 
        N2 <- attr(fit,"sim")$n_save[1L]
        ###
        outputMAT <- c()
        ###
        for (i in seq(extractPAR)) {
            samples[[i]] <- matrix(nrow=N2-N1, ncol=chains)
            colnames(samples[[i]]) <- paste(extractPAR[i], ".chain.", seq(chains), sep="")
            ###
            for (j in seq(chains)) {
                msp <- mcmcSamples[[j]]
                samples[[i]][,j] <-  (msp[[match(extractPAR[i],names(msp))]])[(N1+1L):N2]
            } # end for.
            ###
            if (chains==1L) samples[[i]] <- as.numeric(samples[[i]])
            ###
            outputMAT <- cbind(outputMAT, samples[[i]])
        } # end for. 
        ###
        names(samples) <- extractPAR
        ###
        if (!is.null(outfile)) write.csv(outputMAT, file=paste(outfile, ".csv", sep=""))
        ###
        ###
        if (!is.null(outpdf)) { 
            colors <- c("deepskyblue", "purple", "red", "blue", 
                        "springgreen", "orange", "brown")
            ###
            pdf(paste(outpdf, ".pdf", sep="")) 
            layout(rbind(c(1L,2L),c(3L,4L)), respect=TRUE)
            ###
            for (i in seq(extractPAR)) {
                isap <- samples[[i]]
                ###
                ### Draw a trace plot. 
                plot(NA, NA, main="A: trace plot", xlim=c(N1+1L,N2), ylim=range(isap), 
                     xlab="Number of iterations (after thining)", 
                     ylab=paste("Value of ", extractPAR[i], sep=""))
                ###
                if (chains>=2L) {
                    for (j in seq(chains)) {                  
                        points(x=(N1+1L):N2, y=isap[,j], type="l", col=colors[j], lwd=1.1)
                    } # end for.
                    points(x=(N1+1L):N2, y=rowMeans(isap), type="l", col="grey50", lwd=1.5)
                } else {
                    points(x=(N1+1L):N2, y=isap, type="l", col=colors[j], lwd=1.1)
                } # end if.
                ###
                ### Draw a probability density plot.
                disap <- density(c(isap))
                plot(disap, type="n", main="B: density plot", xlab=extractPAR[i], 
                     ylab=paste("Density of ", extractPAR[i],sep=""),lwd=2.0)
                polygon(cbind(disap$x,disap$y), col="skyblue")  
                ###
                ### Draw a acf plot.
                if (chains>=2L) {
                    acfv <- rowMeans( apply(isap, MARGIN=2L, function(x) 
                                      as.numeric(acf(x, plot=FALSE)$acf) ) )
                } else {
                    acfv <- acf(isap, plot=FALSE)$acf
                } # end if.
                ###
                lagv <- seq(acfv) - 1L
                plot(x=lagv, y=acfv, type="h", lwd=2.0,  main="C: acf plot", 
                     xlab="Lag", ylab=paste("Average ACF of ",extractPAR[i], sep=""))
                ###
                ### Diagnostic statistics. 
                if (chains>=2L) {                
                    GR <- function(x) {
                        n <- nrow(x)
                        W <- mean(apply(x, MARGIN=2L, var))                 
                        B <- n/(chains-1.0)* sum((colMeans(x)-mean(x))^2L)
                        VarTheta <- (1.0 - 1.0/n)*W + 1.0/n*B
                        return( sqrt(VarTheta/W))
                    } # end function GR.
                    ###
                    xvec <- seq(from=50L, to=N2-N1, by=(N2-N1-50.0)/299.0)
                    yvec <- vector(length=length(xvec))
                    for (j in seq(xvec))  yvec[j] <- GR(isap[1L:xvec[j],,drop=FALSE])
                    ###
                    plot(xvec+N1, yvec, ylim=c(0.99*min(yvec),1.02*max(yvec)), type="o", lwd=1.5, cex=0.3, 
                         pch=21, bg="grey", main="D: diagnostic plot", xlab="Number of iterations (after thining)",
                         ylab=paste("Shrink factor of ", extractPAR[i], sep=""))
                    abline(h=1.0, col="red", lwd=2.0, lty="dashed")
                } else {
                    plot(x=0.5, y=0.5, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
                } # end if.
                ###
            } # end for.
            ###
            dev.off() 
            ###
        } # end if.
        ###
        print(fit, pars=extractPAR, digits_summary=4)
        ###
        pars <- summary(fit, pars=extractPAR, digits_summary=4)$summary
        ###
        output <- list("samples"=samples, "pars"=pars)
        ###
        invisible(output)
        ###
    } # end function mcmcSAMseq.
    ###------------------------------------------------------------------------------------------
    ###