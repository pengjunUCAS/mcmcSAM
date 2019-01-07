
 
    ### 
    library(rstan)
    ###
    ###------------------------------------------------------------------------------------------------------------
    ### Function mcmcSAM() is used for obtaining posterior distributions of parameters of statistical age models  
    ### (SAM) for a set of equivlalent dose (ED) values from a sample using a Markov Chain Monte Carlo (MCMC) method.
    ###
    ### SAM that can be implemented included: (1) Central Age Model, 
    ### (2) Minimum Age Model(3-parameter), (3) Maximum Age Model (3-parameter).
    ###
    ### Arguments in function mcmcSAM():
    ### [EDdata]    a two-column matrix containing equivalent dose (ED) and associated errors.
    ###
    ### [DRdata]    a two-element vector containing dose rate (DR) and associated error.
    ###             For example, if the dose rate of the sample is 3.0¡À0.21 Gy/ka, then DRdata=c(3.0, 2.1).      
    ###
    ### [priorAge]  a two-element vector containing ranges (i.e.,the lower and upper bounds of a uniform distribution)
    ###             of the prior for the age. For example, if the prior of age is assummed to follow a uniform
    ###             distribution in the interval [1,100], then priorAge=c(0,100).
    ###
    ### [Sdata]     a real value containing additional sigma value. For example, if the extra error that want
    ###             be added in quadrature to the relative standard errors of ED values is 20%,then Sdata=0.2.
    ###
    ### [model]     a character string indicating the model to be fitted£¬i.e., 
    ###             (1) minimum¡ªthe 3-parameter minimum age model;
    ###             (2) central¡ªthe central age model; 
    ###             (3) maximum¡ªthe 3-parameter maximum age model.
    ###
    ### [iflog]     a logical value indicating if the model should be fitted using log-scaled EDdata or not. 
    ###
    ### [iter]      a positive integer specifying the number of iterations for each chain.
    ### 
    ### [warmup]    a positive integer specifying the number of warmup (i.e.,burn-in) iterations for each chain.
    ###             For example, if warmup=1e4 then samples simulated from the first 1e4 iterations for each 
    ###             chain will be disgarded. 
    ###
    ### [thin]      a positive integer specifying the period for saving samples. For example, if thin=5
    ###             then only samples 1 out of 5 iterations will be reserved from the remaining samples.
    ### 
    ### [chains]    a positive integer specifying the number of Markov chains. The number of chain is 
    ###             allowed to be chosen from 1 to 7.
    ###
    ### [cores]     a positive integer specifying the number of cores to use when executing the MCMC sampling.
    ###
    ### [outpdf]    a character string specifying the name of a file in which MCMC sampling results are saved.
    ###             The saved results will be written to a file in the current working directory. 
    ###
    ### [outfile]   a character string specifying the name of a file in which posterior samples are saved.
    ###             The saved samples will be written to a CSV file in the current working directory. 
    ###
    ### Output of function mcmcSAM():
    ### a list containing posterior draws for each parameter.
    ### 
    ### Dependence: the R software; the R package "rstan".
    ###
    ### Author: Peng Jun, 2019.01.07.
    ###------------------------------------------------------------------------------------------------------------
     
    mcmcSAM <- function(EDdata, DRdata, priorAge, Sdata, model,  
                        iflog=TRUE, iter=1e4, warmup=2e3, thin=1, chains=3,  
                        cores=1, outpdf="mcmcSAM", outfile="mcmcSAM") {
        ###
        if (ncol(EDdata)!=2L) 
            stop("Argument [EDdata] should be a two-column matrix!")
        ###
        if (any(as.numeric(EDdata[,2L])<=0)) 
            stop("The second column of argument [EDdata] should not contain values below zero!")
        ###
        if (length(DRdata)!=2L) 
            stop("Argument [DRdata] should be a two-element vector!")
        ###
        if (DRdata[2L]<=0) 
            stop("The second element in argument [DRdata] should exceed zero!")
        ###
        if (length(priorAge)!=2L) 
            stop("Argument [priorAge] should be a two-element vector matrix!")
        ###
        if (length(Sdata)!=1L) 
            stop("Argument [Sdata] should be an one-element vector!")
        ###
        if (!model %in% c("minimum","central","maximum")) 
            stop("Argument [model] should be chosen from 'minimum', 'central', or 'maximum'!")
        ###
        if (!is.logical(iflog)) 
            stop("Argument [iflog] should be of type logical!")
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
        if (!is.numeric(chains)) 
            stop("Argument [chains] should be of type numeric!")
        ###
        if (chains<1L || chains>7L) 
            stop("Argument [chains] should lie between 1 and 7!")
        ###
        if (warmup>=floor(iter/thin)) 
            stop("Argument [warmup] should not exceed the ratio of [iter] to [thin]!")
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
        ed <- as.numeric(EDdata[,1L])
        ###
        ###
        fSAM <- function(chainID) {
            p <- runif(n=1L, min=0.001, max=0.999)
            ###
            age <- runif(n=1L, min=priorAge[1L], max=priorAge[2L])
            ###
            sigma <- runif(n=1L, min=0.0, max=5.0)
            ###
            DR <- rnorm(n=1L, mean=DRdata[1L], sd=DRdata[2L])
            ###
            output <- list("p"=p, "age"=age, "sigma"=sigma, "DR"=DR)
            ###
            return(output)
            ###
        } # end function fSAM.
        ###
        ###
        startPAR <- lapply(seq(chains), function(id) fSAM(chainID=id))
        ###
        mdl <- if (model=="minimum") {
            -1L
        } else if (model=="central") {
            1L
        } else if (model=="maximum") {
            -3L
        } # end if.
        ###
        iflog <- ifelse(iflog==TRUE, 1L, 0L)
        ###
        ###
        fit <- stan(file="SAM0.stan", 
                    data=list(n=n,EDdata=EDdata,DRdata=DRdata,priorAge=priorAge, 
                    S=Sdata,mdl=mdl,iflog=iflog), iter=iter, warmup=warmup, 
                    thin=thin, chains=chains, init=startPAR, cores=cores, 
                    control=list(adapt_delta=0.85))
        ###
        ###
        if (model=="minimum" || model=="maximum") {
            extractPAR <- c("p","cDose", "sigma", "age")
        } else if (model=="central") {
            extractPAR <- c("cDose", "sigma", "age")
        } # end if.
        ###
        mcmcSamples <- attr(fit,"sim")$samples
        samples <- vector(mode="list", length=length(extractPAR))
        ###
        N1 <- attr(fit,"sim")$warmup2[1L] 
        N2 <- attr(fit,"sim")$n_save[1L]
        ###
        outputMAT <- c()
        ###
        for (i in seq(extractPAR)) {
            samples[[i]] <- matrix(nrow=N2-N1, ncol=chains)
            colnames(samples[[i]]) <- paste(extractPAR[i],".chain.", seq(chains), sep="")
            ###
            for (j in seq(chains)) {
                msp <- mcmcSamples[[j]]
                samples[[i]][,j] <-  msp[[match(extractPAR[i],names(msp))]][(N1+1L):N2]
            } # end for.
            ###
            if (chains==1L) samples[[i]] <- as.numeric(samples[[i]])
            ###
            outputMAT <- cbind(outputMAT, samples[[i]])
        } # end for. 
        ###
        names(samples) <- extractPAR
        ###
        if (!is.null(outpdf)) write.csv(outputMAT, file=paste(outfile,".csv",sep=""))
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
                ###
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
                    plot(xvec+N1, yvec, ylim=c(0.99*min(yvec),1.02*max(yvec)), type="o", 
                         lwd=1.5, cex=0.3, pch=21, bg="grey", main="D: diagnostic plot", 
                         xlab="Number of iterations (after thining)",
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
        ###
        print(fit, pars=extractPAR, digits_summary=4)
        ###
        pars <- summary(fit, pars=extractPAR, digits_summary=4)$summary
        ###
        output <- list("samples"=samples, "pars"=pars)
        ###
        invisible(output)
        ###
    } # end function mcmcSAM.
    ###------------------------------------------------------------------------------------------
    ###