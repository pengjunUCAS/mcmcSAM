


    // Estimate parameters in a sequence of De distributions according to:
    // (1) Central Age Model, 
    // (2) Minimum Age Model, 
    // (3) Maximum Age Model,
    // With stratigraphic constraint. 
    //
    // Peng Jun, 2019.01.07.

    functions{

        real SAMseq0_log(matrix EDdata, vector p, vector age, vector sigma, vector DR, 
                         vector S, int n, int m, int[,] idx, int[] mdl, int iflog){

            vector[n] x;
            vector[n] y;
            vector[n] gamma0;
            vector[n] sigma0;
            vector[n] prb1;
            vector[n] prb2;
            vector[n] prb3;
            vector[n] fp;
            vector[m] slfp;
            real mu;
            real gamma;
            real pi;
            int ia;
            int ib;

            pi = 3.141592657;

            for (j in 1:m) {

                slfp[j] = 0;

                ia = idx[j,1];
                ib = idx[j,2];

                if (mdl[j]==1) {

                    // Implement the Central Age Model (CAM).

                    if (iflog==0) {

                        mu = age[j]*DR[j];
 
                    } else if (iflog==1) {
 
                        mu = log(age[j]*DR[j]);

                    } // end if.

                    for (i in ia:ib){

                        if (iflog==0) {

                            x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2)) * EDdata[i,1];
                            y[i] = EDdata[i,1];                        

                        } else if (iflog==1) {

                            x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2));
                            y[i] = log(EDdata[i,1]);

                        } // end if.

                        fp[i] = exp(-pow(p[j]-0.5,2))/sqrt(2.0*pi*(pow(x[i],2)+pow(sigma[j],2)))*
                                exp(-0.5*pow(y[i]-mu,2)/(pow(x[i],2)+pow(sigma[j],2)));

                    } // end for.

                } else if (mdl[j]==-1) {

                    // Implement the Minimum Age model (MAM3).

                    if (iflog==0) {

                        gamma = age[j]*DR[j];
                    
                    } else if (iflog==1) {

                        gamma = log(age[j]*DR[j]);

                    } // end if.
 
                    for (i in ia:ib){

                        if (iflog==0) {

                            x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2)) * EDdata[i,1];
                            y[i] = EDdata[i,1];

                        } else if (iflog==1) {

                            x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2));
                            y[i] = log(EDdata[i,1]);

                        } // end if.

                        gamma0[i] = (gamma/pow(sigma[j],2)+y[i]/pow(x[i],2))/ 
                                    (1.0/pow(sigma[j],2)+1.0/pow(x[i],2));
     
                        sigma0[i] = 1.0/sqrt(1.0/pow(sigma[j],2)+1.0/pow(x[i],2));

                        prb1[i] = p[j]/sqrt(2.0*pi)/x[i]*exp(-0.5*pow(y[i]-gamma,2)/pow(x[i],2));

                        prb2[i] = (1.0-p[j])/sqrt(2.0*pi)/sqrt(pow(sigma[j],2)+pow(x[i],2))*
                                  (1.0-Phi((gamma-gamma0[i])/sigma0[i]))*2.0;

                        prb3[i] = exp(-0.5*pow(y[i]-gamma,2)/(pow(sigma[j],2)+pow(x[i],2)) );

                        fp[i] = prb1[i] + prb2[i]*prb3[i];

                    } // end for.

                } else if (mdl[j]==-3) {

                         // Implement the Maximum Age model (MAM3).

                    if (iflog==0) {

                        gamma = age[j]*DR[j];
                    
                    } else if (iflog==1) {

                        gamma = log(age[j]*DR[j]);

                    } // end if.
 
                    for (i in ia:ib){

                        if (iflog==0) {

                            x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2)) * EDdata[i,1];
                            y[i] = EDdata[i,1];

                        } else if (iflog==1) {

                            x[i] = sqrt( pow(EDdata[i,2]/EDdata[i,1],2) + pow(S[j],2) );
                            y[i] = log(EDdata[i,1]);

                        } // end if.

                        gamma0[i] = (gamma/pow(sigma[j],2)+y[i]/pow(x[i],2))/ 
                                    (1.0/pow(sigma[j],2)+1.0/pow(x[i],2));
     
                        sigma0[i] = 1.0/sqrt(1.0/pow(sigma[j],2)+1.0/pow(x[i],2));

                        prb1[i] = p[j]/sqrt(2.0*pi)/x[i]*exp(-0.5*pow(y[i]-gamma,2)/pow(x[i],2));

                        prb2[i] = (1.0-p[j])/sqrt(2.0*pi)/sqrt(pow(sigma[j],2)+pow(x[i],2))*
                                  (Phi((gamma-gamma0[i])/sigma0[i]))*2.0;

                        prb3[i] = exp(-0.5*pow(y[i]-gamma,2)/(pow(sigma[j],2)+pow(x[i],2)) );

                        fp[i] = prb1[i] + prb2[i]*prb3[i];

                    } // end for.

                } // end if.

                for (i in ia:ib) { 

                    if( fp[i]<1e-13 ) {  fp[i] = 1e-13; } // end if.

                    slfp[j] = slfp[j] + log(fp[i]); 

                } // end for.

            } // end for.
       
            return(sum(slfp));

        } // end function SAMseq0_log.

    } // end <functions> block.



    data{
     
        int n;
        int m;
        matrix[n,2] EDdata;
        matrix[m,2] DRdata;
        int idx[m,2];
        vector[2] priorAge;
        vector[m] S;
        int mdl[m];
        real sDRc;
        int iflog;

    } // end <data> block.



    parameters{

        vector<lower=0, upper=1>[m] p;
        positive_ordered[m] age;
        vector<lower=0>[m] sigma; 
        vector<lower=0>[m] DR;

    } // end <parameter> block.



    transformed parameters{

        vector[m] cDose;

        for (j in 1:m) {

            cDose[j] = age[j]*DR[j];

        } // end for.

    } // end <transformed parameters> block.



    model{
  
        p ~ uniform(0.001,0.999); 

        for (j in 1:m) { 

            age[j] ~ uniform(priorAge[1],priorAge[2]);

            sigma[j] ~ uniform(0, 5); 

            DR[j] ~ normal(DRdata[j,1], sqrt(pow(DRdata[j,2],2)+pow(sDRc,2))); 

        } // end for.     
        
        EDdata ~ SAMseq0(p, age, sigma, DR, S, n, m, idx, mdl, iflog);

    } // end <model> block.



