


    // Estimate parameters for a set of De distribution according to: 
    // (1) Central Dose (Age) Model, 
    // (2) Minimum Dose (Age) Model,
    // (3) Maximum Dose (Age) Model.
    //
    // Peng Jun, 2019.01.07.

    functions{

        real SAM0_log(matrix EDdata, real p, real age, real sigma, 
                      real DR, real S, int n, int mdl, int iflog){

            vector[n] x;
            vector[n] y;
            vector[n] gamma0;
            vector[n] sigma0;
            vector[n] prb1;
            vector[n] prb2;
            vector[n] prb3;
            vector[n] fp;
            real slfp;   
            real mu;
            real gamma;   
            real pi;

            pi = 3.141592657;

            slfp = 0;

            if (mdl==1) {

                // Implement the Central Age Model (CAM).

                if (iflog==0) {

                    mu = age*DR;

                } else if (iflog==1) {

                    mu = log(age*DR);

                } // end if.

                for (i in 1:n){

                    if (iflog==0) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2)) * EDdata[i,1];
                        y[i] = EDdata[i,1];
                        
                    } else if (iflog==1) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2));
                        y[i] = log(EDdata[i,1]);

                    } // end if.

                    fp[i] = exp(-pow(p-0.5,2))/sqrt(2.0*pi*(pow(x[i],2)+pow(sigma,2)))*
                            exp(-0.5*pow(y[i]-mu,2)/(pow(x[i],2)+pow(sigma,2)));

                } // end for.

            } else if (mdl==-1) {

                // Implement the Minimum Age model (MAM3).

                if (iflog==0) {

                    gamma = age*DR;
                    
                } else if (iflog==1) {

                    gamma = log(age*DR);

                } // end if.

                for (i in 1:n){

                    if (iflog==0) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2))*EDdata[i,1];
                        y[i] = EDdata[i,1];

                    } else if (iflog==1) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2));
                        y[i] = log(EDdata[i,1]);

                    } // end if.

                    gamma0[i] = (gamma/pow(sigma,2)+y[i]/pow(x[i],2))/ 
                                (1.0/pow(sigma,2)+1.0/pow(x[i],2));
     
                    sigma0[i] = 1.0/sqrt(1.0/pow(sigma,2)+1.0/pow(x[i],2));

                    prb1[i] = p/sqrt(2.0*pi)/x[i]*exp(-0.5*pow(y[i]-gamma,2)/pow(x[i],2));

                    prb2[i] = (1.0-p)/sqrt(2.0*pi)/sqrt(pow(sigma,2)+pow(x[i],2))*
                              (1.0-Phi((gamma-gamma0[i])/sigma0[i]))*2.0;

                    prb3[i] = exp(-0.5*pow(y[i]-gamma,2)/(pow(sigma,2)+pow(x[i],2)) );

                    fp[i] = prb1[i] + prb2[i]*prb3[i];

                } // end for.

            } else if (mdl==-3) {

                // Implement the Maximum Age model (MAM3).

                if (iflog==0) {
               
                    gamma = age*DR;
                    
                } else if (iflog==1) {

                    gamma = log(age*DR);    

                } // end if.

                for (i in 1:n){

                    if (iflog==0) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2))*EDdata[i,1];
                        y[i] = EDdata[i,1];

                    } else if (iflog==1) {

                        x[i] = sqrt(pow(EDdata[i,2]/EDdata[i,1],2) + pow(S,2));
                        y[i] = log(EDdata[i,1]);

                    } // end if.

                    gamma0[i] = (gamma/pow(sigma,2)+y[i]/pow(x[i],2))/ 
                                (1.0/pow(sigma,2)+1.0/pow(x[i],2));
     
                    sigma0[i] = 1.0/sqrt(1.0/pow(sigma,2)+1.0/pow(x[i],2));

                    prb1[i] = p/sqrt(2.0*pi)/x[i]*exp(-0.5*pow(y[i]-gamma,2)/pow(x[i],2));

                    prb2[i] = (1.0-p)/sqrt(2.0*pi)/sqrt(pow(sigma,2)+pow(x[i],2))*
                              (Phi((gamma-gamma0[i])/sigma0[i]))*2.0;

                    prb3[i] = exp(-0.5*pow(y[i]-gamma,2)/(pow(sigma,2)+pow(x[i],2)) );

                    fp[i] = prb1[i] + prb2[i]*prb3[i];

                } // end for.

            } // end if.

            for (i in 1:n) { 

                if( fp[i]<1e-13 ) {  fp[i] = 1e-13; } // end if.

                slfp = slfp + log(fp[i]); 

            } // end for.
       
            return(slfp);

        } // end function SAM0_log.

    } // end <functions> block.


    data{
     
        int n;
        matrix[n,2] EDdata;
        vector[2] DRdata;
        vector[2] priorAge;
        real S;
        int mdl;
        int iflog;

    } // end <data> block.


    parameters{

        real<lower=0, upper=1> p;
        real<lower=0> age;
        real<lower=0> sigma; 
        real<lower=0> DR;

    } // end <parameter> block.


    transformed parameters{

        real cDose;

        cDose = age*DR;

    } // end <transformed parameters> block.


    model {
  
        p ~ uniform(0.001, 0.999); 

        age ~ uniform(priorAge[1], priorAge[2]);

        //sigma ~ cauchy(0, 5) T[0,10]; 

        sigma ~ uniform(0, 5); 

        DR ~ normal(DRdata[1], DRdata[2]); 
        
        EDdata ~ SAM0(p, age, sigma, DR, S, n, mdl, iflog);

    } // end <model> block.


