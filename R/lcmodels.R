risk.bach.10 <- function(data){     
   bach.formula <- function(outcome){
    
      # age = continuous
      # female = binary indicator
      # cpd = cigarettes per day
      # smkyears = years smoked
      # qtyears = years quit
      # asb = asbestos exposure binary indicator
          
      #Formula from appendix of Bach's paper minus coefficients
      formula.string <- "female+
      asb+
      age+
      I(I((age-53.459001)^3)*I(age>53))+
      I(I((age-61.954825)^3)*I(age>61))+
      I(I((age-70.910335)^3)*I(age>70))+
      qtyears+
      I((qtyears)^3)+
      I(I((qtyears-0.50513347)^3)*I(qtyears>0))+
      I(I((qtyears-12.295688)^3)*I(qtyears>12))+
      smkyears+
      I(I((smkyears-27.6577)^3)*I(smkyears>27))+
      I(I((smkyears-40)^3)*I(smkyears>40))+
      I(I((smkyears-50.910335)^3)*I(smkyears>50))+
      cpd+
      I(I((cpd-15)^3)*I(cpd>15))+
      I(I((cpd-20.185718)^3)*I(cpd>20))+
      I(I((cpd-40)^3)*I(cpd>40))"
      				       
      formula(paste(outcome,"~",formula.string,collapse=""))
   }
   # Coefficients from diagnosis of lung cancer appendix
   bach.coef <- c(intercept =-9.7960571,
                female =-0.05827261,
                asb =0.2153936,
                age=0.070322812,
                age2=-0.00009382122,
                age3=0.00018282661,
                age4=-0.000089005389,
                qtyears=-0.085684793,
                qtyears2=0.0065499693,
                qtyears3=-0.0068305845,
                qtyears4=0.00028061519,
                smkyears=0.11425297,
                smkyears1=-0.000080091477,
                smkyears2=0.00017069483,
                smkyears3=-0.000090603358,
                cpd =0.060818386,
                cpd1 =-0.00014652216,
                cpd2 =0.00018486938,
                cpd3 =-0.000038347226)
   
   # Coefficients from probability of death from something other than lung cancer
   bach.coef.mort <- c(Intercept=-7.2036219, female=-0.49042298, asb=0.06084611, 
   						         age=0.099168033, age1=6.2174577e-06, age2=-1.2115774e-05, age3=5.8983164e-06, 
   						         qtyears=-0.023358962, qtyears2=0.0019208669, qtyears3=-0.0020031611, qtyears4=8.2294194e-05, 
   						         smkyears=0.020041889, smkyears1=6.5443781e-06, smkyears2=-1.3947696e-05, smkyears3=7.4033175e-06, 
   						         cpd=0.015490665, cpd1=-1.737645e-05, cpd2=2.1924149e-05, cpd3=-4.5476985e-06)

  vars <- c("female",
  				   "age",
  				   "cpd",
            "smkyears",
            "qtyears",
            "asb")

  var.not.in <- names(data)[!(names(data)%in%vars) & sapply(data, function(x) sum(is.na(x)))==0][1] # SELECT ARBITRARY VAR FOR OUTCOME
  
  P.lung <- matrix(0, nrow(data), ncol=10)
  P.death <- matrix(0, nrow(data), ncol=10)
  
  for(i in 1:10){
    # Calculate RRs and probabilities
  	rr.lung <- model.matrix(bach.formula(var.not.in), data)%*%bach.coef
	  rr.death <- model.matrix(bach.formula(var.not.in), data)%*%bach.coef.mort
	  P.lung[complete.cases(subset(data,select=vars)),i] <- .99629^exp(rr.lung)
	  P.death[complete.cases(subset(data,select=vars)),i] <- 0.9917663^exp(rr.death)
    # Update age, smoking and quit years for next cycle
	  data$age <- ifelse(is.na(data$age),NA,data$age+1)
	  data$smkyears <- ifelse(is.na(data$smkyears),NA,data$smkyears+ifelse(data$qtyears==0,1,0)) # Add smoking year for current smokers
	  data$qtyears <- ifelse(is.na(data$qtyears),NA,data$qtyears+ifelse(data$qtyears==0,0,1)) # Add quit year for former smokers
   }
	
	# CUMULATIVE PROBABILITIES OF NO CANCER AND NO DEATH
   
   C.lung <- P.lung
   # Delete last column, put 1s in first column and shift other columns over by one
   C.lung[,2:10] <- C.lung[,1:9]
   C.lung[,1] <- 1
   C.lung <- t(apply(C.lung, 1, cumprod))
   
   C.mort <- P.death
   C.mort[,2:10] <- C.mort[,1:9]
   C.mort[,1] <- 1
   C.mort <- t(apply(C.mort, 1, cumprod))
	 
   #The final risk calculation
   ifelse(rowSums((1-P.lung)*P.death*C.lung*C.mort)==0,NA,rowSums((1-P.lung)*P.death*C.lung*C.mort))
} 



risk.hoggart <- function(data){
   weibull.one.year <- function(t,lambda, gamma, HR=1){     
      shape <- exp(gamma)
      scale <- exp(lambda)
     
      #Time zero set to 35
      t1 <- t-35+1
      t0 <- t-35
      #Calculate HR times shape/scale difference for one time point
      diff <- HR*((t0/scale)^shape-(t1/scale)^shape)
      #Return exponentiated difference  
      exp(diff)
   }
   
   hoggart.coef <- function(current.age, age, former, smkyears, hr.lung, hr.mort){
      #For current smokers, cut points for age started smoking
      age.current <- c(-1,18,20,22,24,26,28,30,100)
      #For former smokers, cut points for age started smoking
      age.former <- c(-1,22,26,100)
      #For former smokers, cut points for duration of smoking
      duration <- c(-1,20,30,36,42,100)
      
      #Shape and scale parameters from Table 2 of Hoggart paper
      lambda.lung.current <- c(3.819, 4.056, 4.23, 4.339, 4.38, 4.567, 4.506, 4.504)
      gamma.lung.current <- c(0.999, 1.071, 1.298, 1.518, 1.679, 1.517, 1.615, 1.684)
      lambda.mort.current <- c(3.69, 3.774, 3.859, 3.944, 3.979, 4.049, 4.096, 4.069)
      gamma.mort.current <- c(1.22, 1.312, 1.56, 1.775, 1.925, 1.854, 1.838, 1.912)
      
      lambda.lung.former <- c(4.987, 4.723, 4.321, 5.179, 4.786, 4.651, 4.654, 4.366, 5.563, 4.68, 4.491, 4.241, 4.177)
      gamma.lung.former <- c(0.75, 0.819, 1.032, 1.165, 1.353, 1.46, 1.318, 1.662, 1.088, 1.705, 1.879, 2.336, 2.574)
      lambda.mort.former <- c(3.754, 3.75, 3.8, 4.008, 3.975, 3.954, 3.928, 3.982, 4.123, 4.058, 4.02, 3.999, 4.041)
      gamma.mort.former <- c(1.21, 1.511, 1.669, 1.621, 1.742, 1.835, 2.129, 2.15, 1.769, 1.865, 2.148, 2.36, 2.49)
      
      
      if(former==1){#For former smokers
      	strata.age <- as.numeric(cut(age, age.former, lab=1:3)) # <=22, 22-26, >26
      	strata.duration <- as.numeric(cut(smkyears, duration, lab=1:5)) # Note that strata.age=1 only sees cut points at 20&30
      	strata <- ifelse(strata.age==1,ifelse(strata.duration>=3,3,strata.duration),
      									ifelse(strata.age==2,strata.duration+3,strata.duration+8)) # ifelse(test, yes, no)
        #Now assign correct lambdas and gammas to people in each strata
      	lambda.lung <- lambda.lung.former[strata]
      	gamma.lung <- gamma.lung.former[strata]
      	lambda.mort <- lambda.mort.former[strata]
      	gamma.mort <- gamma.mort.former[strata]
      }
      else{#For current smokers
      	strata <- as.numeric(cut(age, age.current, lab=1:8)) #Cut points for age started smoking
        #Assign correct lambdas and gammas to people in each strata
      	lambda.lung <- lambda.lung.current[strata]
      	gamma.lung <- gamma.lung.current[strata]
      	lambda.mort <- lambda.mort.current[strata]
      	gamma.mort <- gamma.mort.current[strata]
      }
      	#Return weibull estimates for probability of getting lung cancer after one year
      (1-weibull.one.year(current.age,lambda=lambda.lung, gamma=gamma.lung,HR=hr.lung))*weibull.one.year(current.age,lambda=lambda.mort, gamma=gamma.mort,HR=hr.mort)
   }

	data$former <- 1*(data$qtyears>0)
	#Cigarettes per day = CPD
	lung.cpd.beta <- ifelse(data$former==0,0.105300795832441,0.0421816490792434) #This is kindof like log(HR), but only matches
	mort.cpd.beta <- ifelse(data$former==0,0.0497334248502077,0.0146138105962148) # 3 decimal points reliably
	#Representation of CPD in the model
	data$cpd <- ifelse(data$former==1, 
								ifelse(data$cpd<15, data$cpd-10.27171,15-10.27171), #former smokers
								ifelse(data$cpd<15, data$cpd-11.27674,15-11.27674) #current smokers
	)
	
	HR.lung <- exp(data$cpd*lung.cpd.beta)
	HR.mort <- exp(data$cpd*mort.cpd.beta)
	
	ifelse(complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears"))),
	       mapply(hoggart.coef, 
	              current.age=data$age[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))], 
                age=data$smoke.age.start[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))], 
                former=data$former[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))], 
                smkyears = data$smkyears[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))],
                hr.lung = HR.lung[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))], 
	              hr.mort = HR.mort[complete.cases(subset(data,select=c("age","smoke.age.start","former","smkyears")))]),NA)
}



risk.llp <- function(data){
   llp.formula <- function(outcome){      
      # pneu = any diagnosis of pneumonia
      # asb = asbestos exposure binary indicator
      # prior.cancer = any prior cancer
      # fam.cancer.onset = 1 or more first degree ; early onset (3 categories)
      # smkyears.cat = 5 categories      
      formula.string <- "pneu+asb+prior.cancer+fam.cancer.onset+smkyears.cat"
      formula(paste(outcome,"~",formula.string,collapse=""))
   }
   
   llp.coef <- c(intercept =1,
                 pneu=0.602,
                 asb=0.634,
                 prior.cancer=0.675,
                 fam.cancer.onset1=0.703,
                 fam.cancer.onset2=0.168,
                 smkyears.cat1=0.769,
                 smkyears.cat2=1.452,
                 smkyears.cat3=2.507,
                 smkyears.cat4=2.724)
   
   # To make the probabilities interpretable, they used Liverpool data to estimate the beta0, intercept here is from model
   llp.intercepts <- data.frame(
      incidence = c(15.5, 37.87, 88.65, 172.26, 329.02, 487.42, 616.45, 950.61, 1096.42,
                    5.97, 37.34, 68.14,175.24, 230.6, 288.06, 464.99, 594.19 ,497.09)/100000,
      intercept = -c(9.06,8.16, 7.31,6.63,5.97,5.56,5.31,4.83,4.68,
                     9.90,8.06,7.46,6.50,6.22,5.99,5.49,5.23,5.42),
      age = rep(1:9,2),
      female = rep(c(0,1),each=9))
   
   llp.intercept <- function(age, female){     
      #Assign age groups
      age.group <- function(age){
        if(age<45) 1
        else if(age<50) 2
        else if(age<55) 3
        else if(age<60) 4
        else if(age<65) 5
        else if(age<70) 6
        else if(age<75) 7
        else if(age<80) 8
        else 9
      }
     
      value <- function(age, age.cat, female){       
        if(age.cat==1)
          llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age==age.cat]
        else if(age.cat==2)
          sum(c(50-age-0.5, 0.5+age-45)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=2&llp.intercepts$age<=3])/5
        else if(age.cat==3)
          sum(c(55-age-0.5, 0.5+age-50)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=3&llp.intercepts$age<=4])/5
        else if(age.cat==4)
          sum(c(60-age-0.5, 0.5+age-55)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=4&llp.intercepts$age<=5])/5
        else if(age.cat==5)
          sum(c(65-age-0.5, 0.5+age-60)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=5&llp.intercepts$age<=6])/5
        else if(age.cat==6)
          sum(c(70-age-0.5, 0.5+age-65)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=6&llp.intercepts$age<=7])/5										
        else if(age.cat==7)
          sum(c(75-age-0.5, 0.5+age-70)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=7&llp.intercepts$age<=8])/5	
        else if(age.cat==8)
          sum(c(80-age-0.5, 0.5+age-75)*llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age>=8&llp.intercepts$age<=9])/5	
        else
          llp.intercepts$intercept[llp.intercepts$female==female&
            llp.intercepts$age==age.cat]
      }
      
      value(age, age.group(age), female)	
   }
  
   data$fam.lung.cancer <- 1*(data$fam.lung.trend>0)
   vars <- c("female",
   			     "fam.lung.cancer",
             "age",            
             "asb",            
             "smkyears",            
             "pneu",            
             "prior.cancer",
             "fam.cancer.onset")
   
   var.not.in <- names(data)[!(names(data)%in%vars) & sapply(data, function(x) sum(is.na(x)))==0][1] # SELECT ARBITRARY VAR FOR OUTCOME
   expit <- function(x) exp(x)/(1+exp(x))
   data$fam.cancer.onset <- factor(data$fam.lung.cancer,lev=c(0,1,2))  
   data$smkyears.cat <- cut(data$smkyears, c(-1,1,21,41,60,max(data$smkyears,na.rm=TRUE)+1), right=FALSE, lab=1:5)
   #Manually construct the intercept
   intercept <- mapply(llp.intercept, age=data$age, female=data$female)
   X <- model.matrix(llp.formula(var.not.in), data)
   X[,1] <- intercept[complete.cases(subset(data,select= c("fam.lung.cancer",
                                                           "asb",            
                                                           "smkyears",            
                                                           "pneu",            
                                                           "prior.cancer",
                                                           "fam.cancer.onset")))]

   ifelse(complete.cases(subset(data,select=c("fam.lung.cancer",
                                              "asb",            
                                              "smkyears",            
                                              "pneu",            
                                              "prior.cancer",
                                              "fam.cancer.onset"))),expit(X%*%llp.coef),NA)
} 



risk.llpi <- function(data){
   data$male <- as.numeric(data$female==0)
   #family history of lung cancer treated as late onset
   lp <- 0.036*data$age+0.391*data$male+0.043*data$smkyears+0.890*data$copd+1.044*data$prior.cancer+0.521*I(data$fam.cancer.onset==1)+0.071*I(data$fam.cancer.onset==2)

   #0.9728386 is the baseline survival at 8.7 years in the liverpool data
   #3.556 is the mean linear predictor in the liverpool data
   as.numeric(1-0.9728386**(exp(lp-3.556)))
}      



risk.spitz <- function(data){
   spitz.former.formula <- function(outcome){  
     # emp = emphysema (binary indicator)
     # dust = dust exposure  (binary indicator)
     # fam.cancer = 2 or more first degree relatives
     # fam.smoke.cancer = 1 or more first degree
     # no.hayfever = no hay fever (binary indicator)
     # age.stopped = categorical 3 groups	
     
     formula.string <- "emp+dust+fam.cancer+no.hayfever+age.stopped"   
     formula(paste(outcome,"~",formula.string,collapse=""))
   }
   
   spitz.current.formula <- function(outcome){
     # emp = emphysema (binary indicator)
     # dust = dust exposure  (binary indicator)
     # asb = asbestos exposure binary indicator
     # fam.smoke.cancer = 1 or more first degree
     # no.hayfever = no hay fever (binary indicator)
     # packyears.cat = categorical 4 groups	
     
     formula.string <- "emp+dust+asb+fam.smoke.cancer+no.hayfever+packyears.cat"  
     formula(paste(outcome,"~",formula.string,collapse=""))
   }
   
   spitz.baseline <- function(rr, age, female, former){
     age.group <- function(age){
       if(age<45) 1
       else if(age<50) 2
       else if(age<55) 3
       else if(age<60) 4
       else if(age<65) 5
       else if(age<70) 6
       else if(age<75) 7
       else if(age<80) 8
       else if(age<85) 9
       else 10
     }
     
     age.cat <- sapply(age,age.group)
     
     #I is incidence rate, h2 is mortality from something not lung cancer
     I <- spitz.rates$incidence[!is.na(spitz.rates$female)&!is.na(female)&spitz.rates$female==female&!is.na(spitz.rates$age)&!is.na(age.cat)&spitz.rates$age==age.cat] #This selects one gender and one agecat
     h2 <- spitz.rates$mortality[!is.na(spitz.rates$female)&!is.na(female)&spitz.rates$female==female&!is.na(spitz.rates$age)&!is.na(age.cat)&spitz.rates$age==age.cat]
     
     #Generate h1 (hazard from lung cancer) for men/women and former/current smokers
     h1 <- ifelse(is.na(female) | is.na(former),NA,
                  ifelse(female==1&former==1,I*3.76*(1-0.45352),
                         ifelse(female==1&former==0,I*4.17*(1-0.51404),
                                ifelse(female==0&former==1,I*3.17*(1-0.45352),I*3.88*(1-0.51404)))))
     #Model output
     h1*rr/(h1*rr+h2)*(1-exp(-(h1*rr+h2)))
   }
   
   spitz.former.coef <- c(intercept =0,
                          emp=0.9734,
                          dust = 0.4654,
                          fam.cancer = 0.4636,
                          no.hayfever = 0.3711,
                          age1 = 0.2130,
                          age2 = 0.4080)
   
   spitz.current.coef <- c(intercept =0,
                           emp=0.7561,
                           dust = 0.3067,
                           asb=0.4109,
                           fam.smoke.cancer = 0.3859,
                           no.hayfever = 0.4047,
                           py1 = 0.2219,
                           py2 = 0.3747,
                           py3 = 0.6151)
    
   spitz.rates <- data.frame(
     incidence = c(10.78, 25.49, 56.60, 116.58, 221.18, 346.77, 478.10, 564.36, 532.36, 498.44, 
                   11.03, 23.19, 45.51, 93.93, 164.9, 246.85, 318.69, 344.67, 308.28, 266.72)/100000,
     mortality = c(275.1, 400.7, 560.0, 786.9, 1210.2, 1855.1, 2947.4, 4836.4, 7980.7, 15559.4,
                   153.20, 218.80, 313.40, 479.10, 762.90, 1197.00, 1968.30, 3306.10, 5761.20, 14016.2)/100000,
     age = rep(1:10,2),
     female = rep(c(0,1),each=10))
     #Incidence and mortality rates are each done for women and men separately

   data$age.stopped[is.na(data$age.stopped)] <- 0
   data$former <- 1*(data$qtyears>0)   
   vars <- c("female", 
             "age",
             "asb",                  
             "emp",
             "dust",
             "former",
             "no.hayfever",
             "pkyr.cat",
             "fam.cancer",
             "fam.smoke.cancer",
             "age.stopped")
  
   var.not.in <- names(data)[!(names(data)%in%vars) & sapply(data, function(x) sum(is.na(x)))==0][1] # SELECT ARBITRARY VAR FOR OUTCOME

   data$packyears.cat <- cut(data$pkyr.cat,c(-1,28,42,57.5,max(data$pkyr.cat,na.rm=TRUE)+1),right=FALSE,labels=1:4)
   data$age.stopped <- cut(data$age.stopped, c(-1,42,53.1,max(data$age.stopped,na.rm=TRUE)+1),right=FALSE,labels=1:3)

   # DATA FRAME HAS TO HAVE VARIABLES FOR EACH MODEL
   rr.former <- ifelse(complete.cases(subset(data,select=c("emp","dust","fam.cancer","no.hayfever","age.stopped"))),model.matrix(spitz.former.formula(var.not.in), data)%*%spitz.former.coef,NA)
   rr.current <- ifelse(complete.cases(subset(data,select=c("emp","dust","asb","fam.smoke.cancer","no.hayfever","packyears.cat"))),model.matrix(spitz.current.formula(var.not.in), data)%*%spitz.current.coef,NA)
    
   rr.former <- exp(rr.former)
   rr.current <- exp(rr.current)

   risk.former <- ifelse(complete.cases(data$female)&complete.cases(data$former),
                         mapply(spitz.baseline,
                                rr=rr.former[complete.cases(data$female)&complete.cases(data$former)],
                                age=data$age[complete.cases(data$female)&complete.cases(data$former)],
                                female=data$female[complete.cases(data$female)&complete.cases(data$former)],
                                MoreArgs=list(former=1)),NA)
    
   risk.current <- ifelse(complete.cases(data$female)&complete.cases(data$former),
                          mapply(spitz.baseline,
                                rr=rr.current[complete.cases(data$female)&complete.cases(data$former)],
                                age=data$age[complete.cases(data$female)&complete.cases(data$former)],
                                female=data$female[complete.cases(data$female)&complete.cases(data$former)],
                                MoreArgs=list(former=0)),NA)

   ifelse(data$former==1,risk.former,risk.current)
} 



risk.tammemagi <- function(data){
   tammemagi.formula <- function(outcome){   
       # age = continuous
       # race = categorical
       # edu = ordinal
       # bmi = continuous
       # fam.lung.cancer = binary
       # prior.cancer = binary
       # copd = binary
       # cpd = cigarettes per day
       # current = binary
       # smkyears = years smoked
       # qtyears = years quit
       formula.string <- "I(age-62)+
                    black+hispanic+asian+indian+islander+
                    I(edu6-4)+
                    I(bmi-27)+
                    fam.lung.cancer+
                    prior.cancer+
                    copd+
                    current+
                    I(10/cpd-.402154)+
                    I(qtyears-10)+
                    I(smkyears-27)"              
       
       formula(paste(outcome,"~",formula.string,collapse=""))
   }

   tammemagi.coef <- c(intercept=-4.532506,
                       age = 0.0778868,
                       black = 0.3944778,
                       hispanic = -0.7434744,
                       asian = -0.466585,
                       indian = 1.027152,
                       islander = 0,
                       edu = -0.0812744,
                       bmi = -0.0274194,
                       fam.lung.cancer = 0.587185,
                       prior.cancer = 0.4589971,
                       copd = 0.3553063,
                       current = 0.2597431,
                       cpd = -1.822606,
                       qtyears = -0.0308572,
                       smkyears = 0.0317321)
              
   data$black <- 1*(data$race==1)
   data$hispanic <- 1*(data$race==2)
   data$fam.lung.cancer <- 1*(data$fam.lung.trend>0)
   data$current <- 1*(data$qtyears==0)
   vars <- c("age",
             "black","hispanic","asian","islander","indian",
             "edu6",
             "bmi",
             "fam.lung.cancer",
             "prior.cancer",
             "copd",
             "cpd",
             "current",
             "smkyears",
             "qtyears")
  
   var.not.in <- names(data)[!(names(data)%in%vars) & sapply(data, function(x) sum(is.na(x)))==0][1] # SELECT ARBITRARY VAR FOR OUTCOME
 
   X <- model.matrix(tammemagi.formula(var.not.in), data)
   lp <- ifelse(complete.cases(subset(data,select=vars)),X%*%tammemagi.coef,NA)
  
   exp(lp)/(1+exp(lp))
} 



risk.kovalchik <- function (begin, end, newdata, coxph1, ...) 
  {
    c.coxph.risk <- function (coxph, ...) 
    {
      models <- list(coxph)
      if (!missing(..1)) 
        models <- c(models, list(...))
      models
    }  
  
    coxph.relrisk.uncentered <- function (coxph.object, newdata) 
    {
      center <- coxph.object$means %*% coef(coxph.object)
      if (!missing(newdata)) 
        lp <- predict(coxph.object, newdata, type = "lp")
      else lp <- coxph.object$linear.predictor
      exp(lp + rep(center,length(lp)))
    }  

    projection.relrisk <- function (object, data) 
    {
      if (is.numeric(object)) 
        return(object)
      else if (class(object) == "coxph") 
        if (missing(data)) 
          return(coxph.relrisk.uncentered(object))
        else return(coxph.relrisk.uncentered(object, data))
      else stop(cat("No method for class", class(object)))
    }    
    
    #calculate absolute risk given hazards/survival and relative risks
    risk.fixed.interval <- function(H, RR) {
      if (!is.list(H)) 
        0
      else {
        absrisk <- H[[1]]$surv^RR[1] * H[[1]]$haz * RR[1]
        for (i in 2:length(H)) {
          absrisk <- absrisk * (H[[i]]$surv^RR[i])
        }
        sum(absrisk)
      }
    }
    
    models <- c.coxph.risk(coxph1, ...)

    #calculate relative risk for each subject and store as list
    rr <- sapply(models, projection.relrisk, data = newdata)

    if (is.matrix(rr)) { rr.list <- lapply(1:nrow(rr), function(x) rr[x, ])
    } else { rr.list <- list(rr)}

    #estimate risk
    if (length(begin) == 1) {    
       AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
        in.interval <- function(x, begin, end) x >= begin & x <= end
        if ("lung.cancer.death" %in% AllVars){
          H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end),], models[[2]]$basehaz_LCDRAT[in.interval(models[[1]]$basehaz$time, begin, end),])
        } else if ("case" %in% AllVars){
          H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end),], models[[2]]$basehaz_LCRAT[in.interval(models[[1]]$basehaz$time, begin, end),])        
        }
        #calculate absolute risk given hazards/survival for average covariate values and relative risks
        risks <- mapply(risk.fixed.interval, RR = rr.list, MoreArgs = list(H = H))
    } else {
        risks <- mapply(risk, begin = begin, end = end, RR = rr.list, MoreArgs = list(models = models))
    }
    ifelse(complete.cases(newdata),risks,NA)
}

risk.pittsburgh <- function(data){
  S <- -10*I(data$smkyears<30)+8*I(data$smkyears>=40&data$smkyears<50)+14*I(data$smkyears>=50)+
    4*I(data$smkyears<30&data$age>=57)+4*I(data$smkyears>=30&data$smkyears<40&data$age>=59)+
    4*I(data$smkyears>=40&data$smkyears<50&data$age>=61)+4*I(data$smkyears>=50&data$age>=68)-
    3*I(data$qtyears>=1)-4*I(data$cpd<20)+2*I(data$cpd>=30&data$cpd<40)+5*I(data$cpd>=40)  
  
  as.numeric(1/(1+exp(-1*(-4.2195+0.10*S))))
}      


#' Wrapper to the cut() function
#'
#' A wrapper to the cut() function, so that you can automatically break into quantiles as
#' the default behavior, otherwise if the breakpoints are included, then just break on those.
#' In all cases, include.lowest is set to True      
#' @export
cuts <- function(data,npieces,simple.labels=T,...) {

  if (length(npieces)==0 | any(is.na(npieces)) | any(!is.numeric(npieces)))
    stop("npieces must be a numeric scalar or a vector of breakpoints")
  
  if (length(npieces)==1)
    # Just break into quantiles.  Use quantile labelling: Q1,Q2,Q3,etc. if you want
    if (simple.labels)
      cutdata <- cut(data,breaks=quantile(data,seq(0,1,1/npieces),na.rm=T,...),
                     include.lowest=T,labels=paste("Q",1:npieces,sep=""),...)
    else
      cutdata <- cut(data,breaks=quantile(data,seq(0,1,1/npieces),na.rm=T,...),
                     include.lowest=T,...)
  else
    # Break into pieces as specified by the breakpoints
    if (simple.labels)
      cutdata <- cut(data,breaks=npieces,include.lowest=T,
                     labels=paste("Q",1:(length(npieces)-1),sep=""),...)
    else
      cutdata <- cut(data,breaks=npieces,include.lowest=T,...)
  
  return(cutdata)
}



#' Risk Predictions from Lung Cancer Models
#'
#' The R package provides individual risks of lung cancer and lung cancer death 
#' based on various published papers: Bach et al., 2003; Spitz et al., 2007; 
#' Cassidy et al., 2008 (LLP); Hoggart et al., 2012; Tammemagi et al., 2013 (PLCOm2012); 
#' Marcus et al., 2015 (LLPi); Wilson and Weissfeld, 2015 (Pittsburgh);
#' Katki et al., 2016 (LCRAT and LCDRAT).
#'
#' @import coxph.risk
#' 
#' @section Warning:
#'  VGAM is a required dependency of this package.
#'  VGAM may automatically be installed the first time this package is used.
#'  Inputs must be in numerical format to ensure correct output.
#'  For data frame x, this can be checked using sapply(x,class)
#'
#' @param x A data frame or matrix containing individuals' covariate values.  
#'  Covariates should be in the following columns and numerical formats:
#'
#'  \itemize{
#'  \item column 1 - current age (numeric);
#'  \item column 2 - gender (1=Female, 0=Male);
#'  \item column 3 - years smoked (numeric);
#'  \item column 4 - years quit (numeric with 0 to indicate current smoker);
#'  \item column 5 - cigarettes per day (numeric);
#'  \item column 6 - race (0=Non-hispanic white,
#'                   1=Non-hispanic Black/African American, 
#'                   2=Hispanic, 
#'                   3=Other Ethnicity);
#'  \item column 7 - lung disease (1=COPD or Emphysema, 0=No COPD or Emphysema);
#'  \item column 8 - number of first degree relatives with lung cancer (0,1,2);
#'  \item column 9 - bmi;
#'  \item column 10 - highest education level (1=<12 grade, 
#'                                       2=HS graduate, 
#'                                       3=post hs, no college, 
#'                                       4=associate degree/some college, 
#'                                       5=bachelors degree, 
#'                                       6=graduate school);
#'  \item column 11 - asbestos exposure binary indicator;
#'  \item column 12 - prior history of pneumonia indicator;
#'  \item column 13 - prior history of cancer indicator;
#'  \item column 14 - family history of lung cancer (0=none, 1=early onset, 2=late onset);
#'  \item column 15 - Dust exposure  (binary indicator);
#'  \item column 16 - 2 or more first degree relatives with cancer (binary indicator);
#'  \item column 17 - 1 or more first degree relatives with smoking cancer (binary indicator);
#'  \item column 18 - no hay fever (binary indicator);
#'  \item column 19 - asian ethnicity (binary indicator);
#'  \item column 20 - islander ethnicity (binary indicator);
#'  \item column 21 - American indian ethnicity (binary indicator);
#'  }
#'  
#'  Risk factors used by each model
#'  \tabular{cccccccccc}{
#'  Covariate \tab Bach \tab Spitz \tab LLP \tab Hoggart \tab PLCOm2012 \tab Pittsburgh \tab LLPi \tab LCRAT \tab LCDRAT\cr
#'  Age \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \cr
#'  Gender \tab Yes \tab Yes \tab Yes \tab No \tab No \tab No \tab Yes \tab Yes \tab Yes \cr
#'  Race/ethnicity \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab Yes \tab Yes\cr
#'  Asian \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab No \tab No\cr
#'  Pacific islander \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab No \tab No\cr
#'  America Indian \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab No \tab No\cr
#'  Education \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab Yes \tab Yes\cr 
#'  BMI \tab No \tab No \tab No \tab No \tab Yes \tab No \tab No \tab Yes \tab Yes\cr
#'  Smoking status \tab No \tab Yes \tab No \tab Yes \tab Yes \tab Yes \tab No \tab No \tab No\cr
#'  Years/Age quit \tab Yes \tab Yes \tab No \tab Yes \tab Yes \tab No \tab No \tab Yes \tab Yes\cr
#'  Years smoked \tab Yes \tab No \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes \tab Yes\cr
#'  Cigs per day \tab Yes \tab No \tab No \tab Yes \tab Yes \tab Yes \tab No \tab Yes \tab Yes\cr
#'  Pack-years \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab Yes \tab Yes\cr
#'  Prior cancer \tab No \tab No \tab Yes \tab No \tab Yes \tab No \tab Yes \tab No \tab No\cr
#'  Lung disease \tab No \tab Yes \tab No \tab No \tab Yes \tab No \tab Yes \tab Yes \tab Yes\cr
#'  Pneumonia \tab No \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Hayfever \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Asbestos exposure \tab Yes \tab Yes \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Dust exposure \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Any FDR w/ cancer \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Any FDR w/ smoking-related cancer \tab No \tab Yes \tab No \tab No \tab No \tab No \tab No \tab No \tab No\cr
#'  Any FDR w/ LC \tab No \tab No \tab Yes \tab No \tab Yes \tab No \tab Yes \tab Yes \tab Yes\cr
#'  Num. FDR w/ LC \tab No \tab No \tab No \tab No \tab No \tab No \tab No \tab Yes \tab Yes\cr
#'  Early onset FDR w/ LC \tab No \tab No \tab Yes \tab No \tab No \tab No \tab Yes \tab No \tab No
#'  }
#'
#' @return A numeric matrix containing individuals' predictions:
#' 
#'  \itemize{
#'  \item column 1 - An indicator variable for whether the individual is eligible 
#'              for CT lung screening according to 
#'              US Preventive Services Task Force (USPSTF) recommendations.
#'  \item column 2 - This is the probability of dying from lung cancer 
#'              within 5 years if not undergoing screening (Katki, 2016).
#'  \item column 3 - This is the reduction in the probability of dying from lung cancer in 5 years 
#;              if undergoing 3 yearly CT lung screens as in the NLST (Katki, 2016).
#'  \item column 4 - This is the probability of being diagnosed with lung cancer
#'              within 5 years if not undergoing screening (Katki, 2016).
#'  \item column 5 - This is the extra probability of lung cancer diagnosis in 5 years
#'              if undergoing 3 yearly CT lung screens as in the NLST  (Katki, 2016).
#'  \item column 6 - This is the probability of having at least one false-positive 
#'              CT screen out of 3 screens (Katki, 2016).
#'  \item column 7 - This is the expected number of false-positive CT screens 
#'               after 3 screens (Katki, 2016).
#'  \item column This is the probability of being diagnosed with lung cancer 
#'              within 10 years if not undergoing screening (Bach, 2003).
#'  \item column 9 - This is the probability of being diagnosed with lung cancer 
#'              within 1 years if not undergoing screening (Hoggart, 2012).
#'  \item column 10 - This is the probability of being diagnosed with lung cancer 
#'              within 5 years if not undergoing screening (LLP, 2008).
#'  \item column 11 - This is the probability of being diagnosed with lung cancer 
#'              within 8.7 years if not undergoing screening (LLPi, 2015).
#'  \item column 12 - This is the probability of being diagnosed with lung cancer 
#'              within 1 years if not undergoing screening (Spitz, 2007).
#'  \item column 13 - This is the probability of being diagnosed with lung cancer 
#'              within 6 years if not undergoing screening (Tammemagi, 2013).
#'  \item column 14 - This is the probability of being diagnosed with lung cancer
#'              within 6 years if not undergoing screening (Pittsburgh, 2015).            
#'  }
#'  
#' @author Li C. Cheung, \email{li.cheung@nih.gov}, Stephanie A. Kovalchik, Hormuzd A. Katki
#'
#' @references 
#' \itemize{
#'  \item Bach PB, Kattan MW, Thornquist MD, et al. 
#'        Variations in lung cancer risk among smokers. J Natl Cancer Inst 2003;95:470-8.
#'  \item Spitz MR, Hong WK, Amos CI, et al. 
#'        A risk model for prediction of lung cancer. J Natl Cancer Inst 2007;99:715-26.
#'  \item Cassidy A, Myles JP, van Tongeren M, et al. 
#'        The LLP risk model: an individual risk prediction model for lung cancer. 
#'        Br J Cancer 2008;98:270-6.
#'  \item Hoggart C, Brennan P, Tjonneland A, et al. 
#'        A risk model for lung cancer incidence. Cancer Prev Res (Phila) 2012;5:834-46.
#'  \item Tammemagi MC, Katki HA, Hocking WG, et al. 
#'        Selection criteria for lung-cancer screening. N Engl J Med 2013;368:728-36.
#'  \item Marcus MW, Chen Y, Raji OY, Duffy SW, Field JK.
#'        LLPi: Liverpool lung project risk prediction model for lung cancer incidence.
#'        Cancer Prev Res (Phila) 2015;8:570-5.
#'  \item Wilson DO, Weissfeld J.
#'        A simple model for predicting lung cancer occurrence in a lung cancer screening program:
#'        the Pittsburgh predictor. Lung Cancer 2015;89:31-37.      
#'  \item Katki HA, Kovalchik SA, Berg CD, Cheung LC, Chaturvedi AK.
#'        Development and validation of risk models to select ever-smokers
#'        for CT lung cancer screening. JAMA. 2016;315:2300-11.
#'        doi: 10.1001/jama.2016.6255.
#'  \item Katki HA, Kovalchik SA, Petito LC, Cheung LC, Jacobs E, Jemal A, Berg CD, Chaturvedi AK.
#'        Implications of nine risk prediction models for selecting ever-smokers for computed
#'        tomography lung cancer screening.  Ann Intern Med. 2018;doi: 10.7326/M17-2701. 
#'           
#' }
#'             
#'
#' @export
#' @examples
#' age <- c(66,58,75,72,56)
#' bmi <- c(23,28,26,27,24)
#' cpd <- c(36,36,40,24,40)
#' emp <- c(0,1,1,0,1)
#' fam.lung.trend <- c(0,2,0,2,0)
#' female <- c(0,1,0,1,0)
#' smkyears <- c(43,37,45,42,29)
#' qtyears <- c(0,0,9,6,6)
#' race <- c(0,1,2,2,3)
#' edu6 <- c(3,5,4,5,5)
#' asb <- c(0,0,0,0,0)
#' pneu <- c(0,0,0,0,0)
#' prior.cancer <- c(0,0,0,0,0)
#' fam.cancer.onset <- c(0,1,0,2,0)
#' dust <- c(0,0,0,0,0)
#' fam.cancer <- c(0,1,0,1,0)
#' fam.smoke.cancer <- c(0,1,0,1,0)
#' no.hayfever <- c(1,1,1,1,1)
#' asian <- c(0,0,0,0,1)
#' islander <- c(0,0,0,0,0)
#' indian <- c(0,0,0,0,0)
#' 
#' persons <- data.frame(age,
#'                       female,
#'                       smkyears,
#'                       qtyears,
#'                       cpd,
#'                       race,
#'                       emp,
#'                       fam.lung.trend,
#'                       bmi,
#'                       edu6,
#'                       asb,
#'                       pneu,
#'                       prior.cancer,
#'                       fam.cancer.onset,
#'                       dust,
#'                       fam.cancer,
#'                       fam.smoke.cancer,
#'                       no.hayfever,
#'                       asian,
#'                       islander,
#'                       indian)
#' 
#' persons_predictions <- lcmodels(persons)
#' persons_predictions

lcmodels <- function(x) {
  if (!require("VGAM")) {
    install.packages("VGAM",repos="http://watson.nci.nih.gov/cran_mirror/")
  } 
  library("VGAM")
  library("coxph.risk")   
      
  x <- data.matrix(x)
  
  age <- x[,1]
  female <- x[,2]
  smkyears <- x[,3]
  qtyears <- x[,4]
  cpd <- x[,5]
  race <- as.factor(x[,6])
  emp <- x[,7]
  fam.lung.trend <- x[,8]
  bmi <- x[,9]
  edu6 <- x[,10]
  pkyr.cat <- smkyears*cpd/20
  age.stopped <- age-qtyears
  smoke.age.start <- age-smkyears-qtyears
  asb <- x[,11]
  pneu <- x[,12]
  prior.cancer <- x[,13]
  fam.cancer.onset <- as.factor(x[,14])
  copd <- emp
  dust <- x[,15]
  fam.cancer <- x[,16]
  fam.smoke.cancer <- x[,17]
  no.hayfever <- x[,18]
  asian <- x[,19]
  islander <- x[,20]
  indian <- x[,21]
   
   covar <- data.frame(age=age,
                       bmi=bmi,
                       cpd=cpd,
                       emp=emp,
                       fam.lung.trend=fam.lung.trend,
                       female=female,
                       qtyears=qtyears,
                       smkyears=smkyears,
                       race=race,
                       edu6=edu6,
                       pkyr.cat=pkyr.cat,
                       asb=asb,
                       smoke.age.start=smoke.age.start,
                       pneu = pneu,
                       prior.cancer = prior.cancer,
                       fam.cancer.onset = fam.cancer.onset,
                       copd=copd,
                       dust=dust,
                       fam.cancer=fam.cancer,
                       fam.smoke.cancer=fam.smoke.cancer,
                       no.hayfever=no.hayfever,
                       age.stopped=age.stopped,
                       asian = asian,
                       islander = islander,
                       indian = indian)

   ### Bach 2003 ###   
   r.bach <- risk.bach.10(covar)
   ### Hoggart 2012 ###
   r.hoggart <- risk.hoggart(covar)
   ### LLP 2008 ###
   r.llp <- risk.llp(covar)
   ### LLPi 2015 ###
   r.llpi <- risk.llpi(covar)
   ### Spitz 2007 ###
   r.spitz <- risk.spitz(covar)
   ### Tammemagi 2013 ###
   r.tamme <- risk.tammemagi(covar)
   ### Pittsburgh 2015 ###
   r.pitt <- risk.pittsburgh(covar)
      
   uspstf_elig <- ifelse(age>=55 & age<=80 & pkyr.cat>=30 & qtyears <= 15,1,0)
   
   #calculate predicted risks of lung cancer death within 5 years based on LCDRAT with competing risk of death
   covar$LCDRAT <- risk.kovalchik(0, 5, covar[,1:11], LCDRAT, cox.death)
   
   #calculate predicted risks of lung cancer death using chest xray
   covar$cxLCDRAT <- 0.796*covar$LCDRAT
   
   #calculate predicted risks of lung cancer incidence within 5 years based on LCRAT with competing risk of death
   covar$LCRAT <- risk.kovalchik(0, 5, covar[,1:11], LCRAT, cox.death)
   covar$LCRAT <- pmax(covar$LCRAT,covar$LCDRAT)
   
   #calculate predicted risks of lung cancer incidence using chest xray
   covar$cxLCRAT <- 1.124*covar$LCRAT
   
   #calculate probability of Y, number of false positives, taking values of 0, 1, 2, or 3 
   #responses are reported as log(P(Y=1)/P(Y=0)), log(P(Y=2)/P(Y=0)), and log(P(Y=3)/P(Y=0))
   prob_numfalsepos <- predict(polytmod,type="response",newdata=covar)
   
   #solving for P(Y=0), P(Y=1), P(Y=2), and P(Y=3), we have:
   covar$prob_0falsepos <- ifelse(is.na(covar$cxLCRAT),NA,prob_numfalsepos[,1])
   covar$prob_1falsepos <- ifelse(is.na(covar$cxLCRAT),NA,prob_numfalsepos[,2])
   covar$prob_2falsepos <- ifelse(is.na(covar$cxLCRAT),NA,prob_numfalsepos[,3])
   covar$prob_3falsepos <- ifelse(is.na(covar$cxLCRAT),NA,prob_numfalsepos[,4])
   covar$expected_falsepos <- covar$prob_1falsepos + 2*covar$prob_2falsepos + 3*covar$prob_3falsepos
   
   ### US preventive service task force eligible - uspstf_elig
   
   ### probability of dying from lung cancer death in 5 years without screening - covar$LCDRAT  
   
   ### Reduction in 5 year lung cancer death risk due to screening - (covar$LCDRAT - covar$cxLCDRAT)
   
   ### chance of lung cancer diagnosis without screening - covar$LCRAT
   
   ### increase in lung cancer diagnosis due to screening - (covar$cxLCRAT - covar$LCRAT)
   
   ### chance of false positive CT lung screen - (1-covar$prob_0falsepos)
   
   predicted <- data.frame(uspstf_elig,
                           covar$LCDRAT,
                           (covar$LCDRAT - covar$cxLCDRAT),
                           covar$LCRAT,
                           (covar$cxLCRAT - covar$LCRAT),
                           (1-covar$prob_0falsepos),
                           covar$expected_falsepos,
                           r.bach,
                           r.hoggart,
                           r.llp,
                           r.llpi,
                           r.spitz,
                           r.tamme,
                           r.pitt)
                           
   colnames(predicted) <- c("USPSTF eligible",
                            "Risk of lung cancer death within 5 years without screening  (Katki, 2016)",
                            "Reduction in 5-years risk of lung cancer death with CT screening  (Katki, 2016)",
                            "Risk of lung cancer diagnosis within 5 years without screening  (Katki, 2016)",
                            "Increase in 5-years risk lung cancer diagnosis with CT screening (Katki, 2016)",
                            "Probability of a false-positive CT lung screens (Katki, 2016)",
                            "Expected number of false-positive CT lung screens (Katki, 2016)",
                            "Risk of lung cancer diagnosis within 10 years (Bach, 2003)",
                            "Risk of lung cancer diagnosis within 1 year (Hoggart, 2012)",
                            "Risk of lung cancer diagnosis within 5 years (LLP, 2008)",
                            "Risk of lung cancer diagnosis within 8.7 years (LLPi, 2015)",
                            "Risk of lung cancer diagnosis within 1 year (Spitz, 2007)",
                            "Risk of lung cancer diagnosis within 6 years (Tammemagi, 2013)",
                            "Risk of lung cancer diagnosis within 6 years (Pittsburgh, 2015)")
   predicted                         
}
