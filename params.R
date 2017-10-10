#---------------------------------------------------------------------------------------
# FUNCTION TO GET Parameters of model
#---------------------------------------------------------------------------------------

set_param<-function(){
  
  param=list()
  
  #---------------------------------------------------------------------------------------
  # Parameters that can vary
  #---------------------------------------------------------------------------------------
  
  #refuges <- read.delim("data/FORCE_refuge_ss_per_site.txt", header = TRUE)
  #param$refuge <- refuges[,8]
  
  #Algal dynamics
  param$alr = 109#0.2 #Maximum algal growth rate per year
  param$AK = 6 #Carrying capacity of algae
  param$min.A=0#0.1               #proportion of population with reduced vulnerability (range: 1-0.1, 1 = no complexity)
  
  #Primary productivity and search rates
  param$pp=-0.5 #1.07            #Primary productivity 
  #param$A.u = 10
  param$A.u_reef= 6.4 #10            #Yearly rate volume searched by predators pre settlement
  param$A.u_pel = 10            #Yearly rate volume searched by predators post settlement
  param$A.u = 6.4
  param$A.v= 0.1 #param$A.u_reef / 10               #Yearly rate volume searched by invertebrates (should be 10x less than predators?)
  param$A.h= 0.2 #1.25                   #Yearly rate volume searched by herbivores
  param$flow = 10
  
  #Theoretical complexity
  param$setsp=1000        #slope of decrease (settlement onto reef structure)
  param$emsp=100          #increasing slope (emergence from reef structure)
  param$settle=0.1        #weight when settle to reef structure
  param$emerge=1000       #maximum refuge size in body mass (range explored = 1-1500) 
  param$min.A=0
  
  #Feeding preferences
  param$pref.pel=0.33           #Preference for "pelagic" carnivorous prey 
  param$pref.ben=0.5            #Preference for "benthic" invertebrate prey
  param$pref.herb=0.17          #Preference for herbivorous prey

  param$pref.det = 0
  param$pref.alg = 1
  
  #Fishing
  param$Fmort_pred=0            #Level of fishing mortality rate yr^(-1)
  param$Fmort_herb=0
  param$min.fishing.size=1      #Minimum log10 body size fished
  param$max.fishing.size=3.5    #Maximum log10 body size fished
    
  #---------------------------------------------------------------------------------------
  # Fixed parameters
  #---------------------------------------------------------------------------------------
  
  param$q0=2.0                  #Mean log10 predator prey mass ratio  100:1.(beta in paper)
  param$sd.q=1.0                #0.6 to 1.0 for lognormal prey preference function. 5=B&R ref value (sigma, standard deviation)
  param$qmax=param$q0 + 2*param$sd.q        #Truncation of the feeding kernel 95% distribution between +/- 2*sd.q
  param$qmin=param$q0 - 2*param$sd.q
  
  param$det.coupling=1.0
  param$sinking.rate=0.8        # 35% of PP material reaches the seafloor davies & Payne 1983 Marine biology
  param$alpha=0.75#0.82              # exponent for metabolic requirements plus swimming for predators(Ware et al 1978)
                                # NOTE:this exponent =0.75 for whole organism basal (sedentary) metabolic rate (see growth.v) from Peters (1983) and Brown et al. (2004) for detritivores
  param$alpha.h=0.75#0.6046
  
  param$K.u=0.15                 #Gross growth conversion efficiency for organisms in the "predator" spectrum Ware (1978)
  param$K.v=0.15                 #Gross growth conversion efficiency for organisms in the "detritivore" spectrum
  param$K.h=0.15                 #Gross growth conversion efficiency for organisms in the "herbivore" spectrum
  param$K.d=0.1                 #Gross growth conversion efficiency for detritus
  param$def.high=0.4            #Fraction of ingested food that is defecated (Peters,1983)
  param$def.low=0.4             #Low = low quality (K) food, high = high quality (K) food
  
  param$K.a=0.1
  
  param$mu0=0.2	                #Residual natural mortality
  param$k.sm =0.1               #Constant for senescence mortality 
  param$xs=3                    #Size at sensenscence e
  param$p.s=0.3  			          #Exponent of senescence
  
  #Initial size spectra properties
  #param$ui0=10^(param$pp)      #Initial intercept of plankton size spectrum, from P.McCloghrie ERSEM output.
  #param$vi0=sinking.rate*ui0   #Initial intercept of detritivore size spectrum  - same fraction as sinking rate
  param$r.plank=-1.0            #Initial (fixed) slope of phyto-zooplankton spectrum
  param$pred.slope=-1#-0.75           #Initial slope of fish spectrum, same as that used in B&R
  param$herb.slope= -1#-0.64          #Initial slope of detrtivore spectrum
  param$invert.slope=-0.36
                                  
                                #Proportion of initial plankton intercept that depicts the abundance of larvae / eggs of 
  param$pred.prod = 1#(10^0.56)
  param$herb.prod= 1#(10^-0.04) #1.65           #Herbivorous fish
  param$invert.prod = 10^2.1#(10^-0.49) # 5.4
    #param$sinking.rate         #Invertebrates
  
  param$W.init=5                #Initial detritus density per scenario
  param$A.init=25
  
  #-------------------------------------------
  ##New parameters for reproductive investment
  
  param$r.u=0.1
  param$r.h=0.1
  param$r.v=0.1
  param$r.a=0.1
  param$r.d=0.1
  
  
  #---------------------------------------------------------------------------------------
  # Parameters for numerical integration
  #---------------------------------------------------------------------------------------
  param$dx=0.1                  #Size increment after discretization for integration (log body weight) 
  param$xmin=-12                #Minimum log10 body size of plankton
  param$x1=-1.5#-3.5                   #Minimum log10 body size in predators
  param$x1.det=-4               #Minimum log10 body size in dynamics benthic detritivores
  param$x1.herb= -1.5 #-3.5           #Minimum log10 body size of herbivores
  param$xmax=3.5                #Maximum log10 body size of predators
  param$xmax2=3                 #Maximum log10 body size before senescence kicks in (departure form linearity)
  
##******NEW
  param$inv_max <- 3.5
  param$invmin = -6
  
  
  ## Vector with size bins 
  param$x    = seq(param$xmin, param$xmax, param$dx)
  param$y    = param$x
  param$end  = length(param$x)

  param$refinv <- which(param$x == param$invmin)
  
#******NEW

  param$inv_end <- which(param$x == param$inv_max)

  param$tstepdays=1.0                             #Timestep in days
  param$tmaxyears=200.0                           #Maximum number of years for model run
  
  param$N=as.integer(365.0*param$tmaxyears/param$tstepdays)	  #Total number of steps in a period
  param$dt=param$tmaxyears/param$N                            #Time step (rates per year)   
  
  param$ref=((param$x1-param$xmin)/param$dx)+1                      #Location of smallest size consumer
  #param$phyto.ref=((-6-param$xmin)/param$dx)+1                #Position in x of upper size of phytoplankton
  param$ref.det=((param$x1.det-param$xmin)/param$dx)+1              #Position in x of x1.det
  param$ref.herb=((param$x1.herb-param$xmin)/param$dx)+1            #Position in x of x1.herb
    
  param$ref2=((param$xmax2-param$xmin)/param$dx)+1                  #Position in x of xmax2 (largest benthic detritivores)
  
  #Setting data ranges for field measurement
  
  #Considering fish bigger than 5cm - or log10 weight > 0
  param$dat_start_5=0
  #Considering fish bigger than 10cm - or log10 weight > 1
  param$dat_start_10=1
  param$dat_end = 3.5#4.5
  
  param$ref.dst_5<-((param$dat_start_5-param$xmin)/param$dx)+1      # position in x of fish bigger than ~5cm (surveyable)
  param$ref.dst_10<-((param$dat_start_5-param$xmin)/param$dx)+1      # position in x of fish bigger than ~10cm (surveyable)
  
  param$ref.den<-((param$dat_end-param$xmin)/param$dx)+1                          # position in x of largest surveyed fish
  
  param$Fref=((param$min.fishing.size-param$xmin)/param$dx)+1   # position in F vector corresponding to smallest size fished
  param$Fref2=((param$max.fishing.size-param$xmin)/param$dx)+1
    
  return(param)
}
