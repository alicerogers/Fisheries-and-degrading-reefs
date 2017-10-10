#Source necessary function and parameter files
source('sizemodel_functions.r')
source('params.r')

#Set model parameters using parameter input file
params<-set_param()

#Load data on refuge availability:
#First with branching corals
refuges <- read.delim("data/high_coral_cover_refuges.txt", header = TRUE)
#...and then without
refuges_nb <- read.delim("data/high_coral_cover_refuges_no_branches.txt", header = TRUE)

#Load data on invertebrate 
inverts <- read.delim("data/Kramer_invert_data.txt", header = TRUE)

#==============================
#Initial runs and model checks
#==============================

#Generate model output for default scenario without refuges  
#Set refuge density to zero in all size classes
  params$refuge <- rep(0, length(refuges[,5]))
#Run model and assign to "initial_res"
  initial.res <- try(run_model(params, initial.run = T))

#Generate model output for a test location with refuges to ensure refuge function is working  
#Set refuges as per test location "Pasi_Gusung_North_nonreserve2" 
  params$refuge <- refuges[,5]
#Run_model and ssign to "complex_check"
  complex_check <- run_model(params, initial.run = T)

#------------
#Plot results
#------------

#Plot non-complex size spectra  
  plotsizespectrum_full(initial.res, params)
#Add invertebrate data to check for fit
  points(log_abundance_m2~Mean.log10.weight, data = inverts, pch = 19, col = "brown")
  points(high_log_abundance_m2~Mean.log10.weight, data = inverts, pch = 19, col = "purple")

#Plot complex size spectra  
  plotsizespectrum_full(complex_check, params)
#Add invertebrate data to check for fit
  points(log_abundance_m2~Mean.log10.weight, data = inverts, pch = 19, col = "brown")
  points(high_log_abundance_m2~Mean.log10.weight, data = inverts, pch = 19, col = "purple")
  
#Biomass through time
plotbiomass(initial.res, params)
plotbiomass(complex_check, params)

#Algae and herbivores through time
time <- seq(by = 1, length = params$N)
plot(initial.res$A_time~time, type = "l", ylim = c(0, 110), col = "seagreen", lwd = 2, main = "Algae and herbivores")
points(initial.res$h_bio_time~time, type = "l", col = "green")

plot(complex_check$A_time~time, type = "l", ylim = c(0, 100), col = "seagreen", lwd = 2, main = "Algae and herbivores")
points(complex_check$h_bio_time~time, type = "l", col = "green")

#Detritus and invertebrates through time
plot(initial.res$Det_time~time, type = "l", ylim = c(0, 300), col = "grey", lwd = 2, main = "Detritus and invertivores")
points(initial.res$i_bio_time~time, type = "l", col = "brown")

plot(complex_check$Det_time~time, type = "l", ylim = c(0, 300), col = "grey", lwd = 2, main = "Detritus and invertivores")
points(complex_check$i_bio_time~time, type = "l", col = "brown")

#Growth rates by size
compare_growth(initial.res, params)
compare_growth(complex_check, params)


#Save size spectra data at end of run
saveendspectra(initial.res)
#Can then be read in as initial values for next model run 

#===========================================================================================================================================
#Run model scenarios with changes in algal productivity and invertebrates
#===========================================================================================================================================

#Set refuge list - with branches
refuge <-list(refuges[,2],refuges[,3],refuges[,4],refuges[,5],refuges[,6],refuges[,7],
              refuges[,8],refuges[,9],refuges[,10])

alr <- c(109, 163.5) #Normal versus 50% increase in algal productivity
invert.prod <- c(10^2.1, 10^2.4) #Normal versus increase in invertebrate productivity

simset1<-data.frame(matrix(0,length(refuge)*length(alr)*length(invert.prod),3))
names(simset1)<-c("refuge", "alr", "invert.prod")

#Set combinations of parameters 
simset1$refuge <- rep(refuge, length(alr))
simset1$alr <- rep(alr, each = length(refuge))
simset1$invert.prod <- rep(invert.prod, each = length(refuge)*length(alr))

#Remove combinations that are not relevant to modelled reef scenarios - i.e. high algae, low invert and low algae, high invert
useless <- 10:27
simset1 <- simset1[-useless, ]

#Set refuge list without branches
refuge2 <-list(refuges_nb[,2],refuges_nb[,3],refuges_nb[,4],refuges_nb[,5],refuges_nb[,6],refuges_nb[,7],
              refuges_nb[,8],refuges_nb[,9],refuges_nb[,10])

#Note scenarios without branches always have high algae and inverts - i.e. corals are dead here
alr2 <- 163.5
invert.prod2 <- 10^2.4

#Set new combinations of parameters
simset2<-data.frame(matrix(0,length(refuge2)*length(alr2)*length(invert.prod2),3))
names(simset2)<-c("refuge", "alr", "invert.prod")

simset2$refuge <- rep(refuge2, length(alr2))
simset2$alr <- rep(alr2, each = length(refuge2))
simset2$invert.prod <- rep(invert.prod2, each = length(refuge2)*length(alr2))

#Set refuge list without any structure
refuge3 <-list(rep(0, length(refuges[,2])),rep(0, length(refuges[,2])),rep(0, length(refuges[,2])),rep(0, length(refuges[,2])),
               rep(0, length(refuges[,2])),rep(0, length(refuges[,2])),
               rep(0, length(refuges[,2])),rep(0, length(refuges[,2])),rep(0, length(refuges[,2])))

#Note scenarios without structure always have high algae and inverts - i.e. corals are dead here
alr2 <- 163.5
invert.prod2 <- 10^2.4

#Set new combinations of parameters
simset3<-data.frame(matrix(0,length(refuge3)*length(alr2)*length(invert.prod2),3))
names(simset3)<-c("refuge", "alr", "invert.prod")

simset3$refuge <- rep(refuge3, length(alr2))
simset3$alr <- rep(alr2, each = length(refuge3))
simset3$invert.prod <- rep(invert.prod2, each = length(refuge3)*length(alr2))


#Combine parameter tables for all relevant reef degradation scenarios
simset <- rbind(simset1, simset2, simset3)

#Check number of scenarios
length(simset[,1])

#Option to allow for fishing mortality or not
#Fmort <- c(0,1)
#simset$Fmort <- rep(Fmort, each = length(simset[,1]))

#Create output table to hold results
output <- data.frame(matrix(0, length(simset[,1]), 15))
names(output) <- c("Site_name", "alr","invert.prod", "Mod_pred_biomass","Mod_herb_biomass", "Mod_inv_biomass",
                   "Mod_pred_biomass_5", "Mod_herb_biomass_5", "Pred_prod", "Herb_prod", 
                   "Inv_prod", "Pred_Fprod", "Herb_Fprod", "net_algae", "net_detritus")

#Set time for plots
time <- seq(by = 1, length = params$N)

#Assign count to label plots
count <- 1:length(simset[,1])

#Set site name labels in output
output$Site_name <- rep(names(refuges[2:10]), 6)

#Create lists to hold size spectra data
Pred_ss_dat <- list(1)
Herb_ss_dat <- list(1)
Inv_ss_dat <- list(1)

#Create lists to hold growth data
Pred_growth_data <- list(1)
Herb_growth_data <- list(1)
Inv_growth_data <- list(1)

#Create lists to hold mortality data
Pred_mort_data <- list(1)
Herb_mort_data <- list(1)
Inv_mort_data <- list(1)

#Create list to hold refuge data
Refuge_data <- list(1)

#Create lists to hold catch data
#Pred_yield <- list(1)
#Herb_yield <- list(1)

#Open pdf to capture model output plots- optional
#pdf("results/ss_and_growth_plots.pdf")

#Run model for each scenario
for (i in 1: (length(simset[,1]))){ 
  #i = 1 
  params$refuge <- simset$refuge[[i]]
  params$alr <- simset$alr[i]
  params$invert.prod <- simset$invert.prod[i]
  #params$Fmort_pred <- simset$Fmort[i]
  #params$Fmort_herb <- params$Fmort_pred
  
  #Run model using initial results (without complexity) as starting point
  res <- run_model(params = params, initial.run = F)
  
  #Info data
  output$alr[i] <- res$params$alr
  output$invert.prod[i] <- res$params$invert.prod
  #output$Fmort[i] <- res$params$Fmort_pred
  
  #Total biomass data
  output$Mod_pred_biomass[i] <- res$Pred_gm
  output$Mod_herb_biomass[i] <- res$Herb_gm
  output$Mod_inv_biomass[i]  <- res$Inv_gm
  
  #Biomass data for fish greater than 5cm - or .... log10
  output$Mod_pred_biomass_5[i] <- res$Pred_dat_gm_5
  output$Mod_herb_biomass_5[i] <- res$Herb_dat_gm_5
  
  #Productivity
  output$Pred_prod[i] <- res$Pred_prod
  output$Herb_prod[i] <- res$Herb_prod
  output$Inv_prod[i]  <- res$Inv_prod
  
  #Fisheries sized productivity
  output$Pred_Fprod[i] <- res$Fpred_prod
  output$Herb_Fprod[i] <- res$Fherb_prod
  
  #Time averaged productivity
  #output$Pred_Fprod_time[i] <- mean(res$FPred_prod_time[36500:73000])
  #output$Herb_Fprod_time[i] <- mean(res$FHerb_prod_time[36500:73000])
  #output$Inv_prod_time[i] <- mean(res$Inv_prod_time[36500:73000])
  
  #Algae and detritus
  output$net_algae[i] <- res$A
  output$net_detritus[i] <- res$W
  
  #Fisheries yields
  #output$Tot_annual_pred_yield_gm2[i] <- sum(res$Yield_p) 
  #output$Tot_annual_herb_yield_gm2[i] <- sum(res$Yield_h)
  
#Predator and Herbivore size spectra, growth and mortality for further calculations  
  Pred_ss_dat[[i]] <- res$Preds
  Herb_ss_dat[[i]] <- res$Herbs
  Inv_ss_dat[[i]] <- res$Invs
  
  Pred_growth_data[[i]] <- res$P_grth
  Herb_growth_data[[i]] <- res$H_grth
  Inv_growth_data[[i]] <- res$I_grth
  
  Pred_mort_data[[i]] <- res$P_mrt
  Herb_mort_data[[i]] <- res$H_mrt
  Inv_mort_data[[i]] <- res$I_mrt
  
  Refuge_data[[i]] <- res$Pred_ref
  
  #Pred_yield[[i]] <- res$Yield_p
  #Herb_yield[[i]] <- res$Yield_h

#Plot outputs for checking
  
  #plotsizespectrum_full(res, params)
  #text(2,15, count[i])
  
  #plotbiomass(res, params)
  #text(2,15, count[i])
  
  #compare_growth(res, params)
  #text(-2,4, count[i])
  
  #plot(Pred_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,2), main = "Predator growth"
  #     , col = "red")
  #plot(Herb_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,2), main = "Herbivore growth"
  #     , col = "green")
  #plot(Inv_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,3), main = "Invertebrate growth"
  #     , col = "brown")
  
  #plot(Pred_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,5), main = "Predator mortality"
  #     , col = "red")
  #plot(Herb_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,5), main = "Herbivore mortality"
  #     , col = "green")
  #plot(Inv_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,20), main = "Invertebrate 
  #     mortality", col = "brown")
  
  #plot(Refuge_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,1), main = "Refuge"
  #     , col = "black")
  #text(4,0.8, count[i])
  
  #plot(res$A_time~time, type = "l", col = "seagreen", ylim = c(0, 250))
  #points(res$h_bio_time~time, type = "l", col = "green")
  #text(30000,100, count[i])
}

#Check results
output

#Close pdf
#dev.off()

#Write results table to file
write.table(output, "results/summary.txt", row.names = FALSE)


#Unlist and re-structure size spectra data for further analyses
site <- unique(output$Site_name)
size <- params$x

location_info <- list(1)

for (i in 1: length(site)){
  location_info[[i]] <- rep(as.character(site[i]), length(params$x))
}

location_info <- unlist(location_info[-17])

#Abundance data
pred_data <- as.data.frame(cbind(location_info, unlist(Pred_ss_dat)))
names(pred_data) <- c("Site", "ss_data")
pred_data$size <- rep(size, length(site))
write.table(pred_data, paste("results/preds.txt"))

herb_data <- as.data.frame(cbind(location_info, unlist(Herb_ss_dat)))
names(herb_data) <- c("Site", "ss_data")
herb_data$size <- rep(size, length(site))
write.table(herb_data, paste("results/herbs.txt"))

inv_data <- as.data.frame(cbind(location_info, unlist(Inv_ss_dat)))
names(inv_data) <- c("Site", "ss_data")
inv_data$size <- rep(size, length(site))
write.table(inv_data, paste("results/invs.txt"))

#Growth data
G_pred_data <- as.data.frame(cbind(location_info, unlist(Pred_growth_data)))
names(G_pred_data) <- c("Site", "ss_data")
G_pred_data$size <- rep(size, length(site))
write.table(G_pred_data, paste("results/GRWTH_preds.txt"))

G_herb_data <- as.data.frame(cbind(location_info, unlist(Herb_growth_data)))
names(G_herb_data) <- c("Site", "ss_data")
G_herb_data$size <- rep(size, length(site))
write.table(G_herb_data, paste("results/GRWTH_herbs.txt"))

G_inv_data <- as.data.frame(cbind(location_info, unlist(Inv_growth_data)))
names(G_inv_data) <- c("Site", "ss_data")
G_inv_data$size <- rep(size, length(site))
write.table(G_inv_data, paste("results/GRWTH_invs.txt"))

Refuge_data <- as.data.frame(cbind(location_info, unlist(Refuge_data)))
names(Refuge_data) <- c("Site", "refs")
Refuge_data$size <- rep(size, length(site))
write.table(Refuge_data, paste("results/refuges.txt"))

#Yield data
#pred_yield <- as.data.frame(cbind(location_info, unlist(Pred_yield)))
#names(pred_yield) <- c("Site", "catch")
#pred_yield$size <- rep(size, length(site))
#write.table(pred_yield, paste("results/July2016/Better_inverts/pred_catch.txt"))

#herb_yield <- as.data.frame(cbind(location_info, unlist(Herb_yield)))
#names(herb_yield) <- c("Site", "catch")
#herb_yield$size <- rep(size, length(site))
#write.table(herb_yield, paste("results/July2016/Better_inverts/herb_catch.txt"))


                                                                                                                                                                                                    