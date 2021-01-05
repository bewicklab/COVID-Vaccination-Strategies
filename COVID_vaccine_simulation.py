#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sunday December 27 2020

@author: sbewick
"""

import numpy as np                      #numerical tools package
import matplotlib.pyplot as plt         #plotting package
import random                           #random number generator
import os                               #output package

#######################################################################################################################
#############################################     DEFINE SYSTEM     ###################################################
#######################################################################################################################

#model parameters
beta_s1 = 0.2                               #per capita per day transmission rate of viral strain 1
beta_s2 = 0.2                               #per capita per day transmission rate of viral strain 2
vaccination_1 = 5000                        #overall number of people per day that receive the vaccination for viral strain 1
vaccination_2 = 0                           #overall number of people per day that receive the vaccination for viral strain 2
efficacy_v1_s1=0.95                         #fraction of viral strain 1 infections prevented by infection by AND/OR the vaccination for viral strain 1 (target protection of the vaccination for viral strain 1)
efficacy_v1_s2 = 0.7                        #fraction of viral strain 2 infections prevented by infection by AND/OR the vaccination for viral strain 1 (cross-protection of the vaccination for viral strain 1)
efficacy_v2_s2 = 0.95                       #fraction of viral strain 2 infections prevented by infection by AND/OR the vaccination for viral strain 2 (target protection of the vaccination for viral strain 2)
efficacy_v2_s1=0.7                          #fraction of viral strain 1 infections prevented by infection by AND/OR the vaccination for viral strain 2 (cross-protection of the vaccination for viral strain 2)
death_rate_v1 = 0.0006                      #per capita per day death rate due to infection by viral strain 1 
death_rate_v2 = 0.0006                      #per capita per day death rate due to infection by viral strain 2
recovery_coefficient = 0.1                  #per capita per day recovery rate from either strain (1/recovery_coefficient = infectious period)
people_0=500000                             #total number of people initially living in the focal city/region
vaccination_acceptance = 1                  #fraction of people willing to get a vaccine

#initial conditions
initial_proportion_s1_immune = 0.1          #initial proportion of people in the focal region that have natural immunity to viral strain 1
initial_proportion_s2_immune = 0            #initial proportion of people in the focal region that have natural immunity to viral strain 2
initial_proportion_s1_infected = 0.001      #initial proportion of people in the focal region that are infected with viral strain 1
initial_proportion_s2_infected = 0.00001    #initial proportion of people in the focal region that are infected with viral strain 2
initial_proportion_v1 = 0                   #initial proportion of people in the focal region that have been vaccinated with the vaccine for viral strain 1
initial_proportion_v2 = 0                   #initial proportion of people in the focal region that have been vaccinated with the vaccine for viral strain 2

#simulation parameters
n_trials = 10                   #number of trials (this is a stochastic simulation... each trial will be different, and average behavior is reported over a number of trials)
max_timeperiod = 2000           #maximum number of days over which to run a single trial (ideally the simulation would finish before reaching this time-point... unless the virus remains uncontrolled)


#vary one of the parameters to make a plot...
for slider in range(0,11):
    
    #to co-vary cross-protection of the vaccines/natural immunity, use both of the following lines
    efficacy_v1_s2 = 0.1*slider
    efficacy_v2_s1 = 0.1*slider
    
    #to vary vaccination rates, use one or the other of the following lines
    #vaccination_1 = (slider+1)*500
    #vaccination_2 = (slider+1)*500
    
    #to vary transmissibility of the second viral variant, use the following line
    #beta_s2 = 0.2 + slider*0.02
    
    #to vary the initial proportion with immunity against strain 1 or the relative proportions with immunity against strain 1 vs strain 2, use the first or both of the following lines 
    #initial_proportion_s1_immune = 0.1 - 0.01*slider
    #initial_proportion_s2_immune = 0.01*slider
    
    #to vary the initial proportion infected with strain 1 versus strain 2, use the following two lines
    #initial_proportion_s1_infected = 0.001+slider*0.000002
    #initial_proportion_s2_infected = 0.00001-slider*0.000002
    
    #calculate 4X4X4 matrices of transmission probabilities, death rates and recovery rates based on infection/immunity/vaccination status
    #in all of the matrices that follow, we use the following coding:
    #the first index describes infection status: 0 - not currently infected, 1 - currently infected with viral strain 1, 2 - currently infected with viral strain 2, 3 - currently infected with both viral strain 1 and viral strain 2
    #the second index describes natural immunity status: 0 - no natural immunity, 1 - natural immunity to viral strain 1, 2 - natural immunity to viral strain 2, 3 - natural immunity to both viral strain 1 and viral strain 2
    #the third index describes vaccination status: 0 - no vaccinations, 1 - received vaccination to viral strain 1, 2 - received vaccination to viral strain 2, 3 - received vaccination to both viral strain 1 and viral strain 2
    
    t_rate1=np.zeros((4,4,4))               #initialize the matrix that describes probability of acquiring viral strain 1 
    t_rate2=np.zeros((4,4,4))               #initialize the matrix that describes probability of acquiring viral strain 2
    death_coeff=np.zeros((4,4,4))           #initialize the matrix that describes probability of death due to infection
    recovery_coeff=np.zeros((4,4,4))        #initialize the matrix that describes probability of recovering from infection
    
    for h in range(0,4):                                                #for each possible infection status...
        for i in range(0,4):                                            #for each possible immunity status...
            for j in range(0,4):                                        #for each possible vaccination stauts...
            
                #per capita probability of acquiring viral strain 1
                if h == 0 or h == 2:                                    #if you currently do NOT have viral strain 1...           
                    if i == 1 or j == 1 or i == 3 or j == 3:            #if you are naturally immune or have been vaccinated against viral strain 1 OR if you are naturally immune or have been vaccinated against both viral strains...
                        t_rate1[h,i,j]=beta_s1*(1-efficacy_v1_s1)       #...reduce per capita infection rate based on the efficacy of strain 1 immunity/vaccination against strain 1
                    elif i == 2 or j == 2:                              #if you are naturally immune or have been vaccinated against viral strain 2...
                        t_rate1[h,i,j]=beta_s1*(1-efficacy_v2_s1)       #...reduce per capita infection rate based on the efficacy of strain 2 immunity/vaccination against strain 1
                    else:                                               #if you have no natural immunity and have not been vaccinated...
                        t_rate1[h,i,j] = beta_s1                        #...assume the full per capita infection rate of viral strain 1
                        
                #per capita probability of acquiring viral strain 2
                if h == 0 or h == 1:                                    #if you currently do NOT have viral strain 2  
                    if i == 2 or j == 2 or i == 3 or j == 3:            #if you are naturally immune or have been vaccinated against viral strain 2 OR if you are naturally immune or have been vaccinated against both viral strains...
                        t_rate2[h,i,j]=beta_s2*(1-efficacy_v2_s2)       #...reduce per capita infection rate based on the efficacy of strain 2 immunity/vaccination against strain 2
                    elif i == 1 or j == 1:                              #if you are naturally immune or have been vaccinated against viral strain 1...
                        t_rate2[h,i,j]=beta_s2*(1-efficacy_v1_s2)       #...reduce per capita infection rate based on the efficacy of strain 1 immunity/vaccination against strain 2
                    else:                                               #if you have no natural immunity and have not been vaccinated...
                        t_rate2[h,i,j] = beta_s2                        #...assume the full per capita infection rate of viral strain 2
                        
                #per capita probability of dying
                if h == 1:                                              #if you are currently infected with viral strain 1...
                    death_coeff[h,i,j]=death_rate_v1                    #...assume the death rate associated with viral strain 1 (regardless of vaccination or immunity status)
                elif h == 2:                                            #if you are currently infected with viral strain 2...
                    death_coeff[h,i,j]=death_rate_v2                    #...assume the death rate associated with viral strain 2 (regardless of vaccination or immunity status)
                elif h == 3:                                            #if you are currently infected with both viral strains...
                    death_coeff[h,i,j]=death_rate_v2+death_rate_v1      #...assume the combined death rate (regardless of infection or immunity status)
                
                #probability of recovering
                if h > 0:                                               #if you are currently infected with at least one viral strain...
                    recovery_coeff[h,i,j]=recovery_coefficient          #...assume the recovery rate (regardless of which strain you are infected with, or vaccination or immunity status)
    
    
    #define vectors that will store all the information that will be collected across trials
    infections_summary=[]                   #total number of new infections over the period of simulation
    infections1_summary=[]                  #total number of new infections with viral strain 1 over the period of simulation
    infections2_summary=[]                  #total number of new infections with viral strain 2 over the period of simulation
    infections_trajectories=[]              #timeseries of new infections over the period of simulation
    infections1_trajectories=[]             #timeseries of new infections with viral strain 1 over the period of simulation
    infections2_trajectories=[]             #timeseries of new infections with viral strain 2 over the period of simulation
    time_trajectories=[]
    deaths_summary=[]                       #total number of deaths over the period of simulation
    max_infections_summary=[]
    time_summary=[]                         #time before herd immunity is acquired (no virus present in the focal region)
    final_state=[]                          #final matrix of the system
    
    
    #######################################################################################################################
    #################################     RUN SIMULATION OVER A NUMBER OF TRIALS     ######################################
    #######################################################################################################################
    
    for trials in range(0,n_trials):                    #cycle through each trial, storing output values at the end
    
        print(trials)
        
    #######################################################################################################################
    ##########################################     INITIALIZE THE TRIAL     ###############################################
    #######################################################################################################################
    
        
        #define initial numbers of susceptible, infected, immune and vaccinated individuals 
        #PLEASE NOTE: we do not consider any initial 'joint classes' - i.e. currently infected with viral strain 1 and naturally immune to viral strain 2, currently immune to viral strain 1 and viral strain 2, etc. 
        #This is a good assumption provided infection, immunity and vaccination rates are still quite low at the beginning of the simulation or, at the very least, are low for at one of the two viral strains
        people = people_0                                               #reset the total number of people
        N=np.zeros( (4, 4,4) )                                          #initialize the matrix that stores the current number of individuals in each infection/immunity/vaccination class
        N[1,0,0] = round(initial_proportion_s1_infected*people)         #initialize the number infected with viral strain 1 at the beginning of the simulation          
        N[2,0,0] = round(initial_proportion_s2_infected*people)         #initialize the number infected with viral strain 2 at the beginning of the simulation
        N[0,0,1] = round(initial_proportion_v1*people)                  #initialize the number vaccinated against strain 1 at the beginning of the simulation
        N[0,0,2] = round(initial_proportion_v2*people)                  #initialize the number vaccinated against strain 2 at the beginning of the simulation
        N[0,1,0] = round(initial_proportion_s1_immune*people)           #initialize the number naturally immune to strain 1 at the beginning of the simulation
        N[0,2,0] = round(initial_proportion_s2_immune*people)           #initialize the number naturally immune to strain 2 at the beginning of the simulation
        N[0,0,0] = people-np.sum(N)                                     #initialize the number fully susceptible at the beginning of the simulation (everyone else)
        t=0                                                             #initialize the time at t = 0
    
        #define variables that store the final output data for this simulation trial
        death_total = 0                                                 #initialize the total number of deaths
        infection_total = 0                                             #initialize the total number of infections
        infection1_total = 0                                            #initialize the total number of infections due to viral strain 1
        infection2_total = 0                                            #initialize the total number of infections due to viral strain 2
        
        #define vectors that store the trajectory data for this simulation trial
        t_values=[]             #vector of times at each timepoint (i.e., the times when 'events' occur)
        infected_values=[]      #vector of the total number of infected individuals at each timepoint
        infected1_values=[]     #vector of the total number of individuals infected with viral strain 1 at each timepoint
        infected2_values=[]     #vector of the total number of individuals infected with viral strain 2 at each timepoint
    
        
        #record the initial values at timepoint t = 0 into the trajectory vectors
        t_values.append(t)                                                              #timepoint t = 0
        infected_values.append(np.sum(N[1,:,:])+np.sum(N[2,:,:])+np.sum(N[3,:,:]))      #total number of infected individuals at t = 0
        infected1_values.append(np.sum(N[1,:,:])+np.sum(N[3,:,:]))                      #total number of individuals infected with viral strain 1 at t = 0
        infected2_values.append(np.sum(N[2,:,:])+np.sum(N[3,:,:]))                      #total number of individuals infected with viral strain 2 at t = 0
        
    
    #######################################################################################################################
    #############################################     RUN THE TRIAL     ###################################################
    #######################################################################################################################
    
        
        #run the simulation until either there are no infected individuals or until the maximum allowable time
        while np.sum(N[1:3,:,:]) > 0  and t < max_timeperiod:
    
    
    #######################################################################################################################
    ###################################     UPDATE PROBABILITIES OF EACH EVENT     ########################################
    #######################################################################################################################
    
            
            #calculate the total number of infected individuals
            infected1=np.sum(N[1,:,:])+np.sum(N[3,:,:])                 #number infected with viral strain 1
            infected2=np.sum(N[2,:,:])+np.sum(N[3,:,:])                 #number infected with viral strain 2
            
            #calculate the overall transmission probabilities for each viral strain to each infection/immunity/vaccination class 
            #(per capita probability of acquiring virus for the class) X (current number of individuals in the class) X (current number of infected individuals) / (total number of people)
            #PLEASE NOTE: this assumes frequency-dependent transmission
            transmission1=infected1*np.multiply(t_rate1,N)/people       #overall transmission probabilities for viral strain 1
            transmission2=infected2*np.multiply(t_rate2,N)/people       #overall transmission probabilities for viral strain 2
    
            #calculate the overall recovery probabilities
            #(per capita probability of recoverying for the class) X (current number of individuals in the class) 
            recovery=np.multiply(recovery_coeff,N)
            
            #calculate the overall probability of dying
            #(per capita probability of dying for the class) X (current number of individuals in the class)
            death=np.multiply(death_coeff,N)
            
            #define a matrix of the number of individuals in each unvaccinated class
            unvaccinated=N[:,:,0]   #all infection classes and immunity classes, but only the unvaccinated class
            
            #define a matrix of the number of individuals in each infected class
            sick=N[1:4,:,:]     #all the immunity classes and vaccination classes, but only the classes that are infected with viral strain 1, viral strain 2 or both viral strains
            
            #find overall probabilities for each 'type' of transition (viral transmissions, recovery, death, vaccination)
            transmission1_rate=np.sum(transmission1)    #transmission of viral strain 1
            transmission2_rate=np.sum(transmission2)    #transmission of viral strain 2
            recovery_rate=np.sum(recovery)              #viral recovery
            death_rate=np.sum(death)                    #death due to virus
            
            #if everyone who accepts the vaccine has been vaccinated (the virus may still be spreading, since vaccination efficacy may be low)...
            if np.sum(unvaccinated)<people*(1-vaccination_acceptance)+1:
                vaccination_1_on = 0                #do not perform any more vaccinations
                vaccination_2_on = 0
            #if there are still some unvaccinated people willing to get the vaccine...
            else:
                vaccination_1_on = vaccination_1    #continue with vaccinations
                vaccination_2_on = vaccination_2
                
            
            #define a vector of each of the different process-level rates (transmission rates for viral strains 1 and 2, death rate, recovery rate, vaccination rates for viral strain 1 and 2)
            rates=[transmission1_rate,vaccination_1_on,recovery_rate,transmission2_rate,death_rate,vaccination_2_on]
            #sum ALL of the rates for all of the different processes
            srates=np.sum(rates)
            #define a normalized, cumulative sum of all the process-level rates... the first entry tells you the probability of the first event, the second tells you the probability of the first OR second event, the third tells you the probability of the first OR second OR third... etc.
            rate_thresholds=np.cumsum(rates)/srates
    
    #######################################################################################################################
    ##############################     DETERMINE THE NEXT EVENT AND UPDATE CLASSES     ####################################
    #######################################################################################################################
            
    
            #use an exponentially-distributed random variable to pick the 'time' at which the next event occurs 
            timestep = (np.log(1/np.random.uniform()))/srates;       
            #update the current time by adding the time to the next event
            t=t+timestep
            #add this timepoint to the vector of timepoints for the current simulation
            t_values.append(t)
        
            #decide which event occurs at this time by picking a uniform random number between 0 and 1 and then...
            event_dice=np.random.uniform()
            #...figuring out where the random number falls in the rate_thresholds vector that defines the probabilities of each process happening
            which_event=np.argmax(rate_thresholds > event_dice)
            
            #if the event is a transmission of viral strain 1...
            if which_event == 0:                                            
                nzte=np.nonzero(transmission1)                              #find the elements in the tranmission matrix for viral strain 1 where there are non-zero transmission probabilities
                                                                            #this step gives three vectors of equal length.  The first, second and third vectors specify the infection status, immunity status and vaccination status respectively for each non-zero viral strain 1 transmission probability
                nzt=transmission1[nzte]                                     #define a vector of the non-zero tranmission probabilities for viral strain 1... these will be in the same order as the element vector from the line above
                rate_thresholds_t1=np.cumsum(nzt)/transmission1_rate        #define a normalized, cumulative sum of all the non-zero probabilities associated with transmission of viral strain 1... again, events will be in the same order as the two lines above
                event_t_dice=np.random.uniform()                            #decide which transmission event occurs by picking a uniform random number between 0 and 1 and then...
                which_t_event=np.argmax(rate_thresholds_t1>event_t_dice)    #...figuring out where the random number falls in the rate_thresholds_t1 vector that defines the probabilities of each event associated with transmission of viral strain 1
                #use the element vector from line 255 to decide which class of individuals just got infected
                infection_class=nzte[0][which_t_event]                      #determine the (prior) infection status of the newly infected individual (from the first vector generated on line 255)          
                immunity_class=nzte[1][which_t_event]                       #determine the (prior) immunity status of the newly infected individual (from the second vector generated on line 255)
                vaccination_class=nzte[2][which_t_event]                    #determine the (prior) vaccination status of the newly infected individual (from the third vector generated on line 255)
                N[infection_class,immunity_class,vaccination_class]=N[infection_class,immunity_class,vaccination_class]-1   #remove the newly infected individual from their previous class
                if infection_class == 0 or infection_class == 1:                                    #if they were already infected with viral strain 1, or they were not infected with any virus...
                    N[1,immunity_class,vaccination_class]=N[1,immunity_class,vaccination_class]+1   #...they preserve their prior infection and vaccination statuses, but now become infected with virus 1 (or remain infected with virus 1, if they were already infected)
                else:                                                                               #if they were already infected with viral strain 2, or they were already infected with both viral strains...
                    N[3,immunity_class,vaccination_class]=N[3,immunity_class,vaccination_class]+1   #...they preserve their prior infection and vaccination statuses, but now become infected with both viral strains (or remain infected with both viral strains, if they already had both)
                if infection_class == 0:                                        #if they were not already infected with a strain of the virus...
                    infection_total=infection_total+1                           #...add a new infection to the total infection count
                if infection_class == 2 or infection_class == 0:                #if they were not already infected with strain 1 of the virus...
                    infection1_total=infection1_total+1                         #...add a new infection to the total count of infections with viral strain 1
                    
            #if the event is vaccination with the vaccine for viral strain 1...
            elif which_event == 1:
                nzve=np.nonzero(unvaccinated)                                               #find the elements of the unvacinnated matrix that contain non-zero numbers of unvaccinated individuals
                                                                                            #this step gives two vectors of equal length. The first and second vectors specify the infection status and immunity status respectively for each unvaccinated class with non-zero individuals
                nzv=unvaccinated[nzve]                                                      #find the number of individuals in each unvaccinated class with a non-zero number of individuals... the probability of vaccination going to each unvaccinated class will be proportional to the number in each unvaccinated class
                nzvs=np.sum(nzv)                                                            #find the total number of unvaccinated individuals across all classes
                rate_thresholds_v1=np.cumsum(nzv)/nzvs                                      #define a normalized, cumulative sum of all the non-zero probabilities of each unvaccinated class receiving a vaccination
                event_v_dice=np.random.uniform()                                            #decide which vaccination event occurs by picking a uniform random number between 0 and 1 and then...
                which_v_event=np.argmax(rate_thresholds_v1>event_v_dice)                    #...figuring out where the random number falls in the rate_thresholds_v1 vector that defines the probabilities of each event associated with vaccination using the vaccine for viral strain 1
                #use the element vector from line 277 to decide which class of individuals just got vaccinated
                infection_class=nzve[0][which_v_event]                                      #determine the (prior) infection status of the newly vaccinated individual (from the first vector generated on line 277)
                immunity_class=nzve[1][which_v_event]                                       #determine the (prior) immunity status of the newly vaccinated individual (from the second vector generated on line 277)
                N[infection_class,immunity_class,0]=N[infection_class,immunity_class,0]-1   #remove the newly vaccinated individual from their previous class
                N[infection_class,immunity_class,1]=N[infection_class,immunity_class,1]+1   #preserve their infection status and their natural immunity status, but move them to the class vaccinated with the vaccine for viral strain 1
                
            #if the event is viral recovery...
            elif which_event == 2:
                nzre=np.nonzero(recovery)                                       #find the elements of the recovery matrix that are non-zero
                                                                                #this step gives three vectors of equal length.  The first, second and third vectors specify the infection status, immunity status and vaccination status respectively for each non-zero recovery probability
                nzr=recovery[nzre]                                              #define a vector of the non-zero recovery probabilities
                rate_thresholds_r=np.cumsum(nzr)/recovery_rate                  #define a normalized, cumulative sum of all the non-zero recovery probabilities     
                event_r_dice=np.random.uniform()                                #decide which recovery event occurs by picking a uniform random number between 0 and 1 and then...
                which_r_event=np.argmax(rate_thresholds_r>event_r_dice)         #...figuring out where the random number falls in the rate_threshold_r vector that defines the probabilities of each recovery event
                #use the element vector from line 292 to decide which class of individuals just recovered
                infection_class=nzre[0][which_r_event]                          #determine the (prior) infection status of the newly recovered individual (from the first vector generated on line 292)
                immunity_class=nzre[1][which_r_event]                           #determine the (prior) immunity status of the newly recovered individual (from the second vector generated on line 292)
                vaccination_class=nzre[2][which_r_event]                        #determine the (prior) vaccination status of the newly recovered individual (from the third vector generated on line 292)
                N[infection_class,immunity_class,vaccination_class]=N[infection_class,immunity_class,vaccination_class]-1   #remove the newly recovered individual from their previous class
                if immunity_class == 0 or infection_class == immunity_class:                            #if they were previously susceptible, or else already naturally immune to the viral strain they just recovered from...
                    N[0,infection_class,vaccination_class]=N[0,infection_class,vaccination_class]+1     #... move them to the class with no current infection, but with natural immunity to the viral strain they just recovered from; preserve their vaccination status
                elif immunity_class == 1 and infection_class == 2:                                      #if they were previously naturally immune to viral strain 1, but just recovered from viral strain 2...
                    N[0,3,vaccination_class]=N[0,3,vaccination_class]+1                                 #...move them to the class with no current infection and natural immunity to both viral strains; preserve their vaccination status
                elif immunity_class == 2 and infection_class == 1:                                      #if they were previously naturally immune to viral strain 2, but just recovered from viral strain 1...
                    N[0,3,vaccination_class]=N[0,3,vaccination_class]+1                                 #...move them to the class with no current infection and natural immunity to both viral strains; preserve their vaccination status
                elif immunity_class == 3 or infection_class == 3:                                       #if they were previously naturally immune to both viral strains or just recovered from both viral strains...
                    N[0,3,vaccination_class]=N[0,3,vaccination_class]+1                                 #...move them to the class with no current infection and natural immunity to both viral strains; preserve their vaccination status
                else:                                                                                   #if they don't fall into any category, print a warning.
                    print('some class is misisng')
                    print(infection_class)
                    print(immunity_class)
                
            #if the event is a transmission of viral strain 2...
            elif which_event == 3:
                nzte=np.nonzero(transmission2)                              #find the elements in the tranmission matrix for viral strain 2 where there are non-zero transmission probabilities
                                                                            #this step gives three vectors of equal length.  The first, second and third vectors specify the infection status, immunity status and vaccination status respectively for each non-zero viral strain 2 transmission probability
                nzt=transmission2[nzte]                                     #define a vector of the non-zero tranmission probabilities for viral strain 2... these will be in the same order as the element vector from the line above
                rate_thresholds_t2=np.cumsum(nzt)/transmission2_rate        #define a normalized, cumulative sum of all the non-zero probabilities associated with transmission of viral strain 2... again, events will be in the same order as the two lines above
                event_t_dice=np.random.uniform()                            #decide which transmission event occurs by picking a uniform random number between 0 and 1 and then...
                which_t_event=np.argmax(rate_thresholds_t2>event_t_dice)    #...figuring out where the random number falls in the rate_thresholds_t2 vector that defines the probabilities of each event associated with transmission of viral strain 2
                #use the element vector from line 318 to decide which class of individuals just got infected
                infection_class=nzte[0][which_t_event]                      #determine the (prior) infection status of the newly infected individual (from the first vector generated on line 318)          
                immunity_class=nzte[1][which_t_event]                       #determine the (prior) immunity status of the newly infected individual (from the first vector generated on line 318)          
                vaccination_class=nzte[2][which_t_event]                    #determine the (prior) vaccination status of the newly infected individual (from the first vector generated on line 318)          
                N[infection_class,immunity_class,vaccination_class]=N[infection_class,immunity_class,vaccination_class]-1   #remove the newly infected individual from their previous class
                if infection_class == 0 or infection_class == 2:                                        #if they were already infected with viral strain 2, or they were not infected with any virus...
                    N[2,immunity_class,vaccination_class]=N[2,immunity_class,vaccination_class]+1       #...they preserve their prior infection and vaccination statuses, but now become infected with virus 2 (or remain infected with virus 2, if they were already infected)
                else:                                                                                   #if they were already infected with viral strain 1, or they were already infected with both viral strains...
                    N[3,immunity_class,vaccination_class]=N[3,immunity_class,vaccination_class]+1       #...they preserve their prior infection and vaccination statuses, but now become infected with both viral strains (or remain infected with both viral strains, if they already had both)
                if infection_class == 0:                                        #if they were not already infected with a strain of the virus...
                    infection_total=infection_total+1                           #...add a new infection to the total infection count
                if infection_class == 1 or infection_class == 0:                #if they were not already infected with strain 2 of the virus...
                    infection2_total=infection2_total+1                         #...add a new infection to the total count of infections with viral strain 2
                  
            #if the event is a death...
            elif which_event == 4:
                nzse=np.nonzero(death)                                          #find the elements of the death rate matrix that contain non-zero probabilities
                                                                                #this step gives three vectors of equal length. The first, second vectors and third specify the infection status, immunity status and vaccination status respectively for each infected class with non-zero probabilities of dying
                nzs=death[nzse]                                                 #define a vector of the non-zero probabilities of dying
                rate_thresholds_s=np.cumsum(nzs)/death_rate                     #define a normalized, cumulative sum of all the non-zero probabilities of dying     
                event_s_dice=np.random.uniform()                                #decide which death event occurs by picking a uniform random number between 0 and 1 and then...
                which_s_event=np.argmax(rate_thresholds_s>event_s_dice)         #...figuring out where the random number falls in the rate_threshold_s vector that defines the probabilities of each death event
                #use the element vector from line 340 to decide which class of individuals just died
                infection_class=nzse[0][which_s_event]                          #determine the infection status of the dead individual (from the first vector generated on line 340)          
                immunity_class=nzse[1][which_s_event]                           #determine the immunity status of the dead individual (from the second vector generated on line 340)          
                vaccination_class=nzse[2][which_s_event]                        #determine the vaccination status of the dead individual (from the third vector generated on line 340)          
                N[infection_class,immunity_class,vaccination_class]=N[infection_class,immunity_class,vaccination_class]-1   #remove the dead individual from their previous class
                people=people-1                                                 #decrease the number of people in the region by 1
                death_total=death_total+1                                       #increase the death count by one
                
            #if the event is vaccination with the vaccine for viral strain 2...
            elif which_event == 5:
                nzve=np.nonzero(unvaccinated)                                   #find the elements of the unvacinnated matrix that contain non-zero numbers of unvaccinated individuals
                                                                                #this step gives two vectors of equal length. The first and second vectors specify the infection status and immunity status respectively for each unvaccinated class with non-zero individuals
                nzv=unvaccinated[nzve]                                          #find the number of individuals in each unvaccinated class with a non-zero number of individuals... the probability of vaccination going to each unvaccinated class will be proportional to the number in each unvaccinated class
                nzvs=np.sum(nzv)                                                #find the total number of unvaccinated individuals across all classes
                rate_thresholds_v2=np.cumsum(nzv)/nzvs                          #define a normalized, cumulative sum of all the non-zero probabilities of each unvaccinated class receiving a vaccination
                event_v_dice=np.random.uniform()                                #decide which vaccination event occurs by picking a uniform random number between 0 and 1 and then...
                which_v_event=np.argmax(rate_thresholds_v2>event_v_dice)        #...figuring out where the random number falls in the rate_thresholds_v2 vector that defines the probabilities of each event associated with vaccination using the vaccine for viral strain 2
                #use the element vector from line 356 to decide which class of individuals just got vaccinated
                infection_class=nzve[0][which_v_event]                          #determine the (prior) infection status of the newly vaccinated individual (from the first vector generated on line 356)
                immunity_class=nzve[1][which_v_event]                           #determine the (prior) immunity status of the newly vaccinated individual (from the second vector generated on line 356)
                N[infection_class,immunity_class,0]=N[infection_class,immunity_class,0]-1   #remove the newly vaccinated individual from their previous class
                N[infection_class,immunity_class,2]=N[infection_class,immunity_class,2]+1   #preserve their infection status and their natural immunity status, but move them to the class vaccinated with the vaccine for viral strain 2
            
    #######################################################################################################################
    #################################     RECORD TRAJECTORIES AND FINAL OUTPUTS    ########################################
    #######################################################################################################################
            
            
            infected_values.append(np.sum(N[1,:,:])+np.sum(N[2,:,:])+np.sum(N[3,:,:]))  #sum of total number of individuals infected with viral variant 1, viral variant 2 or both viral variants at the current timepoint
            infected1_values.append(np.sum(N[1,:,:])+np.sum(N[3,:,:]))                  #sum of the total number of individuals infected with viral variant 1 at the current timepoint
            infected2_values.append(np.sum(N[2,:,:])+np.sum(N[3,:,:]))                  #sum of the total number of individuals infected with viral variant 2 at the current timepoint
    
        #record summaries of the timeseries
        deaths_summary.append(death_total)                          #total deaths
        max_infections_summary.append(np.max(infected_values))      #peak infections over the whole trajectory
        infections_summary.append(infection_total)                  #total infections
        infections1_summary.append(infection1_total)                #total infections with viral variant 1
        infections2_summary.append(infection2_total)                #total infections with viral variant 2
        time_summary.append(t)                                      #time until the virus went extinct
        
        #subsample the trajectory
#        vvs=np.where(np.mod(t_values,1)<0.0005)                    #do not record every timepoint (to save space, selectively sample timepoints)
#        t_trunc=[]                                                 #sampled timepoints
#        in_trunc=[]                                                #sampled total infections
#        in1_trunc=[]                                               #sampled infections with variant 1
#        in2_trunc=[]                                               #sampled infections with variant 2
#        for k in range(0,len(vvs[0])):
#            t_trunc.append(t_values[vvs[0][k]])
#            in_trunc.append(infected_values[vvs[0][k]])
#            in1_trunc.append(infected1_values[vvs[0][k]])
#            in2_trunc.append(infected2_values[vvs[0][k]])
         
        #store the subsampled trajectories for the current simulation trial   
 #       time_trajectories.append(t_trunc)
 #       infections_trajectories.append(in_trunc)
 #       infections1_trajectories.append(in1_trunc)
 #       infections2_trajectories.append(in2_trunc)
 #       final_state.append(N)
        
    

    #print output to a series of csv files, one for each parameter value scanned    
    parameter_summary=str(beta_s1)+','+str(beta_s2)+','+str(vaccination_1)+','+str(vaccination_2)+','+str(efficacy_v1_s1)+','+str(efficacy_v1_s2)+','+str(efficacy_v2_s2)+','+str(efficacy_v2_s1)+','+str(death_rate_v1)+','+str(death_rate_v2)+','+str(recovery_coefficient)+','+str(people)+','+str(vaccination_acceptance)+','+str(initial_proportion_s1_immune)+','+str(initial_proportion_s2_immune)+','+str(initial_proportion_s1_infected)+','+str(initial_proportion_s2_infected)+','+str(initial_proportion_v1)+','+str(initial_proportion_v2)
    summary_printer=[]
    summary_printer.append(parameter_summary+'\n')
    summary_printer.append('deaths,max_infections,infections,infections1,infections2,time\n')
    
    for u in range(0,len(deaths_summary)):
        summary_printer.append(str(deaths_summary[u])+','+str(max_infections_summary[u])+','+str(infections_summary[u])+','+str(infections1_summary[u])+','+str(infections2_summary[u])+','+str(time_summary[u])+'\n')
    
    if os.path.exists('COVID_vaccine_data'+str(slider)+'.csv'):
        os.remove('COVID_vaccine_data'+str(slider)+'.csv')
    
    j=open('COVID_vaccine_data'+str(slider)+'.csv','w')
    for u in range(0,len(summary_printer)):
        j.write(summary_printer[u])
    j.close()
    
