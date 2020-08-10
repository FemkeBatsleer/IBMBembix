'''
Created on 6 okt. 2017

@author: Femke Batsleer
'''
import matplotlib.image as mpimg
import matplotlib
matplotlib.use('Agg')
import tkinter as tk
import numpy as np
import math
import colorsys
from copy import deepcopy
import os
import datetime
import time

save_path = 'data'

'''
******************************************************************************
Parameterization of the model
First, needed field data are loaded to add in the model
(like distance histogram, individual distances and periods between nests etc)
or functions defined (density response function)
******************************************************************************
'''
#distance kernel data
load_dist = np.loadtxt('Hist_distances.txt', delimiter='\t', usecols=(1,), skiprows=1)
#histogram is per 2 m
freq_dist = [i/np.sum(load_dist) for i in load_dist] #frequency distribution of distances
freq_dist_norm = [i/freq_dist[0] for i in freq_dist] #normalised frequency distribution: chances will add up to one

#starting day frequencies
p_startingday = [0.0306,0.0036,0.018,0.0503,0.0306,0.0198,0.0288,0.0288,0.0521,0.0234,0.0396,0.0521,0.0521,0.0414,0.0449,0.0144,0.0378,0.0413,0.0432,0.0629,0.036,0.0324,0.0629,0.0198,0.018,0.0306,0.0072,0.0252,0.0288,0.0234]
p_period = [0.0704,0.0916,0.1056,0.1127,0.0916,0.0845,0.0775,0.0493,0.0775,0.0423,0.0211,0.0352,0.007,0.0423,0.007,0.0211,0.0141,0.0211,0,0,0.0141,0,0,0.007,0,0,0.007,0,0,0]

#starting amount of nests chances, already adjusted for the fact that during the modeling, not all number of nests will be finished
#thus, several wasps starting with 2 nests will and up with only 1 nest during running period
p_number_nests = [0.691, 0.216, 0.079, 0.014]

#response function for density: type 2
def response_type2(dens, a=0.2, h=1): #the higher a, the higher the chance gets with lower densities
    resp = a*dens/(1+a*h*dens) +0.05 #type 2 (holling's disc equation)
    return resp if resp <= 1 else 1

#response function for density: parabolic
def response_par(dens, a=math.pi*4*4, h=1):
    parabole = (-4/(a**2))*dens**2 + (4/a)*dens #parabolic function with max y=1 and NP(0,0) and (0,a=max_dens)
    return parabole
    

'''
******************************************************************************
Definition of classes
Second, classes are defined to make Wasp, Nest and Population objects
class Visual is to make a visual output
******************************************************************************
'''
 
# class Visual:
#     '''This class arranges the visual output.'''
#     def __init__(self, max_x, max_y):
#         '''Initialize the visual class'''
#         self.zoom = 2
#         self.max_x = max_x
#         self.max_y = max_y
#         self.root = tk.Tk()
#         self.canvas = tk.Canvas(self.root, 
#                                 width =  self.max_x * self.zoom, 
#                                 height = self.max_y * self.zoom) #create window
#         self.canvas.pack()
#         self.canvas.config(background = 'white')
#         self.squares = np.empty((self.max_x, self.max_y),dtype=object)
#         self.initialize_squares()
#           
#     def create_nest(self,nest, day, amount_days):       
#         '''Create circle for nest'''
#         x, y = nest.x, nest.y
#         radius = 1
#         #generating colour according to day
#         color = float(day/float(amount_days))
#         if color < 0:
#             color = 0
#         elif color > 1:
#             color = 1
#         #color range including all colours is best in hls system    
#         red, green, blue = colorsys.hls_to_rgb(color, 0.5, 1)
#         rgb = int(red*255), int(green*255), int(blue*255)
#         hex_code = '#%02x%02x%02x' % rgb
#           
#         return self.canvas.create_oval((x - radius) * self.zoom,
#                                        (y - radius) * self.zoom,
#                                        (x + radius) * self.zoom,
#                                        (y + radius) * self.zoom,
#                                        outline=str(hex_code), 
#                                        fill=str(hex_code))                                      
#           
#     def color_square(self, resources, x, y):
#         '''Colors the square relative to the suitability of that point'''
#         '''white is 1, black is 0'''
#         #trying gray values: rgb are equal to each other then        
#         color = float(resources)
#         if color < 0:
#             color = 0
#         elif color > 1:
#             color = 1
#   
#         blue = int(255*color)    
#         green = int(255*color)
#         red = int(255*color)     
#           
#         rgb = red, green, blue     
#         hex_code = '#%02x%02x%02x' % rgb
#         #edge errors in the image... value >=204
#         if color >= 0.8:
#             hex_code = '#%02x%02x%02x' % (255,250,205)
#         self.canvas.itemconfigure(self.squares[x, y],fill = str(hex_code))
#           
#     def initialize_squares(self):
#         '''returns a square (drawing object)'''
#         for x in range(self.max_x):
#             for y in range(self.max_y):
#                 self.squares[x, y] = self.canvas.create_rectangle(self.zoom * x,
#                                                      self.zoom * y, 
#                                                      self.zoom * x + self.zoom,
#                                                      self.zoom * y + self.zoom,
#                                                      outline = '', 
#                                                      fill = 'black')

class Nest:
    """
    Creates nest objects with a position and day of initialisation
    """
    def __init__(self, x, y, day):
        self.x = x
        self.y = y
        self.day=day
        
    def __repr__(self):
        return 'Nest{}'.format(id(self))
        
    def __str__(self):
        return '({},{})'.format(self.x, self.y)
        

class Wasp:
    """
    Class wasp, which creates a Bembix rostrata individual
    has function search, which lets the wasp search for a nest according to giving scenario
    this search function uses other help-functions to do the search according to the scenario
    """
    def __init__(self, starting_day, amount_nests, period, scenario, pixel_dist, id_real_wasp):
        #variables important for initialisation
        self.starting_day = starting_day #starting day of the wasp
        self.day_next_nest = starting_day #day of next nest building
        self.amount_nests = amount_nests #amount of nests the wasp will make in total
        #self.distances=distances #list with distances between subsequent nests
        self.periods = period #list with periods between the subsequent nests
        
        #variable important for complete population
        self.pixel_dist = pixel_dist #conversion factor to change pixels into metres; [amount of metres/pixels]
        self.scenario=scenario #scenario of the population: random, bottom-up, prev-exp...
#         self.dist_hist=dist_hist #boolean, True if distances (and periods) between subsequent nests are extracted from the distance histogram
# #                                False if distances (and periods) of one individual are drawn from field data
        self.id_real_wasp=id_real_wasp #saves the key to enter self.distances and self.periods, drawn from field data
        
        #variables that are properties of the wasp
        self.nests = [] #list with objects Nests for this individual
        
    def search(self, environment, population, amount_days, active, CA_w_env, CA_day, response_func, bi_CA, weight_day):
        '''
        Function that selects correct function that searches by a specific method, according to the scenario
        '''
        if self.scenario == 'random':
            self.search_random(environment)
        elif self.scenario == 'bottom-up':
            self.search_BU(environment)
        elif self.scenario == 'prev-exp':
            self.search_PE(environment)
        elif self.scenario == 'prev-exp-BU':
            self.search_PE_BU(environment)
        elif self.scenario == 'conspec-attraction':
            self.search_CA(environment, population,  amount_days=amount_days,active=active, CA_w_env=CA_w_env, CA_day=CA_day,
                           response_func=response_func, bi_CA=bi_CA, weight_day=weight_day)
        elif self.scenario == 'bimodal':
            self.search_bimodal(environment, population,  amount_days=amount_days,
                                active=active, CA_w_env=CA_w_env, CA_day=CA_day,
                                response_func=response_func, bi_CA=bi_CA, weight_day=weight_day)
            
    def calc_kernel_PE(self, environment):
        '''
        Helper function for PE BU scenario
        Here the distance kernel is combined with the environment
        New environment is created
        '''
        x_w, y_w = self.nests[-1].x, self.nests[-1].y
        new_environment = deepcopy(environment) 

        for i in range(len(environment[0])): #run over x-values
            for j in range(len(environment)): #run over y-values
                dist = math.sqrt((i-x_w)**2 + (j-y_w)**2) #calculate distance between last nest and looped position/pixel: pythagoras
                dist_m = dist*self.pixel_dist #metres
                if new_environment[j,i] >=0.8: #if larger than 0.8, outside environment
                    new_environment[j,i] = 1
                elif dist_m > 82:#distance kernel is only up to 82 metres
                    new_environment[j,i] = 1
                elif self.scenario=='bimodal' and dist_m > 5:#if it is larger than bimodial distance, not taken into account
                        new_environment[j,i] = 1
                else:
                    #find back frequency in freq_dist_norm for current dist; normalised so chances add up to 1
                    freq_here = freq_dist_norm[int(dist_m/2)] #divided by two because histogram is per 2 m
                    if environment[j,i] < 3/255:
                        environment[j,i] = 3/255
                    new_environment[j,i] = environment[j,i]*freq_here #new environment suitability is kernel_value*(actual suitability)
        return new_environment
    
    def calc_kernel_CA(self, x, y, population, active, amount_days, CA_day, response_func, weight_day):
        '''
        helper function for conspecific attraction
        This calculates the amount of individuals present in radius 2m
        Returns a chance, depending on amount of individuals present
        '''
        #calculating distance to every nest (active nests only, or not active, all nests present
        #calculate distance to each wasp, is smaller than 2m, add to density
        density = 0 #init density: sum of amount of nests within 2m
        for wasp in population: #loop over all wasps
            if wasp.nests:#if the wasps in the population already have nests
                if active: #if active is true, only active nests are taken into account
                    nests = [wasp.nests[-1]]
                else:
                    nests = wasp.nests
                #loop over every nest    
                for nest in nests:  
                    dist = math.sqrt((nest.x - x)**2 + (nest.y-y)**2) #calculate distance between given position and looped nest
                    dist_m = dist*self.pixel_dist #in metres
                    if dist_m < 5:
                        if CA_day:
                            diff_day = abs(self.day_next_nest - wasp.nests[-1].day)/(amount_days) #calculate difference in days, and normalised
                            diff_day_weighted = weight_day*math.exp(-diff_day) #inverse diff_day, lowers expontentially
                            density += diff_day_weighted
                        else:
                            density += 1

        chance = response_func(density, a=0.8 if response_func==response_type2 else math.pi*(5**2)*5) #chance is calculated with density response function
        
        if density==0:
            chance=0.05 #if then bottom-up is taken into account: only that will be taken into account; otherwise it's random
        
        #print(chance)
        return chance if chance > 0 else 0

    def search_random(self, environment):
        '''
        Here the wasp search for a good nesting place
        Here it is random in the available space.
        '''
        y_rn, x_rn = np.random.uniform(0,len(environment)), np.random.uniform(0,len(environment[0])) #select random y and x values
        while environment[int(y_rn), int(x_rn)] >= 0.8: #check if it is inside the focal landscape, while not, choose another random position
            y_rn, x_rn = np.random.uniform(0,len(environment)), np.random.uniform(0,len(environment[0]))
        
        self.dig(x_rn, y_rn)
        
    def search_BU(self, environment):
        '''
        Here the wasp search for a good nesting place
        It searches by picking a rondam pixel, and by chance relative to the suitability,
        chooses to stay or not
        The BU_force indicates the strength of the bottom-up choosing: the higher this factor, the more choosy the wasp is
        --> first runs showed that this parameter did not really matter (prior was same as posterior)
        '''
        chosen = False #init chosen, if chosen is True, while loop is exited
        while not chosen:
            y_rn, x_rn = np.random.uniform(0,len(environment)), np.random.uniform(0,len(environment[0])) #select random position
            if environment[int(y_rn), int(x_rn)] < 0.8: #has to be inside focal field
                chance = environment[int(y_rn), int(x_rn)]
                chosen = bool(np.random.choice((1,0), p=[chance**2, 1-chance**2])) #BU_force = 3
        self.dig(x_rn, y_rn)
    
    def search_PE(self, environment):
        '''
        Here the wasp search for a good nesting place
        Here it is random in the available space. But taking into account prev experience.
        Environment suitability is not taken into account
        '''
        #if previous nest was build, distance to next nest by kernel from field study
        if self.nests:
            #choose direction
            dir_rn = np.random.uniform(0, 2*math.pi)
            #choose distance, second term in numerator is to make dist_ld continuous
            dist_ld = (np.random.choice(np.arange(1,42), p=freq_dist) - #term is subtracted which will make distance continuous
                        np.random.uniform(0,1))*2/self.pixel_dist #in pixels

            #new pixel coordinates: old pixel coordinates + sin or cos fo direction * distance (elementary goniometry)
            i, j = self.nests[-1].x + math.sin(dir_rn)*dist_ld, self.nests[-1].y + math.cos(dir_rn)*dist_ld
            
            #next positions should be in the range of the image & inside field itselfn otherwise, new direction is generated+ positions calculated
            while (
                (0 > i or i >= len(environment[0])) or #i and j must be inside environment
                (0 > j or j >= len(environment)) or
                environment[int(j)][int(i)] > 0.8 #must be in focal field
                ):
                
                dir_rn = np.random.uniform(0, 2*math.pi)
                
                dist_ld = (np.random.choice(np.arange(1,42), p=freq_dist) - 
                            np.random.uniform(0,1))*2/self.pixel_dist
                
                               
                i, j = self.nests[-1].x + math.sin(dir_rn)*dist_ld, self.nests[-1].y + math.cos(dir_rn)*dist_ld
                
            self.dig(i, j)
        
        #if first nest, random search
        else:   
            #select random y value and x value, check if it is really in the landscape
            self.search_random(environment)
         
    def search_PE_BU(self, environment):
        '''
        Good nesting place: taking into account previous experience + bottom up
        '''
        if self.nests:
            new_environment = self.calc_kernel_PE(environment)
            self.search_BU(new_environment)
        else:
            self.search_BU(environment)
            
    def search_CA(self, environment, population, amount_days, active, CA_w_env, CA_day,
                  response_func, bi_CA, weight_day, chance_inside=None):
        '''
        Conspecific attraction search
        Taking into account density of individuals present + BU
        The higher the CA_force, the higher the attraction at low densities is
        '''
        chosen=False
        while not chosen:
            y_rn, x_rn = np.random.uniform(0,len(environment)), np.random.uniform(0,len(environment[0]))
            if chance_inside and bi_CA and self.nests: #if bimodal choice inside cluster, y_rn and x_rn are searched within circle around nest
                while 5 > math.sqrt((self.nests[-1].x-x_rn)**2 + (self.nests[-1].y-y_rn)**2):
                    x_rn = np.random.uniform(self.nests[-1].x - 5, self.nests[-1].x+5)
                    y_rn = np.random.uniform(self.nests[-1].y - 5, self.nests[-1].y+5)
                    
            if environment[int(y_rn)][int(x_rn)] < 0.8:
                chance = self.calc_kernel_CA(x_rn, y_rn, population, active, amount_days=amount_days,
                                            CA_day=CA_day, response_func=response_func, weight_day=weight_day) #calculate chance depending on other individuals present
                if CA_w_env:
                    chance = chance*(environment[int(y_rn)][int(x_rn)]**2) #combine chance with environmental suitability
                chosen = bool(np.random.choice((1,0), p=[chance, 1-chance]))
            
        self.dig(x_rn, y_rn)

    def search_bimodal(self, environment,population, amount_days, active, CA_w_env, CA_day,
                        response_func, bi_CA, weight_day):
        '''
        Bimodal choice
        Takes into account the chance a wasp builds its nest close to its previous one (in the same cluster)
        otherwise, it will select with conspecific attraction, randomly in the study area
        '''
        chance_inside = bool(np.random.choice((1,0), p=[0.559, 1-0.559])) #amount of internal loops in network from data
        #so chance_inside is the chance to make a next nest in the same cluster
        if chance_inside and self.nests and not bi_CA:#if searching inside is not according to CA
            self.search_PE_BU(environment)
        elif chance_inside and self.nests and bi_CA:#if searching is according to CA
            self.search_CA(environment, population, amount_days, active, CA_w_env, CA_day, response_func, bi_CA, weight_day, chance_inside)
            #the chance_inside is given to the search_CA(), which uses this to choose a nest in a certain radius
        else:
            self.search_CA(environment, population, amount_days, active, CA_w_env, CA_day, response_func, bi_CA, weight_day, chance_inside)
        
    def dig(self, x, y):
        """
        A wasp digs a nest
        amount of nests it will make is minus one
        the day of the next nest is added with period between nests
        """
        self.amount_nests -= 1
        self.nests.append(Nest(x, y, self.day_next_nest))
        if self.amount_nests !=0:
            self.day_next_nest += self.periods[-self.amount_nests]
        
        #self.periods = np.random.choice(np.arange(1,28), p=[0.3545,0.1485,0.1060,0.103,0.0576,0.0333,0.0364,0.0212,0.0364,0.0212,0.0091,0.0152,0.003,0.0152,0.003,0.0091,0.0061,0.0091,0,0,0.0061,0,0,0.003,0,0,0.003])

                
    def __str__(self):
        return "{} with starting day={}, next nest day={}, nests to make ={}".format(
            self.__class__.__name__, self.starting_day, self.day_next_nest, self.amount_nests)
    
    def __repr__(self):
        return "Wasp{}".format(id(self))
                    
class Population:
    """
    Class population, which will contain all the wasps and change in time
    """
    
    def __init__(self, amount_ind=20, amount_days=30, scenario='CF', pixel_dist = 0.21,
                active=True, CA_w_env=True, CA_day=False, choosers_on_followers=None,
                response_func=response_type2, bi_CA=False,weight_day=3, environment=mpimg.imread('suitabilityINLAv2.jpg')):
        #variables important for cover all population/scenario
        self.amount_ind = amount_ind
        self.amount_days = amount_days
        self.pixel_dist = pixel_dist
        self.choosers_on_followers = choosers_on_followers #the amount (0-1) of choosers/followers
        self.response_func=response_func
        
        #parameters, important for search functions
        self.scenario= scenario
        self.active=active #if only active nests are taking into account in CA
        self.CA_w_env = CA_w_env #if conspecific attraction, the BU is also taken into account
        self.CA_day = CA_day #defines if day is taking into account as a weight during density-calculations
        self.bi_CA=bi_CA #defines, when a wasp makes its next nest close by, if it also takes into account the surrounding nests there
        self.weight_day = weight_day #defines weight for CA_day
        
        #properties of the population that will change
        self.population = [] #init the list with wasps in the population
        self.amount_nests = 0 #init total amount of nests of the population
        self.day_count = 1

        self.initialise(self.amount_ind) #wasps created and added to the population

        #reading the image of suitability, turning it into an np.array
        #getting first layer of suitability map, was greyscale but R saved it as RGB with 3 times same value
        self.environment = environment[:,:,0]/255

        #visualisation of the landscape
#         self.x_max, self.y_max = self.environment.shape[1], self.environment.shape[0]
#         self.visual = Visual(self.x_max, self.y_max)
#         self.visual.initialize_squares()
#         for y in range(self.y_max):
#             for x in range(self.x_max):
#                 self.visual.color_square(self.environment[y,x], x, y)
       
    def initialise(self, amount_ind):
        '''Wasps created and added to the population'''
        #if histogram of distances and periods are included in the model
        if not self.choosers_on_followers:
            for _ in range(amount_ind):
                amount_nests_chosen = np.random.choice(np.arange(1,5), p=p_number_nests)
                self.population.append(Wasp(starting_day = np.random.choice(np.arange(1,31), p=p_startingday),
                                            amount_nests = amount_nests_chosen,
                                            period = [np.random.choice(np.arange(1,31), p=p_period) for _ in range(1,amount_nests_chosen)],
                                            scenario=self.scenario, pixel_dist=self.pixel_dist,
                                            id_real_wasp=None))
                        
        else:
            amount_choosers = int(self.choosers_on_followers*amount_ind)
            amount_followers = amount_ind - amount_choosers
            
            for _ in range(amount_choosers):
                amount_nests_chosen = np.random.choice(np.arange(1,5), p=p_number_nests)
                self.population.append(Wasp(starting_day = np.random.choice(np.arange(1,31), p=p_startingday),
                                            amount_nests = amount_nests_chosen,
                                            period = [np.random.choice(np.arange(1,31), p=p_period) for _ in range(1,amount_nests_chosen)],
                                            scenario='prev-exp-BU', pixel_dist=self.pixel_dist,
                                            id_real_wasp=None))
                
            for _ in range(amount_followers):
                amount_nests_chosen = np.random.choice(np.arange(1,5), p=p_number_nests)
                self.population.append(Wasp(starting_day = np.random.choice(np.arange(1,31), p=p_startingday),
                                            amount_nests = amount_nests_chosen,
                                            period = [np.random.choice(np.arange(1,31), p=p_period) for _ in range(1,amount_nests_chosen)],
                                            scenario='conspec-attraction', pixel_dist=self.pixel_dist,
                                            id_real_wasp=None))
            
            
            
    def __str__(self):
        return '{} with {} wasps'.format(self.__class__.__name__, self.amount_ind) #+ '\n' + '\n'.join(
            #str(repr(ind)) for ind in self.population)
        
    def __repr__(self):
        return "Population object"
    
    def day(self):
        for ind in self.population:
            #print(ind.scenario) 
            if ind.day_next_nest == self.day_count and ind.amount_nests > 0:#if the day of building its next nest is equal to current day and it has still remaining nests to build
                ind.search(self.environment, self.population, amount_days=self.amount_days,
                           active=self.active, CA_w_env=self.CA_w_env, CA_day=self.CA_day,
                           response_func=self.response_func, bi_CA=self.bi_CA, weight_day=self.weight_day)
#                 self.visual.create_nest(ind.nests[-1], self.day_count, self.amount_days)#nest is visualised
                self.amount_nests += 1 #amount of nests of the population +1
        self.day_count +=1
        
    def let_it_run(self):
        for i in range(1, self.amount_days+1):
            #♣print('DAY {}'.format(i))
            self.day()
            
           
    def create_output(self, number, savepath = save_path):
        '''
        **********************************************************************
        Function to create output
        3 files will be created: parameters, output and distances
        ************************************************************************
        '''
        
        #output file with parameters
        file_name_par = 'Parameters {} {} {}.txt'.format(self.scenario, number, datetime.datetime.now().strftime("%d_%m_%Y %H_%M"))
        writer = open(os.path.join(savepath, file_name_par), 'w')
        writer.write('scenario\t active\tCA_w_env\tCA_day\tresponse_func\tchoosers_on_followers\tbi_CA\tweight_day\n')
        writer.write(self.scenario + '\t'+ str(self.active) +  '\t'+ str(self.CA_w_env) + '\t'+ str(self.CA_day) + '\t'+ self.response_func.__name__ +
                     '\t'+ str(self.choosers_on_followers) + '\t' + str(self.bi_CA)+ '\t' + str(self.weight_day) + '\n')
        writer.close()
        #output file with waspID, nestID, x, y, day
        file_name = 'Output {} {} {}.txt'.format(self.scenario, number, datetime.datetime.now().strftime("%d_%m_%Y %H_%M"))
        writer = open(os.path.join(savepath, file_name), 'w')
        writer.write('waspID\tnestID\tx\ty\tday\n')
        for wasp in self.population:
            for nest in wasp.nests:
                writer.write(repr(wasp) + '\t' + repr(nest) + '\t' + str(nest.x*self.pixel_dist) + '\t' + str(nest.y*self.pixel_dist) + '\t' + str(nest.day) + '\n')
        writer.close()
        
        #output file with waspID, nest1, nest2, distance
        file_name_w = 'Distances {} {} {}.txt'.format(self.scenario, number, datetime.datetime.now().strftime("%d_%m_%Y %H_%M"))
        writer_dist = open(os.path.join(savepath, file_name_w), 'w')
        writer_dist.write('waspID\tnest1\tnest2\tdistance\n')
        for wasp in self.population:
            if len(wasp.nests) > 1:
                for i in range(len(wasp.nests)-1):
                    dist = math.sqrt((wasp.nests[i].x - wasp.nests[i+1].x)**2 + (wasp.nests[i].y - wasp.nests[i+1].y)**2)*self.pixel_dist
                    writer_dist.write(repr(wasp) + '\t' + repr(wasp.nests[i]) + '\t' + repr(wasp.nests[i+1]) + '\t' + str(dist) + '\n')
        writer_dist.close()
        


'''
# ******************************************************************************
# Models are ran and output is created
# ******************************************************************************
# '''

#     #sample params
# active = False#np.random.choice([True, False])
# CA_w_env = True#np.random.choice([True, False])
# CA_day = True#np.random.choice([True, False])
# response_func = response_par#np.random.choice([response_type2, response_par])
# bi_CA = False#np.random.choice([True, False])#defines, when a wasp makes its next nest closeby, if it also takes into account the surrounding nests there
# weight_day = 12#np.random.normal(10,3)
# while weight_day <= 0:
#     weight_day = np.random.normal(10,3)
#     
# pop6 = Population(amount_ind=455, amount_days=30, scenario='bimodal',
#                                active=active, CA_w_env=CA_w_env, CA_day=CA_day, response_func=response_func,
#                                bi_CA=bi_CA, weight_day=weight_day)
# pop6.let_it_run()
# # pop5.create_output(i, savepath='data/Outputs/bimodal/20-11')
#      
# print('Bimodal population {}'.format(1))
# pop6.create_output(1)
# tk.mainloop()
