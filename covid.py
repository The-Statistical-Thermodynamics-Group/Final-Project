import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import random
from space_ops import *

class Person(object):
    PERSON_CONFIG = {
    "infection_radius": 0.5,
    "probability_of_infection": 0.5,
    "counted_as_infector": False,
    "incubation_period": 1.04,
    "infection_duration": 4.0,
    "dl_bound": [0, 0],
    "ur_bound": [10, 10],
    "wander_step_size": 1,
    "wander_step_duration": 1,
    "gravity_well": None,
    "gravity_strength": 1,
    "social_distance_factor": 1,
    "probability_of_social_distancing": 0.7,
    "n_repulsion_points":10,
    "wall_buffer": 1,
    "max_speed": 1,
    "dt": 0.1
    }
    

    def __init__(self):
        self.time = 0
        self.last_step_change = -1
        self.velocity = np.zeros(2)
        self.status = "S"
        self.repulsion_points = []
        self.num_infected = 0
        self.getting_infected_time = 0
        
        self.infection_radius = self.PERSON_CONFIG["infection_radius"]
        self.probability_of_infection = self.PERSON_CONFIG["probability_of_infection"]
        self.counted_as_infector = self.PERSON_CONFIG["counted_as_infector"]
        self.incubation_period = self.PERSON_CONFIG["incubation_period"]
        self.infection_duration = self.PERSON_CONFIG["infection_duration"]
        self.dl_bound = self.PERSON_CONFIG["dl_bound"]
        self.ur_bound = self.PERSON_CONFIG["ur_bound"]
        self.wander_step_size = self.PERSON_CONFIG["wander_step_size"]
        self.wander_step_duration = self.PERSON_CONFIG["wander_step_duration"]
        self.gravity_well = self.PERSON_CONFIG["gravity_well"]
        self.gravity_strength = self.PERSON_CONFIG["gravity_strength"]
        self.social_distance_factor = self.PERSON_CONFIG["social_distance_factor"]
        self.probability_of_social_distancing = self.PERSON_CONFIG["probability_of_social_distancing"]
        self.n_repulsion_points = self.PERSON_CONFIG["n_repulsion_points"]
        self.wall_buffer = self.PERSON_CONFIG["wall_buffer"]
        self.max_speed = self.PERSON_CONFIG["max_speed"]
        self.dt = self.PERSON_CONFIG["dt"]
        
        self.set_point()
    
    
    def set_point(self):
        right_bound = float(self.ur_bound[0])
        left_bound = float(self.dl_bound[0])
        upper_bound = float(self.ur_bound[1])
        lower_bound = float(self.dl_bound[1])
        
        x = random.uniform(left_bound, right_bound)
        y = random.uniform(lower_bound, upper_bound)
        self.point = np.array([x, y])
        return self.point
        
    

    def update_status(self, status):
        self.status = status
    
    
    def update_position(self, dt=None):
        if dt is None:
            dt = self.dt
        point = self.point
        total_force = np.zeros(2)

        # Gravity
        if self.wander_step_size != 0:
            if (self.time - self.last_step_change) > self.wander_step_duration:
                vect = rotate_vector(np.array([1, 0]), angle = 2*np.pi * random.random())
                self.gravity_well = point + self.wander_step_size * vect
                self.last_step_change = self.time

        if self.gravity_well is not None:
            to_well = self.gravity_well - point
            dist = get_norm(to_well)
            if dist != 0:
                total_force += self.gravity_strength * to_well / dist**3
        
        
        # Avoid walls
        wall_force = np.zeros(2)
        for i in range(2):
            # When i=0, to_lower and to_upper mean the horizontal distances
            # When i=1, to_lower and to_upper mean the vertical distances
            to_lower = point[i] - self.dl_bound[i]
            to_upper = self.ur_bound[i] - point[i]
            
            # Bounce
            if to_lower < 0:
                self.point[i] = self.dl_bound[i]
                self.velocity[i] = abs(self.velocity[i])
            if to_upper < 0:
                self.point[i] = self.ur_bound[i]
                self.velocity[i] = -abs(self.velocity[i])
            
            # Repelling force
            wall_force += max((-1 / self.wall_buffer + 1 / to_lower), 0)
            wall_force -= max((-1 / self.wall_buffer + 1 / to_upper), 0)
        total_force += wall_force
        
        
        # Potentially avoid neighbors (Optional)
        if self.social_distance_factor > 0 and random.random() < self.probability_of_social_distancing:
            repulsion_force = np.zeros(2)
            min_dist = np.inf
            for repulsion_point in self.repulsion_points:
                to_repulsion_point = repulsion_point - point
                dist = get_norm(to_repulsion_point)
                if 0 < dist < min_dist:
                    min_dist = dist
                if dist > 0:
                    repulsion_force -= self.social_distance_factor * to_repulsion_point / dist**3

            total_force += repulsion_force
        
        # Apply force
        self.velocity += total_force * dt

        # Limit speed
        speed = get_norm(self.velocity)
        if speed > self.max_speed:
            self.velocity *= self.max_speed / speed
        
        # Update position
        self.point += self.velocity * dt
    
    
    def update_infection_ring(self, infection_radius):
        self.infection_radius = infection_radius
    

    def update_time(self, dt=None):
        if dt is None:
            dt = self.dt
        self.time += dt


class SIRSimulation(Person):
	CONFIG = { #not finished yet
	"n_cities":6, #the number of the cities in this sim
	"row_cities":2,
	"col_cities":3, #define the row and column of the cities
	"city_pop":100, #the population in each city
	"box_size":10, #the size of each city(box)
	"travel_rate":0.3, #decide how much people in the city will travel
    "trigger_case":30,
    "num_of_total_infected":0,
    "percentage_of_no_symptom": 0.2,
    "percentage_of_quarantine": 0.7,
    "quarantine_mode": False
	}
	def __init__(self): #need to figure out how to update the status
            super().__init__()
            self.n_cities=self.CONFIG["n_cities"]
            self.row_cities=self.CONFIG["row_cities"]
            self.col_cities=self.CONFIG["col_cities"]
            self.city_pop=self.CONFIG["city_pop"]
            self.box_size=self.CONFIG["box_size"]
            self.travel_rate=self.CONFIG["travel_rate"]
            self.trigger_case = self.CONFIG["trigger_case"]
            self.num_of_total_infected = self.CONFIG["num_of_total_infected"]
            self.percentage_of_no_symptom = self.CONFIG["percentage_of_no_symptom"]
            self.percentage_of_quarantine = self.CONFIG["percentage_of_quarantine"]
            self.quarantine_mode = self.CONFIG["quarantine_mode"]
            
            self.add_box()
            self.add_people()
    
    
	def add_box(self):
            boxes = []
            for n in range(self.n_cities):
                box = [] #a single city
                boxes.append(box)
            self.quarantine_zone = []
            self.boxes = boxes
            return self.boxes, self.quarantine_zone
    
    
	def add_people(self):
            for city in self.boxes:
                for i in range(self.city_pop):
                    person = Person()
                    city.append(person)
            for city in self.boxes:
                #choose the patient zero
                random.choice(city).update_status("I")
                self.CONFIG["num_of_total_infected"] += 1
    
    
	def update_statuses(self):
            cities = self.boxes
            for city in cities:
                s_group=[]
                i_group=[]
                for citizen in city:
                    if citizen.status == "S":
                        s_group.append(citizen)
                    elif citizen.status == "I":
                        i_group.append(citizen)
                    # H person 
                    elif citizen.status == "A":
                        i_group.append(citizen)
                
                for s_person in s_group:
                    for i_person in i_group:
                        #the distance between s_person and i_person
                        dist = get_norm(s_person.point - i_person.point)
                        if dist < s_person.infection_radius and random.random() < s_person.probability_of_infection:
                            #s_person have a chance to become an A
                            if random.random() < self.percentage_of_no_symptom:
                                s_person.update_status("A")
                                s_person.getting_infected_time = s_person.time
                                self.CONFIG["num_of_total_infected"] += 1
                                i_person.num_infected += 1
                            #if s_person got an infection than update its status to I
                            else:
                                s_person.update_status("I")
                                s_person.getting_infected_time = s_person.time
                                self.CONFIG["num_of_total_infected"] += 1
                                i_person.num_infected += 1
                
                for i_person in i_group:
                    #need to add a timing when the infection start
                    if (i_person.time - i_person.getting_infected_time) > i_person.infection_duration:
                        #if the infection last longer than infection_duration, update its status to R
                        i_person.update_status("R")
            
            #quarantine
            for city in self.boxes:
                for person in city:
                    if person.time - person.getting_infected_time > self.incubation_period:
                        if person.status is "I" and self.quarantine_mode is True:
                            if random.random() < self.percentage_of_quarantine:
                                self.quarantine_zone.append(person)
                                del city[city.index(person)]
            
            for i_person in self.quarantine_zone:
                    #need to add a timing when the infection start
                    if (i_person.time - i_person.getting_infected_time) > i_person.infection_duration:
                        #if the infection last longer than infection_duration, update its status to R
                        i_person.update_status("R")
            
            #travel
            if self.travel_rate > 0:
                for city in self.boxes:
                    if random.random() < self.travel_rate:
                        # dest for destination
                        dest = round(self.n_cities*random.random()) - 1
                        traveler = round(len(city)*random.random())
                        self.boxes[dest].append(city[traveler])
                        del city[traveler] # need to make sure what del_count exactly do then add it
            
            #social distancing
            if self.CONFIG["num_of_total_infected"] > self.trigger_case:
                for city in self.boxes:
                    points = np.array([person.point for person in city])
                    repelled_points = points
                    for point, person in zip(points, city):
                        if person.social_distance_factor > 0:
                            diffs = np.linalg.norm(repelled_points - point, axis=1)
                            person.repulsion_points = repelled_points[np.argsort(diffs)[1:person.n_repulsion_points + 1]]
            
	def get_status_count(self):
            status_counter = np.zeros(4)
            for city in self.boxes:
                for person in city:
                    if person.status is "S":
                        status_counter[0] += 1
                    elif person.status is "I":
                        status_counter[1] += 1
                    elif person.status is "A":
                        status_counter[1] += 1
                        status_counter[2] += 1
                    elif person.status is "R":
                        status_counter[3] += 1
            for person in self.quarantine_zone:
                if person.status is "I":
                        status_counter[1] += 1
                elif person.status is "R":
                        status_counter[3] += 1
            return status_counter
        


class RunSimpleSimulation(SIRSimulation):
    def __init__(self):
        super().__init__()
        self.last_update_time = 0
        self.effect_reproduction_num = 0
        self.particles = dict()
        self.real_world_time_list = []
        self.percentage_S_list = []
        self.percentage_I_list = []
        self.percentage_R_list = []
        self.percentage_accum_cases_list = []

        self.setup()
        self.anim = animation.FuncAnimation(self.fig, func= self.run_until_zero_infection, interval = 100)
        self.plot_anim = animation.FuncAnimation(self.plot_fig, func= self.run_until_zero_infection, interval = 100)
        
        
    def setup(self):
        for city in range(self.n_cities):
            self.particles[city+1] = dict()
            self.particles[city+1]["S"] = np.zeros([0,2])
            self.particles[city+1]["I"] = np.zeros([0,2])
            self.particles[city+1]["R"] = np.zeros([0,2])
            self.particles[city+1]["A"] = np.zeros([0,2])
            for n in range(len(self.boxes[city])):
                person = self.boxes[city][n]
                if person.status is "S":
                    self.particles[city+1]["S"] = np.append(self.particles[city+1]["S"],[person.point],0)
                if person.status is "I":
                    self.particles[city+1]["I"] = np.append(self.particles[city+1]["I"],[person.point],0)
                if person.status is "R":
                    self.particles[city+1]["R"] = np.append(self.particles[city+1]["R"],[person.point],0)
                if person.status is "A":
                    self.particles[city+1]["A"] = np.append(self.particles[city+1]["A"],[person.point],0)
        
        # Scatter plot
        city_plot = dict()
        group_plot = dict()
        fig = plt.figure(figsize=(12,8))
        for city in range(self.n_cities):
            city_plot[city+1] = plt.subplot(5,3,city+1)
            city_plot[city+1].set_xlim(left=0,right=self.box_size)
            city_plot[city+1].set_ylim(bottom=0,top=self.box_size)
            group_plot[city+1] = dict()
            group_plot[city+1]["S"] = city_plot[city+1].scatter(self.particles[city+1]["S"][:,0],self.particles[city+1]["S"][:,1], color = "blue")
            group_plot[city+1]["I"] = city_plot[city+1].scatter(self.particles[city+1]["I"][:,0],self.particles[city+1]["I"][:,1], color = "red")
            group_plot[city+1]["A"] = city_plot[city+1].scatter(self.particles[city+1]["A"][:,0],self.particles[city+1]["A"][:,1], color = "yellow")
            group_plot[city+1]["R"] = city_plot[city+1].scatter(self.particles[city+1]["R"][:,0],self.particles[city+1]["R"][:,1], color = "grey")
            
        #quarantine
        if self.quarantine_mode is True:
            self.Q_particles = dict()
            self.Q_particles["I"] = np.zeros([0,2])
            self.Q_particles["R"] = np.zeros([0,2])
            for person in self.quarantine_zone:
                if person.status is "I":
                    self.Q_particles["I"] = np.append(self.Q_particles["I"],[person.point],0)
                elif person.status is "R":
                    self.Q_particles["R"] = np.append(self.Q_particles["R"],[person.point],0)
                
            quarantine_plot = plt.subplot(5,3,7)
            plt.tick_params(length=0,labelsize=0)
            quarantine_plot.annotate("Quarantine zone",xy = (0.5,8.6))
            quarantine_plot.set_xlim(left=0,right=self.box_size)
            quarantine_plot.set_ylim(bottom=0,top=self.box_size)
            self.quarantine_scatter = dict()
            self.quarantine_scatter["I"] = quarantine_plot.scatter(self.Q_particles["I"][:,0], self.Q_particles["I"][:,1], color = "red")
            self.quarantine_scatter["R"] = quarantine_plot.scatter(self.Q_particles["R"][:,0], self.Q_particles["R"][:,1], color = "grey")
        
       # Show the reproduction number
        tex = plt.subplot(5,3,10)
        plt.tick_params(length=0,labelsize=0)
        
        tex.annotate("Time in code:",xy=(0.1,0.76),fontsize=17)
        time_label = tex.annotate(self.boxes[1][1].time,xy=(0.8,0.76),fontsize=17)
        
        tex.annotate("Days:",xy=(0.1,0.5),fontsize=17)
        Days = tex.annotate(int(self.boxes[1][1].time/0.2),xy=(0.8,0.5),fontsize=17)
        
        tex.annotate("R:",xy=(0.1,0.16),fontsize=17)
        R_label = tex.annotate(self.effect_reproduction_num,xy=(0.8,0.16),fontsize=17)
        
        self.totol_people = self.get_status_count()[0]+self.get_status_count()[1]+self.get_status_count()[3]
        group_area = plt.subplot(5,3,11)
        plt.tick_params(length=0,labelsize=0)
        group_area.annotate("Susceptible:",xy=(0.1,0.85),fontsize=10)
        S_displayer = group_area.annotate(self.get_status_count()[0],xy=(0.45,0.85),fontsize=10)
        S_percent = group_area.annotate('{:.0%}'.format(self.get_status_count()[0]/(self.totol_people)),xy=(0.7,0.85),fontsize=10)
        
        group_area.annotate("Infectious:",xy=(0.1,0.67),fontsize=10)
        I_displayer = group_area.annotate(self.get_status_count()[1],xy=(0.45,0.67),fontsize=10)
        I_percent = group_area.annotate('{:.0%}'.format(self.get_status_count()[1]/(self.totol_people)),xy=(0.7,0.67),fontsize=10)
         
        group_area.annotate("Asymptotic:",xy=(0.1,0.45),fontsize=10)
        A_displayer = group_area.annotate(self.get_status_count()[2],xy=(0.45,0.45),fontsize=10)
        A_percent = group_area.annotate('{:.0%}'.format(self.get_status_count()[2]/(self.totol_people)),xy=(0.7,0.45),fontsize=10)
        
        group_area.annotate("Removed:",xy=(0.1,0.25),fontsize=10)
        R_displayer = group_area.annotate(self.get_status_count()[3],xy=(0.45,0.25),fontsize=10)
        R_percent = group_area.annotate('{:.0%}'.format(self.get_status_count()[3]/(self.totol_people)),xy=(0.7,0.25),fontsize=10)
        
        self.S_displayer = S_displayer
        self.I_displayer = I_displayer
        self.A_displayer = A_displayer
        self.R_displayer = R_displayer
        
        self.S_percent = S_percent
        self.I_percent = I_percent
        self.A_percent = A_percent
        self.R_percent = R_percent
        
        # Show the CONFIG info
        info = plt.subplot(5,3,13)
        plt.tick_params(length=0,labelsize=0)
        #info of infection radius
        info.annotate("#Infection radius:",xy=(0.03,0.85),fontsize=10)
        info.annotate(self.infection_radius,xy=(0.8,0.85),fontsize=10)
        #info of probability of infection
        info.annotate("#Probability of infection:",xy=(0.03,0.6),fontsize=10)
        info.annotate(self.probability_of_infection,xy=(0.8,0.6),fontsize=10)
        #info of incubation period
        info.annotate("#Incubation period:",xy=(0.03,0.36),fontsize=10)
        info.annotate(self.incubation_period,xy=(0.8,0.36),fontsize=10)
        # info of the infection duration
        info.annotate("#Infection duration:",xy=(0.03,0.1),fontsize=10)
        info.annotate(self.infection_duration,xy=(0.8,0.1),fontsize=10)
        
        
        info_2 = plt.subplot(5,3,14)
        plt.tick_params(length=0,labelsize=0)
        # info of the travel rate
        info_2.annotate("#Travel rate:",xy=(0.03,0.85),fontsize=10)
        info_2.annotate(self.travel_rate,xy=(0.8,0.85),fontsize=10)
        # info of the social distance factor
        info_2.annotate("#Social distance factor:",xy=(0.03,0.67),fontsize=10)
        info_2.annotate(self.social_distance_factor,xy=(0.8,0.67),fontsize=10)
        #info of the percentage of social distancing
        info_2.annotate("#Percentage of \nsocial distancing:",xy=(0.03,0.3),fontsize=10)
        info_2.annotate(self.probability_of_social_distancing,xy=(0.8,0.3),fontsize=10)
        #info of trigger case
        info_2.annotate("#Trigger case:",xy=(0.03,0.1),fontsize=10)
        info_2.annotate(self.trigger_case,xy=(0.8,0.1),fontsize=10)


        info_3 = plt.subplot(5,3,15)
        plt.tick_params(length=0,labelsize=0)
        #info of trigger case
        info_3.annotate("#Percentage of no symptom:",xy=(0.03,0.8),fontsize=10)
        info_3.annotate(self.percentage_of_no_symptom,xy=(0.8,0.8),fontsize=10)
        #info of Quarantine mode
        info_3.annotate("#Quarantine mode:",xy=(0.03,0.5),fontsize=10)
        info_3.annotate(self.quarantine_mode,xy=(0.8,0.5),fontsize=10)
        #info of Quarantine mode
        info_3.annotate("#Percentage of quarantine:",xy=(0.03,0.2),fontsize=10)
        info_3.annotate(self.percentage_of_quarantine,xy=(0.8,0.2),fontsize=10)
        
        
        # Initiate the plot
        self.plot_fig, self.plot_ax = plt.subplots(figsize=(8,5))
        self.plot_ax.set_title("Real-world Time - Percentage of Population", fontweight= "bold", fontsize= 16)
        self.plot_ax.set_xlabel("Real-world Time (days)", fontsize= 12)
        self.plot_ax.set_ylabel("Percentage of Population", fontsize= 12)
        self.plot_ax.set_prop_cycle(color= ['blue', 'red', 'gray', 'orange'])
        self.plot_ax.set_xlim(0, auto= True)
        self.plot_ax.set_ylim(0, 1)
        # Parameter
        real_world_time = self.time / 0.2
        percentage_of_S = self.get_status_count()[0] / self.totol_people
        percentage_of_I = self.get_status_count()[1] / self.totol_people
        percentage_of_R = self.get_status_count()[3] / self.totol_people
        percentage_of_accumulated_cases = percentage_of_I + percentage_of_R
        # time axis (x axis)
        self.real_world_time_list.append(real_world_time)
        # y axis (number of people)
        self.percentage_S_list.append(percentage_of_S)
        self.percentage_I_list.append(percentage_of_I)
        self.percentage_R_list.append(percentage_of_R)
        self.percentage_accum_cases_list.append(percentage_of_accumulated_cases)
        self.plot_ax.plot(self.real_world_time_list, self.percentage_S_list, label= "Susceptible")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_I_list, label= "Infectious")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_R_list, label= "Removed")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_accum_cases_list, label= "Accumulated cases")
        self.plot_ax.legend()
        
        
        '''
        stack_plot = plt.plot(x= self.time_list, y= self.percentage_of_individuals,
        labels= ["Susceptible", "Infectious", "Removed"], colors= ["blue", "red", "grey"], baseline='zero')
        '''

        self.time_label = time_label
        self.R_label = R_label
        self.Days = Days
        self.fig = fig
        self.city_plot = city_plot
        self.group_plot = group_plot


        '''
        self.stack_fig = stack_fig
        self.stack_plot = stack_plot
        '''
        return self.fig, self.city_plot, self.group_plot, self.R_label, self.time_label
        
    
    def run_until_zero_infection(self, frame):
        for city in range(self.n_cities):
            self.particles[city+1] = dict()
            self.particles[city+1]["S"] = np.zeros([0,2])
            self.particles[city+1]["I"] = np.zeros([0,2])
            self.particles[city+1]["R"] = np.zeros([0,2])
            self.particles[city+1]["A"] = np.zeros([0,2])
            for n in range(len(self.boxes[city])):
                person = self.boxes[city][n]
                person.update_position()
                person.update_time()
                if person.status is "S":
                    self.particles[city+1]["S"] = np.append(self.particles[city+1]["S"],[person.point],0)
                if person.status is "I":
                    self.particles[city+1]["I"] = np.append(self.particles[city+1]["I"],[person.point],0)
                if person.status is "A":
                    self.particles[city+1]["A"] = np.append(self.particles[city+1]["A"],[person.point],0)
                if person.status is "R":
                    self.particles[city+1]["R"] = np.append(self.particles[city+1]["R"],[person.point],0)
            
            print("status[S, I, A, R]", self.get_status_count())
        
        #update quarantine zone
        if self.quarantine_mode is True:
            self.Q_particles = dict()
            self.Q_particles["I"] = np.zeros([0,2])
            self.Q_particles["R"] = np.zeros([0,2])
            for person in self.quarantine_zone:
                person.update_position()
                person.update_time()
                if person.status is "I":
                    self.Q_particles["I"] = np.append(self.Q_particles["I"],[person.point],0)
                elif person.status is "R":
                    self.Q_particles["R"] = np.append(self.Q_particles["R"],[person.point],0)
            self.quarantine_scatter["I"].set_offsets(self.Q_particles["I"])
            self.quarantine_scatter["R"].set_offsets(self.Q_particles["R"])
        
        for city in range(self.n_cities):
            self.group_plot[city+1]["S"].set_offsets(self.particles[city+1]["S"])
            self.group_plot[city+1]["I"].set_offsets(self.particles[city+1]["I"])
            self.group_plot[city+1]["A"].set_offsets(self.particles[city+1]["A"])
            self.group_plot[city+1]["R"].set_offsets(self.particles[city+1]["R"])
        
        
        # Collect the population list.
        real_world_time = self.time / 0.2
        percentage_of_S = self.get_status_count()[0] / self.totol_people
        percentage_of_I = self.get_status_count()[1] / self.totol_people
        percentage_of_R = self.get_status_count()[3] / self.totol_people
        percentage_of_accumulated_cases = percentage_of_I + percentage_of_R
        # time axis (x axis)
        self.real_world_time_list.append(real_world_time)
        # y axis (number of people)
        self.percentage_S_list.append(percentage_of_S)
        self.percentage_I_list.append(percentage_of_I)
        self.percentage_R_list.append(percentage_of_R)
        self.percentage_accum_cases_list.append(percentage_of_accumulated_cases)
        self.plot_ax.plot(self.real_world_time_list, self.percentage_S_list, label= "Susceptible")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_I_list, label= "Infectious")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_R_list, label= "Removed")
        self.plot_ax.plot(self.real_world_time_list, self.percentage_accum_cases_list, label= "Accumulated cases")


        '''
        # Update stacked area plot
        # time axis (x axis)
        self.time_list.append(self.time)
        # y axis (number of people)
        self.num_S_list = np.append(self.num_S_list, self.get_status_count()[0])
        self.num_I_list = np.append(self.num_I_list, self.get_status_count()[1])
        self.num_R_list = np.append(self.num_R_list, self.get_status_count()[2])
        self.num_of_individuals = np.vstack([self.num_S_list, self.num_I_list, self.num_R_list])
        
        self.stack_plot.set_data(x= self.time_list, y1= self.num_S_list, y2= self.num_I_list, y3= self.num_R_list)
        '''

        # If number of infected people is zero, stop the simulation.
        if (self.get_status_count()[1] + self.get_status_count()[2]) == 0:
            self.anim.event_source.stop()
            self.plot_anim.event_source.stop()


        self.update_time()
        self.update_statuses()
        self.add_R_label()
        self.R_label.set_text(round(self.effect_reproduction_num,1))
        self.Days.set_text(int(self.boxes[1][1].time/0.2))
        self.time_label.set_text(round(self.boxes[1][1].time,1))
        self.S_displayer.set_text(int(self.get_status_count()[0]))
        self.I_displayer.set_text(int(self.get_status_count()[1]))
        self.A_displayer.set_text(int(self.get_status_count()[2]))
        self.R_displayer.set_text(int(self.get_status_count()[3]))
        
        self.S_percent.set_text('{:.0%}'.format(self.get_status_count()[0]/(self.totol_people)))
        self.I_percent.set_text('{:.0%}'.format(self.get_status_count()[1]/(self.totol_people)))
        self.A_percent.set_text('{:.0%}'.format(self.get_status_count()[2]/(self.totol_people)))
        self.R_percent.set_text('{:.0%}'.format(self.get_status_count()[3]/(self.totol_people)))
        
        return self.group_plot,
    
    
    def add_R_label(self):
        update_period = 1
        all_R0_values = []
        
        if (self.time - self.last_update_time) > update_period:
            return
        self.last_update_time = self.time
        # cities
        for city in self.boxes:
            for person in city:
                if person.status == "I" or person.status == "A":
                    prop = (person.time - person.getting_infected_time) / person.infection_duration
                    # When time passes beyond the interval of 0.1, update the "reproduction number."
                    if prop > 0.1:
                        all_R0_values.append(person.num_infected / prop)
        # quarantine zone
        for person in self.quarantine_zone:
            if person.status == "I":
                prop = (person.time - person.getting_infected_time) / person.infection_duration
                # When time passes beyond the interval of 0.1, update the "reproduction number."
                if prop > 0.1:
                    all_R0_values.append(person.num_infected / prop)
                    
        if len(all_R0_values) > 0:
            # The "effect_reproduction_num"  is for R_label.
            self.effect_reproduction_num = np.mean(all_R0_values)
    

    def show(self):
        self.fig.savefig("Simulation Animation")
        plt.show()
            

sim = RunSimpleSimulation()
sim.show()