import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import random
from space_ops import *

class Person(object):
    PERSON_CONFIG = {
    "infection_radius": 0.5,
    "probability_of_infection": 1,
    "counted_as_infector": False,
    "incubation_period": 1.0,
    "infection_duration": 3.0,
    "dl_bound": [0, 0],
    "ur_bound": [10, 10],
    "wander_step_size": 1,
    "wander_step_duration": 1,
    "gravity_well": None,
    "gravity_strength": 1,
    "social_distance_factor": 0.6,
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
        if self.social_distance_factor > 0:
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
	"n_cities":3, #the number of the cities in this sim
	"row_cities":2,
	"col_cities":3, #define the row and column of the cities
	"city_pop":100, #the population in each city
	"box_size":10, #the size of each city(box)
	"travel_rate":0.1, #decide how much people in the city will travel
    "trigger_case":15,
    "num_of_total_infected":0,
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
            
            self.add_box()
            self.add_people()
    
    
	def add_box(self):
            boxes = []
            for n in range(self.n_cities):
                box = [] #a single city
                boxes.append(box)
            self.boxes = boxes
            return self.boxes
    
    
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
                
                for s_person in s_group:
                    for i_person in i_group:
                        #the distance between s_person and i_person
                        dist = get_norm(s_person.point - i_person.point)
                        if dist < s_person.infection_radius and random.random() < s_person.probability_of_infection:
                            #if s_person got an infection than update its status to I
                            s_person.update_status("I")
                            s_person.getting_infected_time = s_person.time
                            self.CONFIG["num_of_total_infected"] += 1
                            i_person.num_infected += 1
                
                for i_person in i_group:
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
            status_counter = np.zeros(3)
            for city in self.boxes:
                for person in city:
                    if person.status is "S":
                        status_counter[0] += 1
                    elif person.status is "I":
                        status_counter[1] += 1
                    elif person.status is "R":
                        status_counter[2] += 1
            return status_counter
        


class RunSimpleSimulation(SIRSimulation):
    def __init__(self):
        super().__init__()
        self.last_update_time = 0
        self.effect_reproduction_num = 0
        self.particles = dict()
        self.setup()
        self.ani = animation.FuncAnimation(self.fig, func= self.run_until_zero_infection, interval = 100)
        
    def setup(self):
        for city in range(self.n_cities):
            self.particles[city+1] = dict()
            self.particles[city+1]["S"] = np.zeros([0,2])
            self.particles[city+1]["I"] = np.zeros([0,2])
            self.particles[city+1]["R"] = np.zeros([0,2])
            for n in range(len(self.boxes[city])):
                person = self.boxes[city][n]
                if person.status is "S":
                    self.particles[city+1]["S"] = np.append(self.particles[city+1]["S"],[person.point],0)
                if person.status is "I":
                    self.particles[city+1]["I"] = np.append(self.particles[city+1]["I"],[person.point],0)
                if person.status is "R":
                    self.particles[city+1]["R"] = np.append(self.particles[city+1]["R"],[person.point],0)
        
        city_plot = dict()
        group_plot = dict()
        fig = plt.figure(figsize=(5,7))
        for city in range(self.n_cities):
            city_plot[city+1] = plt.subplot(3,2,city+1)
            group_plot[city+1] = dict()
            group_plot[city+1]["S"] = city_plot[city+1].scatter(self.particles[city+1]["S"][:,0],self.particles[city+1]["S"][:,1], color = "blue")
            group_plot[city+1]["I"] = city_plot[city+1].scatter(self.particles[city+1]["I"][:,0],self.particles[city+1]["I"][:,1], color = "red")
            group_plot[city+1]["R"] = city_plot[city+1].scatter(self.particles[city+1]["R"][:,0],self.particles[city+1]["R"][:,1], color = "grey")
            
        self.fig = fig

        self.city_plot = city_plot
        self.group_plot = group_plot
        return self.fig, self.city_plot, self.group_plot
        
    
    def run_until_zero_infection(self, frame):
        for city in range(self.n_cities):
            self.particles[city+1] = dict()
            self.particles[city+1]["S"] = np.zeros([0,2])
            self.particles[city+1]["I"] = np.zeros([0,2])
            self.particles[city+1]["R"] = np.zeros([0,2])
            for n in range(len(self.boxes[city])):
                person = self.boxes[city][n]
                person.update_position()
                person.update_time()
                if person.status is "S":
                    self.particles[city+1]["S"] = np.append(self.particles[city+1]["S"],[person.point],0)
                if person.status is "I":
                    self.particles[city+1]["I"] = np.append(self.particles[city+1]["I"],[person.point],0)
                if person.status is "R":
                    self.particles[city+1]["R"] = np.append(self.particles[city+1]["R"],[person.point],0)
        
            print("status[S, I, R]", self.get_status_count())
            self.update_statuses()
        
        '''
        # If number of infected people is zero, stop the simulation.
        elif self.get_status_count()[1] == 0:
            break
        '''
        
        for city in range(self.n_cities):
            self.group_plot[city+1]["S"].set_offsets(self.particles[city+1]["S"])
            self.group_plot[city+1]["I"].set_offsets(self.particles[city+1]["I"])
            self.group_plot[city+1]["R"].set_offsets(self.particles[city+1]["R"])
        
        return self.group_plot,
    
    
    def show(self):
        plt.show()
    
    
    def add_R_label(self):
        update_period = 1
        all_R0_values = []
        
        if (self.time - self.last_update_time) > update_period:
            return
        self.last_update_time = self.time

        for city in self.boxes:
            for person in city:
                if person.status == "I":
                    prop = (person.time - person.infection_start_time) / person.infection_duration
                    # When time passes beyond the interval of 0.1, update the "reproduction number."
                    if prop > 0.1:
                        all_R0_values.append(person.num_infected / prop)
        
        if len(all_R0_values) > 0:
            # The "effect_reproduction_num"  is for R_label.
            self.effect_reproduction_num = np.mean(all_R0_values)

sim = RunSimpleSimulation()
sim.show()