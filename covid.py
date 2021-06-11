import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import random
from space_ops import *

class Person(object):
    CONFIG = {
    "infection_radius": 2,
    "probability_of_infection": 0.2,
    "counted_as_infector": False,
    "incubation_period": 1.0,
    "infection_duration": 5.0,
    "dl_bound": [0, 0],
    "ur_bound": [10, 10],
    "wander_step_size": 1,
    "wander_step_duration": 1,
    "gravity_well": None,
    "gravity_strength": 1,
    "switch_for_social_distancing": False,
    "percentage_of_social_distancing": 0.7,
    "social_distance_factor": 0.2,
    "n_repulsion_points":10,
    "wall_buffer": 1,
    "max_speed": 1,
    "dt": 0.01,
    }
    

    def __init__(self):
    
        self.time = 0
        self.last_step_change = -1
        self.velocity = np.zeros(2)
        self.status = "S"
        self.repulsion_points = []
        self.num_infected = 0
        self.getting_infected_time = 0
        
        self.infection_radius = self.CONFIG["infection_radius"]
        self.probability_of_infection = self.CONFIG["probability_of_infection"]
        self.counted_as_infector = self.CONFIG["counted_as_infector"]
        self.incubation_period = self.CONFIG["incubation_period"]
        self.infection_duration = self.CONFIG["infection_duration"]
        self.dl_bound = self.CONFIG["dl_bound"]
        self.ur_bound = self.CONFIG["ur_bound"]
        self.wander_step_size = self.CONFIG["wander_step_size"]
        self.wander_step_duration = self.CONFIG["wander_step_duration"]
        self.gravity_well = self.CONFIG["gravity_well"]
        self.gravity_strength = self.CONFIG["gravity_strength"]
        self.percentage_of_social_distancing = self.CONFIG["percentage_of_social_distancing"]
        self.social_distance_factor = self.CONFIG["social_distance_factor"]
        self.n_repulsion_points = self.CONFIG["n_repulsion_points"]
        self.wall_buffer = self.CONFIG["wall_buffer"]
        self.max_speed = self.CONFIG["max_speed"]
        self.dt = self.CONFIG["dt"]
        
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
        if self.CONFIG["switch_for_social_distancing"] is True:
            if self.social_distance_factor > 0 and random.random() < self.percentage_of_social_distancing:
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
	"city_pop":50, #the population in each city
	"box_size":10, #the size of each city(box)
	"travel_rate":0, #decide how much people in the city will travel
    "triggered_case": 15
    }

    def __init__(self):
        super().__init__()
        self.num_of_total_infected_case = 0

        self.num_of_toal_infected_case = 0
        self.n_cities = self.CONFIG["n_cities"]
        self.row_cities = self.CONFIG["row_cities"]
        self.city_pop=self.CONFIG["city_pop"]
        self.box_size=self.CONFIG["box_size"]
        self.travel_rate=self.CONFIG["travel_rate"]
        self.triggered_case = self.CONFIG["triggered_case"]

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
                city[n].append(person)
        for city in self.boxes:
            #choose the patient zero
            random.choice(city).update_status("I")
    

    def update_statuses(self):
        if self.num_of_total_infected_case > self.triggered_case:
            self.start_social_disdancing()
        
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
                    dist = np.linalg.norm(s_person.point - i_person.point)
                    if dist < s_person.infection_radius and random.random() < s_person.probability_of_infection:
                        #if s_person got an infection than update its status to I
                        s_person.update_status("I")
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
        for city in self.boxes:
            points = np.array([person.point for person in city])
            repelled_points = points
            for point, person in zip(points, city):
                if person.social_distance_factor > 0:
                    diffs = np.linalg.norm(repelled_points - point, axis=1)
                    person.repulsion_points = repelled_points[np.argsort(diffs)[1:person.n_repulsion_points + 1]]
    
    
    def start_social_distancing(self):
        super().CONFIG["switch_for_social_distancing"] = True



class RunSimpleSimulation(SIRSimulation):
    
    def __init__(self):
        super().__init__()
        self.last_update_time = 0
        self.effect_reproduction_num = 0
        
    
    
    def run_until_zero_infection(self):
        while True:
            self.update_statuses()

            for city in self.boxes:
                for person in city:
                    person.update_position()
                    person.update_time()
            
            # If number of infected people is zero, stop the simulation.
            if
    

    
    def run(self):
        animation = FuncAnimation(fig, run_until_zero_infection, interval = self.dt)
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
                    prop = (self.time - self.infection_start_time) / self.infection_duration
                    # When time passes beyond the interval of 0.1, update the "reproduction number."
                    if prop > 0.1:
                        all_R0_values.append(person.num_infected / prop)
        
        if len(all_R0_values) > 0:
            all_R0_values.append(np.mean(values))
            # The "effect_reproduction_num"  is for R_label.
            self.effect_reproduction_num = np.mean(all_R0_values)


    def get_status_counts(self):
        num_I_people = np.array([] )


















print("**Example for SIR_Simulation**")
sim = SIRSimulation()

for n in range(len(sim.boxes)):
    for i in range(len(sim.boxes[n])):
        person = sim.boxes[n][i]
        while person.time < 1:
            sim.update_statuses()
            person.update_position()
            person.update_time()
            
            
    print("city:",n+1)


for n in range(len(sim.boxes)):
    s_counter = 0
    for m in range(len(sim.boxes[n])):
        if sim.boxes[n][m].status is "S":
            s_counter += 1
    print("There's ",s_counter,"S_person in the ",n+1,"city\n")

for n in range(len(sim.boxes)):
    i_counter = 0
    for m in range(len(sim.boxes[n])):
        if sim.boxes[n][m].status is "I":
            i_counter += 1
    print("There's ",i_counter,"I_person in the ",n+1,"city\n")

for n in range(len(sim.boxes)):
    r_counter = 0
    for m in range(len(sim.boxes[n])):
        if sim.boxes[n][m].status is "R":
            r_counter += 1
    print("There's ",r_counter,"R_person in the ",n+1,"city")

'''
for n in range(len(sim.boxes)):
    for i in range(len(sim.boxes[n])):
        if sim.boxes[n][i].status is "I":
            print("The number ",i+1,"person in number ",n+1,"city is",sim.boxes[n][i].status)
'''



