import numpy as np
import random
from space_ops import *
import matplotlib

COLOR_MAP = {
"S": 'blue',
"I": 'red',
"R": 'grey',
}



class Person(object):
    CONFIG = {
    "infection_radius": 0.5,
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
    "social_distance_factor": 0.2,
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
        self.social_distance_factor = self.CONFIG["social_distance_factor"]
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
        
        '''
        # Avoid walls
        wall_force = np.zeros(2)
        for i in range(2):
            # When i=0, to_lower and to_upper mean the horizontal distances
            # When i=1, to_lower and to_upper mean the vertical distances
            to_lower = point[i] - self.dl_bound[i]
            to_upper = point[i] - self.ur_bound[i]
            
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
        '''
        '''
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
        '''
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


# Example of usage
person = Person()
print("initial time: ", person.time, "inital point: ", person.point)
print("initial status: ", person.status)
print("initial infection_radius: ", person.infection_radius)

print("\n")
person.update_time()
person.update_position()
print("new time: ", person.time, "new point: ", person.point)

person.update_status("I")
print("new status", person.status)

person.update_infection_ring(1)
print("new infection_radius: ", person.infection_radius)
