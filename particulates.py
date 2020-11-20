from typing import Union
from math import sqrt, nan, isnan, inf, ceil, floor
import numpy as np
from numpy import dot
from numpy.linalg import norm as length
from time import time

TYPE_ARRAY = Union[list, tuple]


def quad_formula(a: float, b: float, c: float, neg=False) -> float:
    """
    Returns the real solution to the quadratic equation formed by the coefficients given
    Returns nan if solution is complex, or a is zero
    Set neg to True to calculate the negative solution
    """
    if a == 0:
        return nan

    determinant = b * b - 4 * a * c
    if determinant < 0:
        return nan

    return (- b + (neg * - 2 + 1) * sqrt(determinant)) / 2 / a


def vector(*args) -> np.ndarray:
    """
    Constructs a numpy vector of a variable dimensions
    Also ensures that the data type of the vector is a float

    eg. a = vector(1, 2)
        b = vector(2, 5, 1)

        print(a)    # [1. 2.]
        print(b)    # [2. 5. 1.]
    """
    return np.array(args, float)


def normalized(vector: np.ndarray) -> np.ndarray:
    """Returns the unit/normalized vector of the vector given"""
    return vector / length(vector)


def generate_pairs(array: TYPE_ARRAY) -> tuple:
    """A generator which yields pairs of objects from the array, such that all pairs are unique"""
    array_len = len(array)
    for i, obj1 in enumerate(array):
        for j in range(i, array_len):
            yield obj1, array[j]


class Particle:
    def __init__(self, radius=1, mass=1, origin=np.zeros(2, float), velocity=np.zeros(2, float), name=None):
        """
        A particle shaped like a circle or sphere
        Can participate in Particle-Particle collisions

        :param radius: the distance from the origin to the collidable surface
        :param mass: affects the amount of momentum the particle has
        :param origin: the centre of the particle
        :param velocity: the speed and magnitude of the particle

        origin and velocity are both vectors of any dimension (as long as they are in same dimension)
        origin and velocity can both be lists, tuples, or ndarrays (vectors)
        """
        self.radius = radius
        self.mass = mass
        self.name = None

        if type(origin) != np.ndarray:
            self.origin = vector(*origin)
        else:
            self.origin = origin

        if type(origin) != np.ndarray:
            self.velocity = vector(*velocity)
        else:
            self.velocity = velocity

    def check_collision(self, other: "Particle") -> bool:
        """Returns True if this Particle is intersecting another particle"""
        return length(self.origin - other.origin) <= self.radius + other.radius

    def predict_collision(self, other: "Particle") -> float:
        """
        Returns the time in the future when self collides with the Particle given
        Returns nan if there is no collision
        """
        v = other.velocity - self.velocity
        o = other.origin - self.origin
        return quad_formula(length(v) ** 2, 2 * dot(v, o), length(o) ** 2 - (self.radius + other.radius) ** 2, True)

    def elastic_collide(self, other: "Particle") -> None:
        """Calculates the new velocities of self and other if they were to elastically collide"""
        momentum_vector = (self.velocity - other.velocity) * self.mass
        normal_vector = other.origin - self.origin
        impulse_vector = dot(momentum_vector, normalized(normal_vector)) * normalized(normal_vector)
        other.velocity += impulse_vector / other.mass
        self.velocity -= impulse_vector / self.mass

    def inelastic_collide(self, other: "Particle") -> None:
        """Calculates the new velocities of self and other if they were to inelastically collide"""
        self.velocity = other.velocity = (self.velocity * self.mass + other.velocity * other.mass) / (self.mass + other.mass)

    def handle_collision(self, other: "Particle") -> None:
        """An overridable function to determine the collision that should occur with another particle"""
        self.elastic_collide(other)

    def next_event(self, system: "System") -> float:
        """This function is called by the system to check if this particle has any events happening soon"""
        return inf


class BondingParticle(Particle):
    """A particle which can bond with another Bonding Particle"""
    def handle_collision(self, other: "Particle") -> None:
        if type(other) == BondingParticle:
            self.inelastic_collide(other)

        else:
            self.elastic_collide(other)


class System:
    def __init__(self, particles: TYPE_ARRAY, bounds=((-100, 100), (-100, 100))):
        """
        Defines a space in which Particles can move around

        :param bounds: The coordinates of the walls in each axes
        Warning! Particles spawned outside of bounds will cause time reversals
        """
        self.particles = particles
        self.bounds = bounds
        self._collisions = 0
        self.time = 0

    def check_system(self) -> dict:
        """Returns a report of any errors in the system"""
        out_of_bounds = []      # a list of particles outside of the bounds
        intersecting = []       # a list of pairs of particles which are intersecting

        for particle, other in generate_pairs(self.particles):
            if length(other.origin - particle.origin) < other.radius + particle.radius:
                intersecting.append((particle, other))

        for particle in self.particles:
            for i, bound in enumerate(self.bounds):
                if particle.origin[i] > max(bound) - particle.radius or particle.origin[i] < min(bound) + particle.radius:
                    out_of_bounds.append(particle)
                    break

        return {"out_of_bounds": out_of_bounds, "intersecting": intersecting}

    def get_system_statistics(self, params=("origin", "velocity")) -> dict:
        """Returns a dictionary of the params of each particle in the system"""
        return {particle.name: {param: getattr(particle, param) for param in params} for particle in self.particles}

    def process(self, max_time=100, max_collisions=10) -> None:
        """
        Runs the system until the simulation time exceeds max_time, or more collisions than max_collisions have occurred
        """
        t = 0
        for collision in range(max_collisions):
            next_t = []         # stores the time of occurrence of the next events
            collider1 = []      # stores the first particle involved in an event
            collider2 = []      # stores the second particle or the index of the bounds involved in an event

            # check for particle-particle collisions
            for particle, other in generate_pairs(self.particles):
                prediction = particle.predict_collision(other)
                if not isnan(prediction) and prediction > 0:
                    next_t.append(prediction)
                    collider1.append(particle)
                    collider2.append(other)

            # check for particle-bounds collisions
            for particle in self.particles:
                for i, bound in enumerate(self.bounds):
                    speed = particle.velocity[i]
                    if speed == 0:
                        continue

                    elif speed > 0:
                        next_t.append((max(bound) - particle.radius - particle.origin[i]) / speed)

                    else:
                        next_t.append((min(bound) + particle.radius - particle.origin[i]) / speed)

                    collider1.append(particle)
                    collider2.append(i)

            if len(next_t) > 0:
                min_t = min(next_t)

            else:
                min_t = 0

            if len(next_t) == 0 or t + min_t > max_time:
                finish_t = self.time + min_t - max_time
                for particle in self.particles:
                    particle.origin += particle.velocity * finish_t

                self.time += max_time
                return

            for particle in self.particles:
                particle.origin += particle.velocity * min_t

            t += min_t
            min_index = next_t.index(min_t)
            if type(collider2[min_index]) == int:
                collider1[min_index].velocity[collider2[min_index]] *= -1

            else:
                collider1[min_index].handle_collision(collider2[min_index])

            self._collisions += 1
        self.time += t

    def get_total_momentum(self) -> np.ndarray:
        """Returns the sum of the momentum's of the particles in the system"""
        return sum([p.velocity * p.mass for p in self.particles])

    def get_total_collisions(self) -> int:
        """Returns the total number of collisions which have occurred"""
        return self._collisions
