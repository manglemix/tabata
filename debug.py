from particulates import *
from time import time


p1 = BondingParticle(origin=vector(0, sqrt(2)), velocity=vector(10, 1), name="p1")
p2 = BondingParticle(origin=vector(20, 0), velocity=vector(-10, 1), name="p2")
s = EventSystem([p1, p2])
st = time()
print(s.get_total_momentum())
s.process(max_collisions=2)
print("time taken:", round(time() - st, 6))
print(s.get_total_momentum())
print(p1.origin, p1.velocity)
print(p2.origin, p2.velocity)