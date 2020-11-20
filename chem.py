from particulates import *


class Atom(BondingParticle):
    """Defines a particle which bond with other atoms to complete their outer shell"""
    def __init__(self, valence_electrons, radius=1, mass=1, origin=np.zeros(2, float), velocity=np.zeros(2, float), stability= 0.5, name=None):
        super().__init__(radius, mass, origin, velocity, name)
        self.initial_valence_electrons = valence_electrons
        self.current_valence_electrons = valence_electrons
        self.stability = stability

    def handle_collision(self, other: "Particle") -> None:
        if type(other) == Atom:
            shared = other.current_valence_electrons + self.current_valence_electrons

            if shared <= 8:
                self.inelastic_collide(other)
                shared /= 2
                if self.current_valence_electrons > other.current_valence_electrons:
                    self.current_valence_electrons = ceil(shared)
                    other.current_valence_electrons = floor(shared)

                else:
                    self.current_valence_electrons = floor(shared)
                    other.current_valence_electrons = ceil(shared)

        else:
            self.elastic_collide(other)

    def next_event(self, system: "System") -> float:
        """This function is called by the system to check if this particle has any events happening soon"""
        return self.stability * (self.current_valence_electrons - self.initial_valence_electrons) / (8 - self.initial_valence_electrons)
