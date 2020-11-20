from particulates import *
from math import ceil, floor


class Atom(Particle):
    """Defines a particle which bond with other atoms to complete their outer shell"""
    def __init__(self, valence_electrons: int, stability: float = 0.5, electron_impulse: float = 0.1, radius: float = 1, mass: float = 1, origin: TYPE_ARRAY = np.zeros(2, float), velocity: TYPE_ARRAY = np.zeros(2, float), id=None):
        super().__init__(radius, mass, origin, velocity, id)
        self.initial_valence_electrons = valence_electrons
        self.current_valence_electrons = valence_electrons
        self.stability = stability
        self.electron_impusle = electron_impulse
        self._bonds = []
        self._stability_time = 0

    def process(self, t: float) -> None:
        """Calculates the stability of the particle across time"""
        self.origin += self.velocity * t
        self._stability_time -= t

        if self._stability_time == 0:
            for bond in self._bonds:
                pass
                # TODO Finish

    def handle_collision(self, other: Particle) -> None:
        if type(other) == Atom:
            shared = other.current_valence_electrons + self.current_valence_electrons

            if shared <= 8:
                self._bonds.append(other)
                self.inelastic_collide(other)
                shared /= 2
                if self.current_valence_electrons > other.current_valence_electrons:
                    self.current_valence_electrons = ceil(shared)
                    other.current_valence_electrons = floor(shared)

                else:
                    self.current_valence_electrons = floor(shared)
                    other.current_valence_electrons = ceil(shared)

                self._stability_time = self.stability * (self.current_valence_electrons - self.initial_valence_electrons) / (8 - self.initial_valence_electrons)

        else:
            self.elastic_collide(other)

    def next_event(self, system) -> float:
        """Alerts the system to """
        return self._stability_time
