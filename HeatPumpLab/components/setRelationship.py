"""Simple algebraic relationships used during model setup."""

from .stream import Stream


class SetRelationship:
    """Apply linear relationships between source and target stream properties."""

    name: str
    source_obj: list
    target_obj: list
    converged: bool = False

    def __init__(
        self,
        name="",
        source_obj="",
        target_obj="",
        m=1,
        c=0,
        prop="temperature",
        prop2=0,
    ):
        """Configure the relationship and immediately solve it."""
        self.name = name
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.m = m
        self.c = c
        self.prop = prop
        self.prop2 = prop2
        self.solve()

    def solve(self):
        """Apply the selected relationship type."""
        if self.prop == "temperature":
            self.set_temperature()

        elif self.prop == "quality":
            self.set_quality()

        elif self.prop == "Mass Flow":
            self.set_mass_flow()

        self.converged = True

    def set_temperature(self):
        """Propagate a linear temperature relation to the targets."""
        obj: Stream
        for obj in self.target_obj:
            obj.temperature = self.m * self.source_obj.temperature + self.c
            if self.prop2 == "pressure":
                obj.props_in_mPT()
            elif self.prop2 == "quality":
                obj.props_in_mTQ()
            else:
                raise ValueError("Set relationship failed: " + self.name)

    def set_quality(self):
        """Propagate a quality relation to the targets."""
        obj: Stream
        for obj in self.target_obj:
            obj.quality = self.source_obj[0]
            if self.prop2 == "pressure":
                obj.props_in_mPQ()
            elif self.prop2 == "temperature":
                obj.props_in_mTQ()
            else:
                raise ValueError("Set relationship failed: " + self.name)

    def set_mass_flow(self):
        """Scale target mass flows based on the source mass flow."""
        obj: Stream
        for obj in self.target_obj:
            obj.mass_flow = self.m * self.source_obj.mass_flow + self.c
            obj.heat_flow = obj.mass_flow * obj.enthalpy
