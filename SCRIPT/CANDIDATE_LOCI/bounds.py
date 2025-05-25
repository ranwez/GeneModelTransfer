from attrs import define
@define(slots=True)
class Bounds:
    start:int
    end:int # inclusive
    def __init__(self, start:int, end:int):
        self.start = min(start,end)
        self.end = max(start,end)
    def length(self) -> int:
        return self.end - self.start + 1
    def overlap(self, other: "Bounds") -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))
    @classmethod
    def clone(cls, other: "Bounds") -> "Bounds":
        """Copy constructor to create a new instance from another Bounds object."""
        return cls(other.start, other.end)