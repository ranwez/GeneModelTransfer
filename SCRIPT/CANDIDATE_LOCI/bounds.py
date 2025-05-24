from attrs import define, field

@define(slots=True)
class Bounds:
    start: int
    end: int

    def __init__(self, start, end):
        self.start=min(start,end)
        self.end=max(start,end)

    def length(self) -> int:
        return self.end - self.start + 1

    def overlap(self, other: "Bounds") -> int:
        if self.end < other.start or other.end < self.start:
            return 0
        return min(self.end, other.end) - max(self.start, other.start) + 1

    def distance(self, other: "Bounds") -> int:
        """Calculate the distance between two bounds.
        
        Returns 0 if the bounds overlap, otherwise returns the minimum distance
        between their extremities.
        
        Args:
            other: Another Bounds object to compare with
        Returns:
            The distance between the two bounds
        """
        if self.overlap(other) > 0:
            return 0
        return max(self.start - other.end, other.start - self.end)

    @classmethod
    def clone(cls, bounds: "Bounds") -> "Bounds":
        return cls(bounds.start, bounds.end)

@define(slots=True)
class RangeCoverage:
    """Track the number of covered elements for a set of intervals.
    
    This class maintains a cumulative count of covered elements up to each intervalary position.
    Used to compute the total coverage between any two positions in the intervals set.
    """
    coverage_up_to: dict[int, int] = field(factory=dict)
    
    def __init__(self, intervals_list: list[Bounds]):
        """Initialize coverage tracking from a list of intervals.
        
        Args:
            intervals_list: List of bounds (will be sorted by start position)
        Raises:
            ValueError: If list is empty
        """
        if len(intervals_list) == 0:
            raise ValueError("intervals_list should not be empty")
        self.coverage_up_to = {}

        sorted_intervals = sorted(intervals_list, key=lambda b: b.start)
        (current_range_start, current_range_end) = (sorted_intervals[0].start, sorted_intervals[0].end)
        prev_loci_coverage = 0
        
        for interval in sorted_intervals:
            if interval.start > current_range_end:
                prev_loci_coverage = self.coverage_up_to[current_range_end]
                (current_range_start, current_range_end) =(interval.start, interval.end)
            current_range_end = max(current_range_end, interval.end)
            self.coverage_up_to[interval.start] = prev_loci_coverage + (interval.start - current_range_start + 1)    
            self.coverage_up_to[interval.end] = prev_loci_coverage + (interval.end - current_range_start + 1)
    
    def get_coverage_up_to(self, pos: int) -> int:
        """Get the number of covered elements up to a specific position.
        
        Args:
            pos: Position to check coverage
        Returns:
            Number of elements covered up to the specified position
        Raises:
            ValueError: If position is not found in coverage map
        """
        if pos not in self.coverage_up_to:
            raise ValueError(f"Position {pos} not found in cumulative_coverage")
        return self.coverage_up_to[pos]
    
    def get_coverage(self, start: int, end: int) -> int:
        """Get the number of covered elements between start and end positions.
        
        Args:
            start: Start position to check coverage
            end: End position to check coverage
        Returns:
            Number of elements covered between start and end
        Raises:
            ValueError: If start or end positions are not found in coverage map
        """
        if start not in self.coverage_up_to or end not in self.coverage_up_to:
            raise ValueError(f"Start {start} or end {end} not found in cumulative_coverage")
        if start > end:
            return 0
        return self.coverage_up_to[end] - self.coverage_up_to[start] + 1