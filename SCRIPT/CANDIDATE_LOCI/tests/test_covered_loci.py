import pytest
from CANDIDATE_LOCI.bounds import RangeCoverage, Bounds

def create_test_interval(start: int, end: int) -> Bounds:
    """Helper function to create test intervals"""
    return Bounds(start, end)

def test_init_with_valid_intervals():
    """Test initialization with valid intervals"""
    intervals = [
        create_test_interval(10, 20),
        create_test_interval(30, 40),
        create_test_interval(50, 60)
    ]
    coverage = RangeCoverage(intervals)
    assert coverage.get_coverage_up_to(10) == 1  # First position
    assert coverage.get_coverage_up_to(20) == 11  # End of first interval
    assert coverage.get_coverage_up_to(30) == 12  # Start of second interval
    assert coverage.get_coverage_up_to(60) == 33  # End of last interval

def test_init_with_empty_list():
    """Test initialization with empty list raises ValueError"""
    with pytest.raises(ValueError, match="intervals_list should not be empty"):
        RangeCoverage([])

def test_get_coverage_valid_positions():
    """Test get_coverage with valid start and end positions"""
    intervals = [
        create_test_interval(10, 20),
        create_test_interval(30, 40)
    ]
    coverage = RangeCoverage(intervals)
    assert coverage.get_coverage(10, 20) == 11  # Coverage within first interval
    assert coverage.get_coverage(30, 40) == 11  # Coverage within second interval
    assert coverage.get_coverage(10, 40) == 22  # Coverage across intervals

def test_get_coverage_invalid_positions():
    """Test get_coverage with invalid positions raises ValueError"""
    coverage = RangeCoverage([create_test_interval(10, 20)])
    with pytest.raises(ValueError):
        coverage.get_coverage(5, 15)  # Start position not found
    with pytest.raises(ValueError):
        coverage.get_coverage(10, 25)  # End position not found

def test_overlapping_intervals():
    """Test initialization with overlapping intervals"""
    intervals = [
        create_test_interval(10, 25),
        create_test_interval(20, 30)
    ]
    coverage = RangeCoverage(intervals)
    assert coverage.get_coverage(10, 30) == 21  # Should handle overlap correctly

def test_complex_intervals():
    """Test initialization with complex overlapping intervals"""
    intervals = [
        create_test_interval(10, 25),
        create_test_interval(20, 30),
        create_test_interval(25, 28),
        create_test_interval(40, 60),
        create_test_interval(50, 55)
    ]
    coverage = RangeCoverage(intervals)
    assert coverage.get_coverage(25, 30) == 6  # Should handle overlap correctly
    assert coverage.get_coverage(10, 55) == 37  # Should handle multiple overlaps correctly

def test_adjacent_intervals():
    """Test initialization with adjacent intervals"""
    intervals = [
        create_test_interval(10, 20),
        create_test_interval(21, 30)
    ]
    coverage = RangeCoverage(intervals)
    assert coverage.get_coverage(10, 30) == 21  # Should handle adjacent intervals correctly