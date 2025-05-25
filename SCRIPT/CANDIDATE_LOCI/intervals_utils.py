
## This class represents a list of intervals and implements multiple basic
# operations between intervals. It is used to represent and compare lists of
# CDS bounds of compared loci in the mRNA comparison function 
# compute_matches_mismatches_EI_RF
#
# @see compute_matches_mismatches_EI_RF()
#
# @remark This class was adapted from https://stackoverflow.com/a/20062829 
# by Vincent Ranwez, and translated into python using chatGPT
class OrderedIntervals:
    
    ## Initialises the class instance with the given intervals
    #
    # @param intervals List of interval bounds to represent
    #
    # @param include_ub Boolean indicating wether to include the upper bounds
    # of each interval of the list
    #
    # @remark for CDS bounds of GFF annotations, include_ub needs to be 'True'
    def __init__(self, intervals, include_ub=False):
        if include_ub:
            self.intervals = OrderedIntervals.transform_intervals_to_exclude_ub(intervals)
        else:
            ### List of sequence coordinates representing CDS intervals
            self.intervals = intervals
    
    ## Excludes the upper bound of each interval of the given list
    #
    # @param intervals Interval list to be transformed
    #
    # @returns Returns the same list of intervals with each upper bound
    # removed
    @staticmethod
    def transform_intervals_to_exclude_ub(intervals):
        transformed_intervals = []
        for i in range(0, len(intervals), 2):
            lower_bound = intervals[i]
            upper_bound = intervals[i + 1] + 1
            transformed_intervals.extend([lower_bound, upper_bound])
        return(transformed_intervals)

   ## Returns the class instance's intervals list including the upper bound of
   # each interval
   #
   # @returns Returns the list of intervals of the class instance with
   # their upper bounds
    def get_intervals_with_included_ub(self):
        transformed_intervals = []
        for i in range(0, len(self.intervals), 2):
            lower_bound = self.intervals[i]
            upper_bound = self.intervals[i + 1] - 1
            transformed_intervals.extend([lower_bound, upper_bound])
        return(transformed_intervals)

    ## Calculates the total length of all intervals.
    #
    # @returns The sum of the lengths of all intervals.
    def total_length(self):
        return sum(self.intervals[i + 1] - self.intervals[i] for i in range(0, len(self.intervals), 2))

    ## Creates a new OrderedIntervals object.
    #
    # @param intervals List of integers representing the intervals.
    #
    # @param include_ub Boolean indicating if the upper bound should be 
    # included.
    #
    # @returns A new OrderedIntervals object.
    @staticmethod
    def new(intervals, include_ub=False):
        return OrderedIntervals(intervals,include_ub)
    
    ## Calculates the union an interval list with another interval list.
    #
    # @param other An OrderedIntervals object.
    #
    # @returns A new OrderedIntervals object representing the union of the two 
    # interval lists.
    def union(self, other):
        return self.merge(other, lambda a, b: a or b)


    ## Calculates the intersection of this interval list with another interval 
    # list.
    #
    # @param other An OrderedIntervals object.
    #
    # @returns A new OrderedIntervals object representing the intersection of 
    # the two interval lists.
    def intersection(self, other):
        return self.merge(other, lambda a, b: a and b)

    ## Calculates the difference of this interval list with another interval 
    # list and returns the intervals of this list (asymmetric difference)
    #
    # @param other Another OrderedIntervals object.
    #
    # @returns A new OrderedIntervals object representing the difference 
    # between the two interval lists.
    #
    # @remark This method returns the intervals of the first interval list not
    # present in the second, but not those of the second not present in the 
    # first. To get both, use the symmetric_difference method
    #
    # @see symmetric_difference()
    def difference(self, other):
        return self.merge(other, lambda a, b: a and not b)

    ## Calculates the difference of this interval list with another interval 
    # list and returns both intervals list (symmetric difference)
    #
    # @param other Another OrderedIntervals object.
    #
    # @returns A new OrderedIntervals object representing the difference 
    # between the two interval lists.
    #
    # @remark This method returns the intervals of the first interval list not
    # present in the second, and those of the second not present in the 
    # first. To get only the difference for this instance's list, use the 
    # difference method
    #
    # @see difference()
    def symmetric_difference(self, other):
        return self.merge(other, lambda a, b: a ^ b)

    ## Merges this interval list with another interval list based on the 
    # specified operator.
    #
    # @param other Another OrderedIntervals object.
    #
    # @param keep_operator A function that determines whether to keep an 
    # interval based on its presence in either list.
    #
    # @returns A new OrderedIntervals object representing the merged intervals.
    def merge(self, other, keep_operator):
        res = []
        if not self.intervals and not other.intervals:
            return OrderedIntervals(res)

        sentinel = max(self.intervals[-1], other.intervals[-1]) + 1

        # create iterators for each intervals list
        interval_iters = (
            iter(self.intervals + [sentinel]),
            iter(other.intervals + [sentinel]),
        )

        # initialise all values for the first intervals bounds
        current_bound = (
            next(interval_iters[0]),
            next(interval_iters[1]),
        )
        scan = min(current_bound[0], current_bound[1])
        current_is_lb = (True, True)
        next_res_is_lb = True

        # for each interval bounds, determine wether we are in an interval for
        # both lists, use 'keep_operator' to determine wether to keep the 
        # interval, and if yes append to 'res'
        while scan < sentinel:
            in_interval = (
                (scan >= current_bound[0]) == current_is_lb[0],
                (scan >= current_bound[1]) == current_is_lb[1],
            )
            in_res = keep_operator(in_interval[0], in_interval[1])

            if in_res == next_res_is_lb:
                res.append(scan)
                next_res_is_lb = not next_res_is_lb

            if scan == current_bound[0]:
                current_bound = (next(interval_iters[0]), current_bound[1])
                current_is_lb = (not current_is_lb[0], current_is_lb[1])

            if scan == current_bound[1]:
                current_bound = (current_bound[0], next(interval_iters[1]))
                current_is_lb = (current_is_lb[0], not current_is_lb[1])

            scan = min(current_bound[0], current_bound[1])

        return OrderedIntervals(res)
