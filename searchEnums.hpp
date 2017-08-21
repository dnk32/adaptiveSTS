#ifndef SEARCH_ENUMS_HPP
#define SEARCH_ENUMS_HPP

enum searchType{
    FREE_END_TIME = 1,          // list of fixed start times provided
    FREE_START_TIME = -1        // list of fixed end times provided
};

enum pathCostType{
    MIN_DETERMINISTIC = 0,
    MIN_EXPECTED = 1
};

#endif
