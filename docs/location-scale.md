# Location-Scale Transformations
Location-Scale transformations are supported using the `location-scale` function. 

`(location-scale distribution location scale)`

Specifically, if a distribution d is assumed for a random variable Z and one applies the transformations
Y=a+bZ, then the distribution of Y is obtained via `(location-scale d a b)`. As an illustration, lets consider the actual implementation of the non-standard students-t distribution.
