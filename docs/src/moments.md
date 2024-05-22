# Moments and expectation

## Moment

Given a measure ``\mu`` and a monomial ``m``, the moment ``m`` of the measure is defined by the expectation ``\mathbb{E}_\mu[m]``.
Given a monomial and a value for the moment, a moment can be created using the `moment` function
```@docs
moment
```
The `moment` function returns an `AbstractMoment` which is a subtype of `AbstractMomentLike`.
An `AbstractMomentLike` is a type that can act like an `AbstractMoment` (it is similar to MultivariatePolynomials' `AbstractMonomialLike`, `AbstractTermLike` and `AbstractPolynomialLike`),
that is, it implements the following two functions
```@docs
moment_value
```

## Measure

Given a monomials and a values for the moments, a "moment vector" can be created using the `moment_vector` function
```@docs
moment_vector
```
The `moment_vector` function returns an `AbstractMeasure` which is a subtype of `AbstractMeasureLike`.
Note that it does not actually compute the probability density function of a measure having these moments, it simply stores a vector of moments belonging to a hypothetical measure.
However, it acts like a measure when taking its scalar product with a polynomial.

An `AbstractMeasureLike` is a type that can act like an `AbstractMeasure`,
that is, it implements the following two functions
```@docs
MultivariatePolynomials.variables(::MultivariateMoments.MomentVector)
MultivariatePolynomials.maxdegree(::MultivariateMoments.MomentVector)
MultivariatePolynomials.mindegree(::MultivariateMoments.MomentVector)
MultivariatePolynomials.extdegree(::MultivariateMoments.MomentVector)
moments
```

The moments of the dirac measure for a vector of monomials can be obtained by the `dirac` function
```@docs
dirac
```

## Expectation

The expectation of polynomial with respect to a measure can be computed either using `MultivariateMoments.expectation` or using the `Base.dot` scalar product.
```@docs
MultivariateMoments.expectation
MultivariateMoments.dot
```
