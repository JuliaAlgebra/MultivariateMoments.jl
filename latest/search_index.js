var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MultivariateMoments-–-Multivariate-Moments-for-Julia-1",
    "page": "Index",
    "title": "MultivariateMoments –- Multivariate Moments for Julia",
    "category": "section",
    "text": "Extension of MultivariatePolynomials to moments of multivariate measures and their scalar product with polynomials. It also includes the extraction of atomic measures from moment matrices which is related to the moment problem."
},

{
    "location": "index.html#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"moments.md\", \"atom.md\"]\nDepth = 2"
},

{
    "location": "moments.html#",
    "page": "Moments and expectation",
    "title": "Moments and expectation",
    "category": "page",
    "text": ""
},

{
    "location": "moments.html#Moments-and-expectation-1",
    "page": "Moments and expectation",
    "title": "Moments and expectation",
    "category": "section",
    "text": ""
},

{
    "location": "moments.html#MultivariateMoments.moment",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.moment",
    "category": "function",
    "text": "moment(α, m::AbstractMonomial)\n\nCreates the moment of the monomial m of value α.\n\n\n\n"
},

{
    "location": "moments.html#MultivariateMoments.value",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.value",
    "category": "function",
    "text": "value(m::AbstractMomentLike)\n\nReturns the value of the moment m.\n\nExamples\n\nCalling value(moment(3.1, x*y^2)) should return 3.1.\n\n\n\n"
},

{
    "location": "moments.html#MultivariatePolynomials.monomial-Tuple{MultivariateMoments.Moment}",
    "page": "Moments and expectation",
    "title": "MultivariatePolynomials.monomial",
    "category": "method",
    "text": "monomial(m::AbstractMomentLike)\n\nReturns the monomial of the moment m.\n\nExamples\n\nCalling monomial(moment(3.1, x*y^2)) should return x*y^2.\n\n\n\n"
},

{
    "location": "moments.html#Moment-1",
    "page": "Moments and expectation",
    "title": "Moment",
    "category": "section",
    "text": "Given a measure mu and a monomial m, the moment m of the measure is defined by the expectation mathbbE_mum. Given a monomial and a value for the moment, a moment can be created using the moment functionmomentThe moment function returns an AbstractMoment which is a subtype of AbstractMomentLike. An AbstractMomentLike is a type that can act like an AbstractMoment (it is similar to MultivariatePolynomials\' AbstractMonomialLike, AbstractTermLike and AbstractPolynomialLike), that is, it implements the following two functionsvalue\nmonomial(::MultivariateMoments.Moment)"
},

{
    "location": "moments.html#MultivariateMoments.measure",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.measure",
    "category": "function",
    "text": "measure(a, X::AbstractVector{<:AbstractMonomial})\n\nCreates a measure with moments moment(a[i], X[i]) for each i.\n\n\n\n"
},

{
    "location": "moments.html#MultivariatePolynomials.variables-Tuple{MultivariateMoments.Measure}",
    "page": "Moments and expectation",
    "title": "MultivariatePolynomials.variables",
    "category": "method",
    "text": "variables(μ::AbstractMeasureLike)\n\nReturns the variables of μ in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.\n\n\n\n"
},

{
    "location": "moments.html#MultivariatePolynomials.monomials-Tuple{MultivariateMoments.Measure}",
    "page": "Moments and expectation",
    "title": "MultivariatePolynomials.monomials",
    "category": "method",
    "text": "monomials(μ::AbstractMeasureLike)\n\nReturns an iterator over the monomials of μ sorted in the decreasing order.\n\n\n\n"
},

{
    "location": "moments.html#MultivariateMoments.moments",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.moments",
    "category": "function",
    "text": "moments(μ::AbstractMeasureLike)\n\nReturns an iterator over the moments of μ sorted in decreasing order of monomial.\n\n\n\n"
},

{
    "location": "moments.html#MultivariateMoments.dirac",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.dirac",
    "category": "function",
    "text": "dirac(X::AbstractVector{<:AbstractMoment}, s::AbstractSubstitution...)\n\nCreates the dirac measure by evaluating the moments of X using s.\n\nExamples\n\nCalling dirac([x*y, x*y^2], x=>3, y=>2) should the measure with moment x*y of value 6 and moment x*y^2 of value 12.\n\n\n\n"
},

{
    "location": "moments.html#Measure-1",
    "page": "Moments and expectation",
    "title": "Measure",
    "category": "section",
    "text": "Given a monomials and a values for the moments, a \"measure\" can be created using the measure functionmeasureThe measure function returns an AbstractMeasure which is a subtype of AbstractMeasureLike. Note that it does not actually compute the probability density function of a measure having these moments, it simply stores a vector of moments belonging to a hypothetical measure. However, it acts like a measure when taking its scalar product with a polynomial.An AbstractMeasureLike is a type that can act like an AbstractMeasure, that is, it implements the following two functionsvariables(::MultivariateMoments.Measure)\nmonomials(::MultivariateMoments.Measure)\nmomentsThe moments of the dirac measure for a vector of monomials can be obtained by the dirac functiondirac"
},

{
    "location": "moments.html#MultivariateMoments.expectation",
    "page": "Moments and expectation",
    "title": "MultivariateMoments.expectation",
    "category": "function",
    "text": "MultivariateMoments.expectation(μ::AbstractMeasureLike, p::AbstractPolynomialLike)\nMultivariateMoments.expectation(p::AbstractPolynomialLike, μ::AbstractMeasureLike)\n\nComputes the expectation mathbbE_mup.\n\n\n\n"
},

{
    "location": "moments.html#Base.LinAlg.dot",
    "page": "Moments and expectation",
    "title": "Base.LinAlg.dot",
    "category": "function",
    "text": "dot(μ::AbstractMeasureLike, p::AbstractPolynomialLike)\ndot(p::AbstractPolynomialLike, μ::AbstractMeasureLike)\n\nSee expectation\n\n\n\n"
},

{
    "location": "moments.html#Expectation-1",
    "page": "Moments and expectation",
    "title": "Expectation",
    "category": "section",
    "text": "The expectation of polynomial with respect to a measure can be computed either using MultivariateMoments.expectation or using the Base.dot scalar product.MultivariateMoments.expectation\ndot"
},

]}
