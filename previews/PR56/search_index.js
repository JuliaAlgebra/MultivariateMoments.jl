var documenterSearchIndex = {"docs":
[{"location":"moments/#Moments-and-expectation","page":"Moments and expectation","title":"Moments and expectation","text":"","category":"section"},{"location":"moments/#Moment","page":"Moments and expectation","title":"Moment","text":"","category":"section"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"Given a measure mu and a monomial m, the moment m of the measure is defined by the expectation mathbbE_mum. Given a monomial and a value for the moment, a moment can be created using the moment function","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"moment","category":"page"},{"location":"moments/#MultivariateMoments.moment","page":"Moments and expectation","title":"MultivariateMoments.moment","text":"moment(α, m::AbstractMonomial)\n\nCreates the moment of the monomial m of value α.\n\n\n\n\n\n","category":"function"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"The moment function returns an AbstractMoment which is a subtype of AbstractMomentLike. An AbstractMomentLike is a type that can act like an AbstractMoment (it is similar to MultivariatePolynomials' AbstractMonomialLike, AbstractTermLike and AbstractPolynomialLike), that is, it implements the following two functions","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"moment_value\nMultivariatePolynomials.monomial(::MultivariateMoments.Moment)","category":"page"},{"location":"moments/#MultivariateMoments.moment_value","page":"Moments and expectation","title":"MultivariateMoments.moment_value","text":"moment_value(m::AbstractMomentLike)\n\nReturns the value of the moment m.\n\nExamples\n\nCalling moment_value(moment(3.1, x*y^2)) should return 3.1.\n\n\n\n\n\n","category":"function"},{"location":"moments/#MultivariatePolynomials.monomial-Tuple{Moment}","page":"Moments and expectation","title":"MultivariatePolynomials.monomial","text":"monomial(m::AbstractMomentLike)\n\nReturns the monomial of the moment m.\n\nExamples\n\nCalling monomial(moment(3.1, x*y^2)) should return x*y^2.\n\n\n\n\n\n","category":"method"},{"location":"moments/#Measure","page":"Moments and expectation","title":"Measure","text":"","category":"section"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"Given a monomials and a values for the moments, a \"measure\" can be created using the measure function","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"measure","category":"page"},{"location":"moments/#MultivariateMoments.measure","page":"Moments and expectation","title":"MultivariateMoments.measure","text":"measure(a::AbstractVector{T}, X::AbstractVector{<:AbstractMonomial}; rtol=Base.rtoldefault(T), atol=zero(T))\n\nCreates a measure with moments moment(a[i], X[i]) for each i. An error is thrown if there exists i and j such that X[i] == X[j] but !isapprox(a[i], a[j]; rtol=rtol, atol=atol).\n\n\n\n\n\n","category":"function"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"The measure function returns an AbstractMeasure which is a subtype of AbstractMeasureLike. Note that it does not actually compute the probability density function of a measure having these moments, it simply stores a vector of moments belonging to a hypothetical measure. However, it acts like a measure when taking its scalar product with a polynomial.","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"An AbstractMeasureLike is a type that can act like an AbstractMeasure, that is, it implements the following two functions","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"MultivariatePolynomials.variables(::MultivariateMoments.Measure)\nMultivariatePolynomials.monomials(::MultivariateMoments.Measure)\nmoments","category":"page"},{"location":"moments/#MultivariatePolynomials.variables-Tuple{Measure}","page":"Moments and expectation","title":"MultivariatePolynomials.variables","text":"variables(μ::AbstractMeasureLike)\n\nReturns the variables of μ in decreasing order. Just like in MultivariatePolynomials, it could contain variables of zero degree in every monomial.\n\n\n\n\n\n","category":"method"},{"location":"moments/#MultivariatePolynomials.monomials-Tuple{Measure}","page":"Moments and expectation","title":"MultivariatePolynomials.monomials","text":"monomials(μ::AbstractMeasureLike)\n\nReturns an iterator over the monomials of μ sorted in the decreasing order.\n\n\n\n\n\n","category":"method"},{"location":"moments/#MultivariateMoments.moments","page":"Moments and expectation","title":"MultivariateMoments.moments","text":"moments(μ::AbstractMeasureLike)\n\nReturns an iterator over the moments of μ sorted in decreasing order of monomial.\n\n\n\n\n\n","category":"function"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"The moments of the dirac measure for a vector of monomials can be obtained by the dirac function","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"dirac","category":"page"},{"location":"moments/#MultivariateMoments.dirac","page":"Moments and expectation","title":"MultivariateMoments.dirac","text":"dirac(X::AbstractVector{<:AbstractMoment}, s::AbstractSubstitution...)\n\nCreates the dirac measure by evaluating the moments of X using s.\n\nExamples\n\nCalling dirac([x*y, x*y^2], x=>3, y=>2) should the measure with moment x*y of value 6 and moment x*y^2 of value 12.\n\n\n\n\n\n","category":"function"},{"location":"moments/#Expectation","page":"Moments and expectation","title":"Expectation","text":"","category":"section"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"The expectation of polynomial with respect to a measure can be computed either using MultivariateMoments.expectation or using the Base.dot scalar product.","category":"page"},{"location":"moments/","page":"Moments and expectation","title":"Moments and expectation","text":"MultivariateMoments.expectation\nMultivariateMoments.dot","category":"page"},{"location":"moments/#MultivariateMoments.expectation","page":"Moments and expectation","title":"MultivariateMoments.expectation","text":"MultivariateMoments.expectation(μ::AbstractMeasureLike, p::AbstractPolynomialLike)\nMultivariateMoments.expectation(p::AbstractPolynomialLike, μ::AbstractMeasureLike)\n\nComputes the expectation mathbbE_mup.\n\n\n\n\n\n","category":"function"},{"location":"moments/#LinearAlgebra.dot","page":"Moments and expectation","title":"LinearAlgebra.dot","text":"dot(μ::AbstractMeasureLike, p::AbstractPolynomialLike)\ndot(p::AbstractPolynomialLike, μ::AbstractMeasureLike)\n\nSee expectation\n\n\n\n\n\n","category":"function"},{"location":"#MultivariateMoments-–-Multivariate-Moments-for-Julia","page":"Index","title":"MultivariateMoments –- Multivariate Moments for Julia","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Extension of MultivariatePolynomials to moments of multivariate measures and their scalar product with polynomials. It also includes the extraction of atomic measures from moment matrices which is related to the moment problem.","category":"page"},{"location":"#Contents","page":"Index","title":"Contents","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Pages = [\"moments.md\", \"atoms.md\"]\nDepth = 2","category":"page"},{"location":"atoms/#Atoms-extration","page":"Atoms extraction","title":"Atoms extration","text":"","category":"section"},{"location":"atoms/#Vectorized-matrix","page":"Atoms extraction","title":"Vectorized matrix","text":"","category":"section"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"SymMatrix\nVectorizedHermitianMatrix\nsquare_getindex\nsymmetric_setindex!","category":"page"},{"location":"atoms/#MultivariateMoments.SymMatrix","page":"Atoms extraction","title":"MultivariateMoments.SymMatrix","text":"struct SymMatrix{T} <: AbstractMatrix{T}\n    Q::Vector{T}\n    n::Int\nend\n\nSymmetric n times n matrix storing the vectorized upper triangular part of the matrix in the Q vector (this corresponds to the vectorized form of MathOptInterface.PositiveSemidefiniteConeTriangle). It implement the AbstractMatrix interface except for setindex! as it might break its symmetry. The symmetric_setindex! function should be used instead.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.VectorizedHermitianMatrix","page":"Atoms extraction","title":"MultivariateMoments.VectorizedHermitianMatrix","text":"struct VectorizedHermitianMatrix{T} <: AbstractMatrix{Complex{T}}\n    Q::Vector{T}\n    n::Int\nend\n\nHermitian n times n matrix storing the vectorized upper triangular real part of the matrix followed by the vectorized upper triangular imaginary part in the Q vector (this corresponds to the vectorized form of ComplexOptInterface.HermitianPositiveSemidefiniteConeTriangle). It implement the AbstractMatrix interface except for setindex! as it might break its symmetry. The symmetric_setindex! function should be used instead.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.square_getindex","page":"Atoms extraction","title":"MultivariateMoments.square_getindex","text":"square_getindex!(Q::SymMatrix, I)\n\nReturn the SymMatrix corresponding to Q[I, I].\n\n\n\n\n\nsquare_getindex!(Q::VectorizedHermitianMatrix, I)\n\nReturn the VectorizedHermitianMatrix corresponding to Q[I, I].\n\n\n\n\n\n","category":"function"},{"location":"atoms/#MultivariateMoments.symmetric_setindex!","page":"Atoms extraction","title":"MultivariateMoments.symmetric_setindex!","text":"symmetric_setindex!(Q::SymMatrix, value, i::Integer, j::Integer)\n\nSet Q[i, j] and Q[j, i] to the value value.\n\n\n\n\n\nsymmetric_setindex!(Q::VectorizedHermitianMatrix, value, i::Integer, j::Integer)\n\nSet Q[i, j] to the value value and Q[j, i] to the value -value.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#Moment-matrix","page":"Atoms extraction","title":"Moment matrix","text":"","category":"section"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"MomentMatrix\nmoment_matrix","category":"page"},{"location":"atoms/#MultivariateMoments.MomentMatrix","page":"Atoms extraction","title":"MultivariateMoments.MomentMatrix","text":"mutable struct MomentMatrix{T, B<:MultivariateBases.AbstractPolynomialBasis, MT<:AbstractMatrix{T}} <: AbstractMeasureLike{T}\n    Q::MT\n    basis::B\n    support::Union{Nothing, AlgebraicSet}\nend\n\nMeasure nu represented by the moments of the monomial matrix x x^top in the symmetric matrix Q. The set of points that are zeros of all the polynomials p such that mathbbE_nup = 0 is stored in the field support when it is computed.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.moment_matrix","page":"Atoms extraction","title":"MultivariateMoments.moment_matrix","text":"moment_matrix(μ::Measure, x)\n\nCreates a matrix the moment matrix for the moment matrix  x x^top using the moments of μ.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#Atomic-measure","page":"Atoms extraction","title":"Atomic measure","text":"","category":"section"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"WeightedDiracMeasure\nAtomicMeasure","category":"page"},{"location":"atoms/#MultivariateMoments.WeightedDiracMeasure","page":"Atoms extraction","title":"MultivariateMoments.WeightedDiracMeasure","text":"struct WeightedDiracMeasure{T}\n    center::Vector{T}\n    weight::T\nend\n\nRepresents a weighted dirac measure with centered at center and with weight weight.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.AtomicMeasure","page":"Atoms extraction","title":"MultivariateMoments.AtomicMeasure","text":"struct AtomicMeasure{T, AT, V} <: AbstractMeasureLike{T}\n    variables::V                           # Vector/Tuple of variables\n    atoms::Vector{WeightedDiracMeasure{T, AT}} # Atoms of the measure\nend\n\nAn measure is said to be atomic if it is a sum of weighted dirac measures. For instance, eta = 2 delta_(1 0) + 3 delta_(12 12) is an atomic measure since it is a sum of the diracs centered at (1 0) and (12 12) and weighted respectively by 2 and 3. That is, mathbbE_etap = 2p(1 0) + 3p(12 12).\n\nThe AtomicMeasure struct stores an atomic measure where variables contains the variables and atoms contains atoms of the measure.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#Atoms-extraction","page":"Atoms extraction","title":"Atoms extraction","text":"","category":"section"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"Given a MomentMatrix containing the moments of an atomic measure, atomic_measure attempts to recover the dirac centers and weights by first computing an algebraic system with the atom centers as solution with compute_support!.","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"compute_support!\natomic_measure","category":"page"},{"location":"atoms/#MultivariateMoments.compute_support!","page":"Atoms extraction","title":"MultivariateMoments.compute_support!","text":"compute_support!(ν::MomentMatrix, rank_check, [dec])\n\nComputes the support field of ν. The rank_check and dec parameters are passed as is to the low_rank_ldlt function.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#MultivariateMoments.atomic_measure","page":"Atoms extraction","title":"MultivariateMoments.atomic_measure","text":"atomic_measure(ν::MomentMatrix, rank_check::RankCheck, [dec::LowRankLDLTAlgorithm], [solver::SemialgebraicSets.AbstractAlgebraicSolver])\n\nReturn an AtomicMeasure with the atoms of ν if it is atomic or nothing if ν is not atomic. The rank_check and dec parameters are passed as is to the low_rank_ldlt function. By default, dec is an instance of SVDLDLT. The extraction relies on the solution of a system of algebraic equations. using solver. For instance, given a MomentMatrix, μ, the following extract atoms using a rank_check of 1e-4 for the low-rank decomposition and homotopy continuation to solve the obtained system of algebraic equations.\n\nusing HomotopyContinuation\nsolver = SemialgebraicSetsHCSolver(; compile = false)\natoms = atomic_measure(ν, 1e-4, solver)\n\nIf no solver is given, the default solver of SemialgebraicSets is used which currently computes the Gröbner basis, then the multiplication matrices and then the Schur decomposition of a random combination of these matrices. For floating point arithmetics, homotopy continuation is recommended as it is more numerically stable than Gröbner basis computation.\n\n\n\n\n\n","category":"function"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"For the first step of compute_support!, there are two approaches. The first one is to exploit the flat extension to directly obtain the multiplication matrices.","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"FlatExtension","category":"page"},{"location":"atoms/#MultivariateMoments.FlatExtension","page":"Atoms extraction","title":"MultivariateMoments.FlatExtension","text":"struct FlatExtension{\n    MMS<:SemialgebraicSets.AbstractMultiplicationMatricesSolver,\n}\n    multiplication_matrices_solver::MMS\nend\n\nGiven a moment matrix satisfying the flat extension property described in [L09, Section 5.3], computes the multiplication matrices using the formula given in [L09, Lemma 6.21] or [LLR08, Section 4.4.4]. Given the multiplication matrices, the solutions are obtained with multiplication_matrices_solver.\n\n[L09] Laurent, Monique. Sums of squares, moment matrices and optimization over polynomials. Emerging applications of algebraic geometry (2009): 157-270.\n\n[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp. Semidefinite characterization and computation of zero-dimensional real radical ideals. Foundations of Computational Mathematics 8 (2008): 607-647.\n\n\n\n\n\n","category":"type"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"The second approach is to first obtain the image space of the moment matrix, represented as a MacaulayNullspace and to then compute the multiplication matrices from this image space. This image space is obtained from a low rank LDLT decomposition of the moment matrix. This decomposition can either be obtained by a cholesky or SVD decomposition from which we remove the rows corresponding to the negligeable eigen/singular values.","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"LowRankLDLTAlgorithm\nShiftCholeskyLDLT\nSVDLDLT\nlow_rank_ldlt\nLowRankLDLT\nMacaulayNullspace","category":"page"},{"location":"atoms/#MultivariateMoments.LowRankLDLTAlgorithm","page":"Atoms extraction","title":"MultivariateMoments.LowRankLDLTAlgorithm","text":"LowRankLDLTAlgorithm\n\nMethod for computing a r times n matrix U of a n times n rank r psd matrix Q such that Q = U'U.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.ShiftCholeskyLDLT","page":"Atoms extraction","title":"MultivariateMoments.ShiftCholeskyLDLT","text":"ShiftCholeskyLDLT <: LowRankLDLTAlgorithm\n\nShift the matrix by shift times the identity matrix before cholesky.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.SVDLDLT","page":"Atoms extraction","title":"MultivariateMoments.SVDLDLT","text":"SVDLDLT <: LowRankLDLTAlgorithm\n\nUse SVD decomposition.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.low_rank_ldlt","page":"Atoms extraction","title":"MultivariateMoments.low_rank_ldlt","text":"MultivariateMoments.low_rank_ldlt(Q::AbstractMatrix, dec::LowRankLDLTAlgorithm, ranktol)\n\nReturns a r times n matrix U of a n times n rank r positive semidefinite matrix Q such that Q = U^top U. The rank of Q is the number of singular values larger than ranktol cdot sigma_1 where sigma_1 is the largest singular value.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#MultivariateMoments.LowRankLDLT","page":"Atoms extraction","title":"MultivariateMoments.LowRankLDLT","text":"struct LowRankLDLT{T,Tr,C<:RankCheck}\n    L::Matrix{T}\n    singular_values::Vector{Tr}\n    rank_check::C\nend\n\nLow-Rank cholesky decomposition L * Diagonal(singular_values) * L' of size (n, r) of a n x n matrix with singular values singular_values[1] > ... > singular_values[n]. The rank was chosen given singular_values using rank_check via the rank_from_singular_values function.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.MacaulayNullspace","page":"Atoms extraction","title":"MultivariateMoments.MacaulayNullspace","text":"struct MacaulayNullspace{T,MT<:AbstractMatrix{T},BT}\n    matrix::MT\n    basis::BT\nend\n\nThis matrix with rows indexed by basis can either be interpreted as the right null space of a Macaulay matrix with columns indexed by basis (resp. or the image space of a moment matrix with rows and columns indexed by basis). The value of matrix[i, j] should be interpreted as the value of the ith element of basis for the jth generator of the null space (resp. image) space.\n\n\n\n\n\n","category":"type"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"The choice of cutoff between the significant and neglibeable eigen/singular values is parametrized by the following interface:","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"RankCheck\nrank_from_singular_values\naccuracy\ndoubt","category":"page"},{"location":"atoms/#MultivariateMoments.RankCheck","page":"Atoms extraction","title":"MultivariateMoments.RankCheck","text":"abstract type RankCheck end\n\nMethod for computing the rank with rank_from_singular_values based on a list of singular values.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.rank_from_singular_values","page":"Atoms extraction","title":"MultivariateMoments.rank_from_singular_values","text":"rank_from_singular_values(σ, check::RankCheck)\n\nReturn the rank of a matrix with singular values σ (in decreasing order) using check.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#MultivariateMoments.accuracy","page":"Atoms extraction","title":"MultivariateMoments.accuracy","text":"accuracy(σ, r, check::RankCheck)\n\nReturns a value measuring the accuracy of the rank check check returning rank r. This is used by Echelon to determine the accuracy to use for the Gaussian elimination.\n\n\n\n\n\naccuracy(chol::LowRankLDLT)\n\nReturn the ratio rtol = σ[r+1]/σ[1] where σ is the vector of singular values of the matrix for which chol is the rank-r Cholesky decomposition. This is a good relative tolerance to use with the matrix as σ[r+1] is the first singular value that was discarded.\n\n\n\n\n\n","category":"function"},{"location":"atoms/#MultivariateMoments.doubt","page":"Atoms extraction","title":"MultivariateMoments.doubt","text":"doubt(σ, check::RankCheck)\n\nReturns a value measuring the doubt of the rank check check. Lower values means more doubt so less certainty. This is used by FallbackRank for deciding whether the fallback should be used.\n\n\n\n\n\n","category":"function"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"The rank check can be chosen among the following:","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"UserRank\nFixedRank\nFixedRanks\nAbsoluteRankTol\nLeadingRelativeRankTol\nDifferentialRankTol\nLargestDifferentialRank\nFallbackRank","category":"page"},{"location":"atoms/#MultivariateMoments.UserRank","page":"Atoms extraction","title":"MultivariateMoments.UserRank","text":"struct UserRank <: RankCheck\n    pagesize::Int\nend\n\nThe user chooses the rank given the singular values in a REPL.TerminalMenus.RadioMenu of page size pagesize.\n\nExample\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], UserRank())\nChoose the last significant singular value:\n   1.0\n   0.1\n > 0.05\n   1.0e-5\n   5.0e-6\n3\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.FixedRank","page":"Atoms extraction","title":"MultivariateMoments.FixedRank","text":"struct FixedRank <: RankCheck\n    r::Int\nend\n\nThe rank is hardcoded to r, independently of the singular values.\n\nExample\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], FixedRank(3))\n3\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.FixedRanks","page":"Atoms extraction","title":"MultivariateMoments.FixedRanks","text":"mutable struct FixedRanks <: RankCheck\n    r::Vector{Int}\n    current::Int\nend\n\nThe ith rank is hardcoded to r[i], independently of the singular values. The field current indicates how many ranks have already been asked. When current is length(r), no rank can be asked anymore.\n\nExample\n\njulia> check = FixedRanks([2, 3])\nFixedRanks([2, 3], 0)\n\njulia> rank_from_singular_values([1, 1e-1, 5e-5, 1e-5, 5e-6], check)\n2\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], check)\n3\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.AbsoluteRankTol","page":"Atoms extraction","title":"MultivariateMoments.AbsoluteRankTol","text":"struct AbsoluteRankTol{T} <: RankCheck\n    tol::T\nend\n\nThe rank is the number of singular values that are strictly larger than tol.\n\nExample\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], AbsoluteRankTol(1e-4))\n3\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.LeadingRelativeRankTol","page":"Atoms extraction","title":"MultivariateMoments.LeadingRelativeRankTol","text":"struct LeadingRelativeRankTol{T} <: RankCheck\n    tol::T\nend\n\nThe rank is the number of singular values that are strictly larger than tol * maximum(σ) where maximum(σ) is the largest singular value.\n\nExample\n\nWhen the matrix is obtained from a homogeneous problem where the scaling is irrelevant, LeadingRelativeRankTol may be preferable to AbsoluteRankTol as shown below\n\njulia> rank_from_singular_values(1e6 * [1, 1e-1, 5e-2, 1e-5, 5e-6], AbsoluteRankTol(1e-4))\n5\n\njulia> rank_from_singular_values(1e6 * [1, 1e-1, 5e-2, 1e-5, 5e-6], LeadingRelativeRankTol(1e-4))\n3\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.DifferentialRankTol","page":"Atoms extraction","title":"MultivariateMoments.DifferentialRankTol","text":"struct DifferentialRankTol{T} <: RankCheck\n    tol::T\nend\n\nThe rank is the number of singular values before a singular value (not included) is tol times the previous one (included).\n\nExample\n\nIt is sometimes difficult to figure out the tolerance to use in LeadingRelativeRankTol. For instance, choosing 1e-3 will consider 1e-3 in the example below as not part of the rank while DifferentialRankTol would include it because it is close to the previous singular value.\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-6, 5e-7], LeadingRelativeRankTol(1e-3))\n5\n\njulia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-6, 5e-7], DifferentialRankTol(1e-2))\n6\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.LargestDifferentialRank","page":"Atoms extraction","title":"MultivariateMoments.LargestDifferentialRank","text":"struct LargestDifferentialRank <: RankCheck\nend\n\nThe rank is the number of singular values until the singular value that has the largest ratio with the next singular value.\n\nExample\n\nIt is sometimes difficult to figure out the tolerance to use in DifferentialRankTol. For instance, choosing 1e-2 will consider 1e-2, 5e-2 and 1e-3 in the example below as not part of the rank while LargestDifferentialRank would include them because there is a largest gap later.\n\njulia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], DifferentialRankTol(1e-2))\n1\n\njulia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], LargestDifferentialRank())\n4\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.FallbackRank","page":"Atoms extraction","title":"MultivariateMoments.FallbackRank","text":"struct FallbackRank{T,D,F} <: RankCheck\n    tol::T\n    default::D\n    fallback::F\nend\n\nDefaults to checking the rank with default and falls back to fallback if the doubt is strictly larger than tol. By default, fallback is UserRank.\n\nExample\n\nThe advantage of UserRank is that the user get to see if the rank check is ambiguous and act accordingly. The downside is that it might be cumbersome if there are many rank checks to do. With FallbackRank, the user only has to sort out the tricky ones. In the example below, the first example is handled by LargestDifferentialRank. For the second one, the user sees that this is a tricky choice and can manually choose one of the two ranks, then see the result of the rest of his code using this value of the code and then choose the other rank and see the impact of this different choice.\n\njulia> check = FallbackRank(1e-1, LargestDifferentialRank())\nFallbackRank{Float64, LargestDifferentialRank, UserRank}(0.1, LargestDifferentialRank(), UserRank(8))\n\njulia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], check)\n4\n\njulia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 5e-6, 5e-7], check)\nChoose the last significant singular value:\n > 1.0\n   0.01\n   0.05\n   0.001\n   5.0e-6\n   5.0e-7\n1\n\njulia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 5e-6, 5e-7], check)\nChoose the last significant singular value:\n   1.0\n   0.01\n   0.05\n > 0.001\n   5.0e-6\n   5.0e-7\n4\n\n\n\n\n\n","category":"type"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"Given the MacaulayNullspace, there are two approaches implemented to obtain the moment matrices:","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"ShiftNullspace\nEchelon","category":"page"},{"location":"atoms/#MultivariateMoments.ShiftNullspace","page":"Atoms extraction","title":"MultivariateMoments.ShiftNullspace","text":"struct ShiftNullspace end\n\nFrom the MacaulayNullspace, computes multiplication matrices by exploiting the shift property of the rows [DBD12].\n\n[DBD12] Dreesen, Philippe, Batselier, Kim, and De Moor, Bart. Back to the roots: Polynomial system solving, linear algebra, systems theory. IFAC Proceedings Volumes 45.16 (2012): 1203-1208.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.Echelon","page":"Atoms extraction","title":"MultivariateMoments.Echelon","text":"struct Echelon\nend\n\nGiven a MacaulayNullspace, computes its echelon form (corresponding to the Canonical Null Space of [D13]) with Gaussian elimination. From this echelon form, the left null space can easily be computed using using [HL05, (8)]. This left null space forms a system of polynomial equations.\n\nnote: Note\nIn the context of compute_support!, if the MacaulayNullspace was obtained using SVDLDLT, the left null space can easily be obtained from the singular vectors corresponding to the negligeable singular values that were removed. However, as mentioned [LLR08, Section 4.4.5], these will usually give an overdetermined bases. As shown in [S04, Section 10.8.2], it is desirable to avoid overdetermined bases because it could lead to inconsistencies in the basis for numerical reasons. For this reason, this method computes the left null space from the echelon form of the significant singular values instead.\n\nLet B be the set of monomials corresponding to the rows of the pivots of this echelon form.  If the moment matrix satisfies the flat extension property described in [L09, Section 5.3] then all monomials of the border of B (as defined in [LLR08, (2.3)]) will correspond to a row of of the matrix. In that case, the polynomial of the system obtained by [HL05, (8)] form a rewriting family for B [L09, (2.16)] a.k.a. a B-border prebasis (as defined in [LLR08, (2.4)]). Therefore, they can be used to compute multiplication matrices.\n\n[HL05] Henrion, D. & Lasserre, J-B. Detecting Global Optimality and Extracting Solutions of GloptiPoly 2005\n\n[D13] Dreesen, Philippe. Back to the Roots: Polynomial System Solving Using Linear Algebra Ph.D. thesis (2013)\n\n[L09] Laurent, Monique. Sums of squares, moment matrices and optimization over polynomials. Emerging applications of algebraic geometry (2009): 157-270.\n\n[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp. Semidefinite characterization and computation of zero-dimensional real radical ideals. Foundations of Computational Mathematics 8 (2008): 607-647.\n\n[S04] Stetter, Hans J. Numerical polynomial algebra. Society for Industrial and Applied Mathematics, 2004.\n\n\n\n\n\n","category":"type"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"The Echelon uses the RowEchelon package to determine the standard monomials (which is not numerically stable) while the ShiftNullspace uses the following function internally which is based on SVD so it should have better numerical behavior.","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"standard_monomials_and_border","category":"page"},{"location":"atoms/#MultivariateMoments.standard_monomials_and_border","page":"Atoms extraction","title":"MultivariateMoments.standard_monomials_and_border","text":"function standard_monomials_and_border(\n    null::MacaulayNullspace,\n    rank_check,\n)\n\nComputes the set of standard monomials using the greedy sieve algorithm presented in [LLR08, Algorithm 1].\n\n[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp. Semidefinite characterization and computation of zero-dimensional real radical ideals. Foundations of Computational Mathematics 8 (2008): 607-647.\n\n\n\n\n\n","category":"function"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"Once the center of the atoms are determined, a linear system is solved to determine the weights corresponding to each dirac. By default, MomentMatrixWeightSolver is used by atomic_measure so that if there are small differences between moment values corresponding to the same monomial in the matrix (which can happen if these moments were computed numerically by a semidefinite proramming solvers, e.g., with SumOfSquares), the linear system handles that automatically.","category":"page"},{"location":"atoms/","page":"Atoms extraction","title":"Atoms extraction","text":"MomentMatrixWeightSolver\nMomentVectorWeightSolver","category":"page"},{"location":"atoms/#MultivariateMoments.MomentMatrixWeightSolver","page":"Atoms extraction","title":"MultivariateMoments.MomentMatrixWeightSolver","text":"struct MomentMatrixWeightSolver\n    rtol::T\n    atol::T\nend\n\nGiven a moment matrix ν and the atom centers, determine the weights by solving a linear system over all the moments of the moment matrix, keeping duplicates (e.g., entries corresponding to the same monomial).\n\nIf the moment values corresponding to the same monomials are known to be equal prefer MomentVectorWeightSolver instead.\n\n\n\n\n\n","category":"type"},{"location":"atoms/#MultivariateMoments.MomentVectorWeightSolver","page":"Atoms extraction","title":"MultivariateMoments.MomentVectorWeightSolver","text":"struct MomentVectorWeightSolver{T}\n    rtol::T\n    atol::T\nend\n\nGiven a moment matrix ν and the atom centers, first convert the moment matrix to a vector of moments, using measure(ν; rtol=rtol, atol=atol) and then determine the weights by solving a linear system over the monomials obtained.\n\nIf the moment values corresponding to the same monomials can have small differences, measure can throw an error if rtol and atol are not small enough. Alternatively to tuning these tolerances MomentVectorWeightSolver can be used instead.\n\n\n\n\n\n","category":"type"}]
}
