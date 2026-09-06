module Finance

import Statistics


export sig_cumsum, tic_diff1, tic_diff2, isConvertible
export ema, ema_std, ema_stats, sample_std
export entropy_index, pow_n, ewt_mean

# The floating point type used for results computed from values of type `T`.
_ftype(::Type{T}) where {T<:Real} = float(T)
_ftype(::Type{S}, ::Type{T}) where {S<:Real, T<:Real} = float(promote_type(S, T))


"""
    isConvertible(S, T)

Boolean function which returns `true` if a value of type `S` 
can be converted to a value of type `T`.

# Type Constraints
- `S <: Real`
- `T <: Real`

# Arguments
- `::Type{S}` -- A numeric type.
- `::Type{T}` -- A numeric type.

# Return
`::Bool`
"""
function isConvertible(::Type{S}, ::Type{T}) where {S<:Real,T<:Real}
    try
        convert(T, one(S))
    catch _
        return (false)
    end

    return (true)
end


"""
    tic_diff1(t, x; chk_inp=false)

Compute the numerical derivative of a function represented by `x`
with respect to `t` when the values in `t` are possibly irregular.
The three point (second order accurate) formula for a non-uniform grid is used;
the derivative is returned at the interior points `t[2:end-1]`.

# Type Constraints
- `S <: Real`
- `T <: Real`

# Arguments
- `t :: AbstractVector{T}`   -- A vector of times.
- `x :: AbstractVector{S}`   -- A vector of values.


# Keyword Arguments
- `chk_inp=false :: Bool`  -- Check the input contract?

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint                     | Description                                                               |
|:------------------------------:|:------------------------------------------------------------------------- |
| `\\|t\\| = \\|x\\|`            | The length of the time and data vectors match.                            |
| `\\|x\\|` ``\\ge`` `2`          | At least two points (always checked).                                     |
| ``\\forall i, t_{i+1} > t_i``  | The times are increasing; consequently, we have a 1-1 map from `t` to `x`.|

# Return
`:: Vector{F}` with `F = float(promote_type(S, T))`, of length `|x| - 2`.
"""
@noinline function tic_diff1(t::AbstractVector{T},
    						 x::AbstractVector{S};
    					     chk_inp::Bool=false  ) where {S<:Real,T<:Real}
    n = length(x)
    n >= 2 || throw(DomainError(n, "The data series must have at least 2 points."))

    if chk_inp
        n == length(t)          || throw(DomainError(length(t), "The length of the time and data series must match."))
        all(diff(t) .> zero(T)) || throw(DomainError(t, "The time series must be strictly increasing."))
    end

    F  = _ftype(S, T)
    df = Vector{F}(undef, n - 2)
    @inbounds for i in 2:(n-1)
        h1 = F(t[i]) - F(t[i-1])
        h2 = F(t[i+1]) - F(t[i])
        df[i-1] = (h1 * h1 * x[i+1] - h2 * h2 * x[i-1] + (h2 * h2 - h1 * h1) * x[i]) / (h1 * h2 * (h1 + h2))
    end

    return (df)
end

"""
    tic_diff2(t, x; chk_inp=false)

Compute the numerical second derivative of a function represented by `x`
with respect to `t` when the values in `t` are possibly irregular.

# Type Constraints
- `T <: Real`
- `S <: Real`

# Arguments
- `t :: AbstractVector{T}` -- A vector of times.
- `x :: AbstractVector{S}` -- A vector of values.

# Keyword Arguments
- `chk_inp=false :: Bool`  -- Check the input contract?

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint                     | Description                                                               |
|:------------------------------:|:--------------------------------------------------------------------------|
| `\\|t\\| = \\|x\\|`            | The length of the time and data vectors match.                            |
| `\\|x\\|` ``\\ge`` `2`          | At least two points (always checked).                                     |
| ``\\forall i, t_{i+1} > t_i``  | The times are increasing; consequently, we have a 1-1 map from `t` to `x`.|

# Return
`:: Vector{F}` with `F = float(promote_type(S, T))`, of length `|x| - 2`.
"""
@noinline function tic_diff2(t::AbstractVector{T},
    						 x::AbstractVector{S};
    						 chk_inp::Bool=false  ) where {T<:Real,S<:Real}
    n = length(x)
    n >= 2 || throw(DomainError(n, "The data series must have at least 2 points."))

    if chk_inp
        n == length(t)          || throw(DomainError(length(t), "The length of the time and data series must match."))
        all(diff(t) .> zero(T)) || throw(DomainError(t, "The time series must be strictly increasing."))
    end

    F  = _ftype(S, T)
    df = Vector{F}(undef, n - 2)
    @inbounds for i in 2:(n-1)
        h1 = F(t[i]) - F(t[i-1])
        h2 = F(t[i+1]) - F(t[i])
        df[i-1] = 2 * (h1 * x[i+1] - (h1 + h2) * x[i] + h2 * x[i-1]) / (h1 * h2 * (h1 + h2))
    end

    return (df)
end



"""
    sig_cumsum(t, x, w, h; chk_inp=false)

Return a tuple of two vectors: tics, signals.

The signals are the collections of all deviations from 
a running mean (with window length `w`) of the series. 
The threshold of the deviation is `h`.

Deviation is determined by:
- ``S_t^+ = {\\rm max}(0, S_{t-1} + x_t - E[x_{t-1}]; S^+_0 = 0``
- ``S_t^- = {\\rm min}(0, S_{t-1} + x_t - E[x_{t-1}]; S^-_0 = 0``
- ``S_t^{\\hphantom{+}} = {\\rm max}(S^+_t, -S^-_t)``

Collect all ``t, S_t`` where ``S_t \\ge h``.

# Type Constraints
- `S <: Real`
- `T <: Real`

# Arguments
- `t :: AbstractVector{S}` -- The tic series to examine.
- `x :: AbstractVector{T}` -- The series to examine.
- `w :: Int`             -- The width of the moving average.
- `h :: T`                 -- The threshold for the deviation to register.

# Keyword Arguments
- `chk_inp=false :: Bool`  -- Check the input contract?

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint                          | Description                                                               |
|:-----------------------------------:|:--------------------------------------------------------------------------|
|`\\|x\\|` ``\\ge`` `2`               | The length of x is at least ``2``.                                        | 
|`\\|x\\| = \\|t\\|`                  | The length of `x` is equal to the length of `t`.                          |
|`w > 1`                              | The window length is greater than ``1``.                                  |
|`h > 0`                              | The deviation threshold is greater than 0.                                |
|``\\forall i, t_{i+1} > t_{i}``      | The times are increasing; consequently, we have a 1-1 map from `t` to `x`.|

# Output Components
- `td :: Vector{S}` -- Values of `t` where deviations occurred.
- `xd :: Vector{F}` -- Values of the deviations (`F = float(T)`).

# Output Contract
- `|td| = |xd|`

# Return
`(td, xd) :: Tuple{Vector{S}, Vector{F}}`
"""
function sig_cumsum(t::AbstractVector{S},
    				x::AbstractVector{T},
    				w::Int            ,
    				h::Real             ;
    				chk_inp::Bool=false  ) where {S<:Real,T<:Real}
    n = length(x)

    # Input contract.
    if chk_inp
        n >= 2                   || throw(DomainError(n, "Vector length of `x` must be >= 2."))
		n == length(t)           || throw(DomainError(length(t), "Length of time sequence should match length of data."))
        w > 1                    || throw(DomainError(w, "Window length must be > 1."))
        h > 0                    || throw(DomainError(h, "Deviation threshold must be > 0."))
		all(diff(t) .> zero(S))  || throw(DomainError(t, "Sequential differences of time seq must always be > 0."))
    end

    F  = _ftype(T)
    z  = zero(F)
    Sp = z          # Running positive deviation.
    Sn = z          # Running negative deviation.

    xm   = F(x[1])  # Running mean of the input `x` computed based on window, `w`.
    sigs = F[]      # The signals/deviations to be returned. 
    tics = S[]      # The tics where the deviations occurred.

    # Loop over the series and populate, `tics` and `sigs`.
    @inbounds for i in 2:n
        xm = ((w - 1) * xm + x[i-1]) / w
        Sp = max(z, Sp + x[i] - xm)
        Sn = min(z, Sn + x[i] - xm)
        delta = max(Sp, -Sn)
        if delta >= h
            push!(sigs, delta)
            push!(tics, t[i])
        end
    end

    return ((tics, sigs))
end



# The default half-life for a decay window of length `m`.
default_half_life(m::Int) = max(2, div(m, 2))

"""
    ema(x,m; h=max(2, div(m,2)))

Compute the Exponential Moving Average of the sequence `x`.

The weights decay by the factor ``\\lambda = 2^{-1/h}`` per step
(so a weight halves every `h` steps) and are normalized over the window of length `m`.

# Type Constraints
- `T <: Real`

# Arguments
- `x :: AbstractVector{T}` -- The series to work with.
- `m :: Int`             -- The width of the decay window.

# Keyword Arguments
- `h=max(2, div(m,2)) :: Int` -- The exponential decay *half-life*. 

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint | Description                                   | 
|:----------:|:----------------------------------------------|
| `m > 1`    | Averaging window length is greater than ``1``.|
| `h > 1`    | Exponential *half-life* is greater than ``1``.|
| `\\|x\\| > 0` | The series is not empty.                 |

# Output 
- `ema :: Vector{F}` -- The exponential moving average of `x` (`F = float(T)`).

# Output Contract
- `|x| = |ema|`

# Return
`ema::Vector{F}`

"""
@noinline function ema(x::AbstractVector{T},
    				   m::Int              ;
    				   h::Int=default_half_life(m)) where {T<:Real}

    # Check Input Contract
    m > 1 || throw(DomainError(m, "The window length must be > 1."))
    h > 1 || throw(DomainError(h, "The half-life must be > 1."))
    isempty(x) && throw(DomainError(0, "The data series must not be empty."))

    F = _ftype(T)
    N = length(x)
	ma = Vector{F}(undef, N)
	xadj = Vector{F}(undef, N + m)
    @inbounds xadj[1:m] .= x[1]
    @inbounds xadj[(m+1):end] = x
	w = Vector{F}(undef, m)

    # Term by term decay factor.
    l = exp(-log(2 * one(F)) / h)

    w[1] = l
    @inbounds @simd for i in 2:m
        w[i] = l * w[i-1]
    end
    w ./= sum(w)

    # Compute the EMA using the difference equation recursion.
    ma[1] = xadj[m+1]
    @inbounds @simd for i in 2:N
        ma[i] = l * (ma[i-1] - w[m] * xadj[i]) + w[1] * xadj[i+m]
    end

    return (ma)
end


"""
    ema_std(x, m; h=max(2, div(m,2)), init_sig=nothing)

Compute the Moving Exponential Standard Deviation of the sequence `x`.
By default, the initial std is taken to be the standard deviation of 
the first window (of length `m`). However, a user specified value
may be used instead.

# Type Constraints
- `T <: Real`

# Arguments
- `x :: AbstractVector{T}` -- The series to work with.
- `m :: Int`             -- The width of the decay window.

# Keyword Arguments
- `h=max(2, div(m,2)) :: Int`             -- The exponential decay *half-life* (see `ema`). 
- `init_sig=nothing :: Union{Nothing, Real}` -- An optional user supplied initial standard deviation for the start of the series.      

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint           | Description                                             |   
|:--------------------:|:------------------------------------------------------- |
| `m > 1`              | Averaging window length is greater than ``1``.          |
| `h > 1`              | Exponential *half-life* is greater than ``1``.          |
| `\\|x\\| > 1`        | The length of the series is greater than ``1``.         |
| `init_sig` ``\\ge 0``| User supplied starting ``\\sigma`` should be ``\\ge 0``.|

# Output 
- `stda :: Vector{F}` -- The moving exponential standard deviation of `x` (`F = float(T)`).

# Output Contract
- `|x| = |stda|`

# Return
`stda::Vector{F}`
"""
@noinline function ema_std(x::AbstractVector{T}             ,
    					   m::Int                           ;
    					   h::Int=default_half_life(m)      ,
    					   init_sig::Union{Nothing,Real}=nothing) where {T<:Real}

    N = length(x)

    # Check input constraints.
    m > 1 || throw(DomainError(m, "The window length must be > 1."))
    h > 1 || throw(DomainError(h, "The half-life must be > 1."))
    N > 1 || throw(DomainError(N, "The length of the data series must be > 1."))
	(init_sig === nothing || init_sig >= 0) || throw(DomainError(init_sig, "The initial sigma must be `nothing` or non-negative."))

    F = _ftype(T)

    # Compute the ema for `x`.
    ma = ema(x, m, h=h)

    # Variance estimates.
	mvar = Vector{F}(undef, N)

    # Set the initial estimated/supplied variance.
    mvar[1] = init_sig !== nothing ? F(init_sig) : sample_std(view(x, 1:min(m, N)))
    mvar[1] *= mvar[1]

    # Add history (`m` zeros) for variance.
    # We do this by augmenting the length of the "x"'s -- `(x - ma)^2`
    # to have size `N + m`, so we can go "back" m. This means that
    # xadj has to be indexed differently than the way 
    # the formula does indexing.
	xadj = Vector{F}(undef, N + m)
	xadj[1:m] .= zero(F)

    # We don't need the following line (like what we have in the corresponding ema code)
    # as `(x - ma)[1]` = 0, and the xadj array is already set to 0.
    d = x .- ma
    @inbounds xadj[(m+1):end] .= d .* d
	w = Vector{F}(undef, m)

    # Term by term decay factor.
    l = exp(-log(2 * one(F)) / h)

    # Use this to define the weights; then normalize.
    # Weights go from large to small.
    w[1] = l
    @inbounds @simd for i in 2:m
        @fastmath w[i] = l * w[i-1]
    end
    w ./= sum(w)
    w2 = sum(w .* w)

    # Recursive formula for variance.
    @inbounds @simd for n in 1:(N-1)
        mvar[n+1] = l * (mvar[n] - xadj[n+1] * w[m]) + xadj[n+m+1] * w[1]
    end

    # Return corrected variances (unbiased).
    return (sqrt.(mvar ./ (one(F) - w2)))
end


"""
    ema_stats(x, m; h=max(2, div(m,2)), init_sig=nothing)

Compute the Moving Exponential Stats of the 
sequence `x`: `ema`, `ema_std`, `ema_rel_skew`, `ema_rel_kurtosis`.

The recursive formulas for the moving statistics come from the paper:
[exponential\\_moving\\_average.pdf](https://github.com/scottrsm/math/tree/main/pdf/exponential_moving_average.pdf).

The estimates are corrected for the bias of weighted moments of independent samples.
With normalized weights ``w_i`` and ``W_k = \\sum_i w_i^k``, the weighted central moments
``M_k = \\sum_i w_i (x_i - \\bar x)^k`` of independent samples with central moments
``\\mu_k`` and variance ``\\sigma^2`` satisfy
- ``E[M_2] = \\sigma^2 (1 - W_2)``
- ``E[M_3] = \\mu_3 (1 - 3W_2 + 2W_3)``
- ``E[M_4] = \\mu_4 (1 - 4W_2 + 6W_3 - 3W_4) + 3\\sigma^4 (4W_3 - 3W_2^2 - W_4)``

which are inverted to estimate ``\\sigma``, ``\\mu_3`` and ``\\mu_4``. The relative skew is
``\\mu_3 / \\sigma^3`` and the relative (excess) kurtosis is ``\\mu_4 / \\sigma^4 - 3``
(zero for a normal distribution). Both are dimensionless: rescaling `x` does not change them.

Returns these stats as a matrix with four columns, each representing the stats above in the order listed.

# Type Constraints
- `T <: Real`

# Arguments
- `x :: AbstractVector{T}` -- The series to work with.
- `m :: Int`             -- The width of the decay window.

# Keyword Arguments
- `h=max(2, div(m,2)) :: Int` -- The exponential decay *half-life* (see `ema`). 
- `init_sig=nothing:: Union{Real, Nothing}` -- An optional user supplied initial standard deviation for the start of the series.      

# Input Contract
The inputs are assumed to satisfy the constraints below:

| Constraint     | Description                                    | 
|:--------------:|:---------------------------------------------- |
| `m > 1`        | Averaging window length is greater than ``1``. |
| `h > 1`        | Exponential *half-life* is greater than ``1``. |
| `\\|x\\| > 3`  | The length of the series is greater than ``3``.|

# Output 
- `stat :: Matrix{F}` -- A matrix of EMA stats: `ema`, `ema_std`, `ema_rel_skew`, `ema_rel_kurtosis` (`F = float(T)`).

# Output Contract
- `|stat| = (N, 4)` 

# Return
`stat::Matrix{F}`
"""
@noinline function ema_stats(x::AbstractVector{T}           ,
    					   	 m::Int                         ;
    						 h::Int=default_half_life(m)    ,
    						 init_sig::Union{Nothing,Real}=nothing) where {T<:Real}
    N = length(x)

    # Check input constraints.
    m > 1 || throw(DomainError(m, "The window length must be > 1."))
    h > 1 || throw(DomainError(h, "The half-life must be > 1."))
	N > 3 || throw(DomainError(N, "N must be > 3."))
    (init_sig === nothing || init_sig >= 0) || throw(DomainError(init_sig, "The initial sigma must be `nothing` or non-negative."))

    F = _ftype(T)

    # Compute the EMA of `x`.
    ma = ema(x, m, h=h)

	mstat = Matrix{F}(undef, N, 4)
    mstat[1, 1] = x[1]

    # Set the initial estimated/supplied variance.
    mstat[1, 2] = init_sig !== nothing ? F(init_sig) : sample_std(view(x, 1:min(m, N)))
    mstat[1, 2] *= mstat[1, 2]
    mstat[1, 3] = zero(F)
    mstat[1, 4] = zero(F)

    # Add history (`m` zeros) for variance, etc.
    # We do this by augmenting the length of the "x"'s -- `(x - ma)^2`, `(x - ma)^3`, etc.
    # to have size `N + m`, so we can go "back" `m`. This means that
    # `xadj` has to be indexed differently than the way 
    # the formula does indexing.
	xadj = Matrix{F}(undef, N + m, 4)
    d = x .- ma
    v = d .* d
    @inbounds xadj[1:m, 1] .= x[1]
    @inbounds xadj[1:m, 2] .= zero(F)
    @inbounds xadj[1:m, 3] .= zero(F)
    @inbounds xadj[1:m, 4] .= zero(F)
    @inbounds xadj[(m+1):end, 1] .= x
    @inbounds xadj[(m+1):end, 2] .= v
    @inbounds xadj[(m+1):end, 3] .= v .* d
    @inbounds xadj[(m+1):end, 4] .= v .* v
	w = Vector{F}(undef, m)

    # Term by term decay factor.
    l = exp(-log(2 * one(F)) / h)

    w[1] = l
    @inbounds for i in 2:m
        w[i] = l * w[i-1]
    end
    w ./= sum(w)

    # Compute the sums of `w` to powers from 2 to 4.
    W2 = zero(F)
    W3 = zero(F)
    W4 = zero(F)
    @inbounds @simd for i in 1:m
        wt = w[i]
        w2 = wt * wt
        W2 += w2
        W3 += w2 * wt
        W4 += w2 * w2
    end

    # Expressions needed to unbias our estimates (see the docstring).
    B2 = one(F) - W2                              # E[M2] = σ² B2
    B3 = one(F) - 3 * W2 + 2 * W3                 # E[M3] = μ3 B3
    B4 = one(F) - 4 * W2 + 6 * W3 - 3 * W4        # E[M4] = μ4 B4 + 3 σ⁴ A4
    A4 = 4 * W3 - 3 * W2 * W2 - W4

    # Recursion to compute the moving stats.
    for i in 1:4
        @inbounds @simd for n in 1:(N-1)
            mstat[n+1, i] = l * (mstat[n, i] - xadj[n+1, i] * w[m]) + xadj[n+m+1, i] * w[1]
        end
    end

    # Unbias the estimates: variance -> standard deviation, then the dimensionless skew and excess kurtosis.
    @inbounds for n in 1:N
        var = mstat[n, 2] / B2
        sd  = sqrt(var)
        mstat[n, 2] = sd
        mstat[n, 3] = sd == 0 ? zero(F) : (mstat[n, 3] / B3) / (sd * var)
        mstat[n, 4] = sd == 0 ? zero(F) : ((mstat[n, 4] / (var * var) - 3 * A4) / B4) - 3
    end

    return (mstat)
end


"""
    sample_std(x)
Compute the "sample" standard deviation of a series, `x`.
(Also available, un-exported, as `Finance.std`; it is not exported under that
name to avoid a clash with `Statistics.std`.)

# Type Constraints
- `T <: Real`

# Arguments
- `x :: AbstractVector{T}` -- The series to work with.

# Input Contract
The inputs are assumed to satisfy the constraints below.

| Constraint     | Description                                    | 
|:--------------:|:---------------------------------------------- |
| `\\|x\\| > 1`  | The length of the series is greater than ``1``.|

# Return
`std::F` -- The sample standard deviation (`F = float(T)`).
"""
@noinline function sample_std(x::AbstractVector{T}) where {T<:Real}
    F = _ftype(T)
    N = length(x)
    N > 1 || throw(DomainError(N, "The length of the series must be > 1."))
    sd = zero(F)
    mn = zero(F)
    @simd for i in 1:N
        @inbounds mn += x[i]
    end
    mn /= N

    @inbounds @simd for i in 1:N
        sd += (x[i] - mn) * (x[i] - mn)
    end
    return (sqrt(sd / (N - 1)))
end

const std = sample_std

@noinline function WWsum(w::AbstractVector{T})::T where {T<:Real}
    WW = zero(T)
    m = length(w)
    @inbounds for i in 1:(m-1)
        wwi = w[i] * w[i]
        wwj = zero(T)
        @simd for j in (i+1):m
            wwj += w[j] * w[j]
        end
        WW += wwi * wwj
    end
    return (WW)
end



"""
    entropy_index(x; <keyword arguments>)

Computes a (Discounted) Binned Entropy Index.
This is the ratio of entropy of the binned distribution of `x` against
the entropy of the uniform distribution.
The vector `x` is first capped by the lower and upper quantiles; then
binned into `n` number of equal width bins. A distribution is formed 
from the bins and the entropy computed. If `λ` is not 1, then a discounted
entropy is computed. This is an exponentially based discounting of the 
bin distribution based
on their "freshness". In either event, the ratio of this entropy to 
the entropy of the corresponding uniform distribution (of `n` bins) is returned.

The discount is applied by *time*: the last (newest) value of `x` has weight ``1``,
the one before it ``\\lambda``, then ``\\lambda^2`` and so on; each value adds its
weight to the bin it falls in.

# Type Constraints
- `T <: Real`

# Arguments
- `x::AbstractVector{T}`          -- Vector to process.

# Keyword Arguments
- `n=10::Int`                     -- The number of bins.
- `tol=1.0/(100 * n)::Float64`    -- Error tolerance used with equivalency test of number to 0 or 1.
- `probs=[0.01, 0.99]`            -- Vector of quantile min and max.
- `λ=1.0::Float64`                -- Discount value.

# Input Contract
- `n > 2`
- `0 < tol < 0.01` 
- `|probs| == 2`
- ``0 < \\lambda \\le 1``
- `|x| > 0`

# Return
`::Float64` -- The (discounted) binned entropy index. A constant series (all values in one bin) has index `0`.
"""
function entropy_index(x::AbstractVector{T}         ;
                       n::Int=10                    ,
    				   tol::Float64=1.0 / (100 * n) ,
     				   probs::AbstractVector{<:Real}=[0.01, 0.99],
    				   λ::Float64=1.0                ) where {T<:Real}

    # Check Input contract.
    n > 2              || throw(DomainError(n,     "Bad number of bins."))
    0.0 < tol < 0.01   || throw(DomainError(tol,   "Bad tolerance value."))
    length(probs) == 2 || throw(DomainError(probs, "Bad quantile vector, must have length 2."))
    0.0 < λ <= 1.0     || throw(DomainError(λ,     "Bad discount parameter."))
    isempty(x)         && throw(DomainError(0,     "The series must not be empty."))

    # Get the data extrema for the quantile filtered data.
    qmin, qmax = Statistics.quantile(x, probs)

    # This will be the data distribution structure based on the granularity (`n`).
    F = Float64
	bdist = zeros(F, n)
    width = (qmax - qmin) / n

    # A (nearly) constant series: everything falls in one bin, whose entropy is 0.
    width > 0 || return zero(F)

    # For each filtered data point assign it to its bin index.
    idxs = floor.(Int, 1.0 .+ (div.(x .- qmin .- tol, width)))
    idxs .= min.(idxs, n)
    idxs .= max.(idxs, 1)

    # Increment all bins for each occurrence from the series discounted by
    # "freshness": starting from the end (newest value) of the time series.
    lm = 1.0
    @inbounds for j in Iterators.reverse(idxs)
        bdist[j] += lm
        lm *= λ
    end

    # Finish off the binned empirical distribution.
    bdist ./= sum(bdist)

    # Get the discounted entropy of the binned distribution.
	ent = zero(F)
    @inbounds @simd for i in 1:n
        prb = bdist[i]
        ent -= isapprox(prb, 0.0; atol=tol) ? 0.0 : prb * log(prb)
    end

    # Return the normalized discounted binned entropy.
    # Normalize by the entropy of the uniform distribution over `n` values.
    return ent / log(n)
end


"""
    pow_n(x, n)

Fast (non-negative) integer powers: ``x^n``.
Uses repeated squaring in combination with the bit vector
representation of `n`.

# Type Constraints
- `T <: Number`

# Arguments
- `x::T`     -- The base value.
- `n::Int` -- The power.

# Input Contract
- ``n \\ge 0`` 

# Return
`::T`        -- The Power Value.
"""
@noinline function pow_n(x::T, n::Int) where {T<:Number}

    # Check input contract.
    if n < 0
        throw(DomainError(n, "Parameter `n` must be non-negative."))
    end

    o = one(T)

    # Anything to the 0'th power is 1.
    if n == 0
        return (o)
    end

    # -- Do repeated squaring based on the digits of `n-1`. --
    # Initialize values.
    s = x
    n2d = digits(n - 1, base=2)

    # Repeated squaring.
    for d in n2d
        s *= d == 1 ? x : o
        x *= x
    end
    return s
end


"""
    pow_n(x, n, m)

Fast integer (non-negative) powers with modulus: ``x^n \\; {\\rm mod } \\; m`` (using `mod`, so the result has the sign of `m`).
Uses repeated squaring in combination with the bit vector
representation of `n`.

The output will be of type `T^* = typeof(promote(x, m))`.

# Type Constraints
- `T <: Real`
- `S <: Real`

# Arguments
- `x::T`     -- The base value.
- `n::Int` -- The power.
- `m::S`     -- The modulus.

# Input Contract
- ``n \\ge 0`` 

# Return
``::T^*``    -- The Power Value mod `m`.
"""
@noinline function pow_n(x::T, n::Int, m::S) where {T<:Real,S<:Real}

    # Check input contract.
    if n < 0
        throw(DomainError(n, "Parameter `n` must be non-negative."))
    end

    # Promote to a common type, this will be the type of the output.
    x, m, o = promote(x, m, Int8(1))

    # Anything to the 0'th power is 1.
    if n == 0
        return (o)
    end

    # Get modulus value (`mod`, so the result is in `[0, m)` for positive `m`).
    x = mod(x, m)

    # -- Do repeated squaring based on the digits of `n-1`. --
    # Initialize values.
    s = x
    n2d = digits(n - 1, base=2)

    # Repeated squaring.
    for d in n2d
        s *= d == 1 ? x : o
        s = mod(s, m)
        x *= x
        x = mod(x, m)
    end

    return s
end


"""
    ewt_mean(ts, xs, b, lm)
    

Computes the moving (exponential decayed) temporal average of the data `xs` over windows of length `b`.
Temporal averaging over a window means that the time stamps are differenced and we associate the difference
of time points ``t_{i}, t_{i+1}``, ``\\Delta_{i} = t_{i+1} - t_i``, with the data point ``x_i``.
The rationale: The data point ``x_i`` has been around since ``t_i`` until ``t_{i+1}``, so it should be should 
weight it (in an un-normalized way) by this distance.
The temporal decay will adjust the temporal weights by an exponential which puts more weight on recent data within the window.
This is also done in an un-normalized way. Then the weights are normalized and data, `xs`, is averaged within the window.

**NOTE:** Set `lm` to 1.0 to just have temporal weighting *without* decay.

# Arguments
- ts::AbstractVector{<:Real} -- Data time stamps -- ordered from smallest (oldest) to largest (newest).
- xs::AbstractVector{<:Real} -- Data values associated with time stamps.
- b::Int            -- The width of the window
- lm::Real          -- The decay factor: 0.0 < lm <= 1.0

# Input Contract
- |ts| == |xs|
- 0 < b < |xs|
- 0.0 < lm <= 1.0

# Return
::Vector{Float64} -- A vector of length |xs| - b.

"""
function ewt_mean(ts::AbstractVector{<:Real},
                  xs::AbstractVector{<:Real},
                  b::Int                    ,
                  lm::Real                   )
    n = length(ts)

    # Check input contract.
    n == length(xs) || throw(DomainError(n, "The length of the time and data series must match."))
    0 < b < n       || throw(DomainError(b, "The length of the moving window, `b`, must be in the interval (0, |xs|)."))
    0.0 < lm <= 1.0 || throw(DomainError(lm, "The decay factor, `lm`, must be in the interval (0.0, 1.0]."))

    # `wm` will be the weighted mean that is returned.
    wm = Vector{Float64}(undef, n - b)

    # Construct the temporal decay factors -- decay more as we go back in time.
    decayFs = Vector{Float64}(undef, b)
    decayFs[b] = 1.0
    @inbounds for k in (b-1):-1:1
        decayFs[k] = decayFs[k+1] * lm
    end

    # Get temporal weighting.
    dts = diff(ts)

    #= Average over all windows (length b -- bandwidth) the data, xs[i, i+b-1], 
       using the temporal difference weighting but modifying them by a factor
       which is based on how far back in time one goes within the band: data[i] * temporal_weight[i] * decay_factor[i]
       Here decay_factor looks like l^(b-1), l^(b-2), ... l^2, l, 1.
       Finally, we need to normalize these modified weights: ws[i] =  (temporal_weight[i] * decay_factor[i]) 
       so that within the band they sum to 1.
	=#

    @inbounds for i in 1:(n-b)
        num = 0.0
        den = 0.0
        @simd for j in 1:b
            w = decayFs[j] * dts[i+j-1]                # Modified (un-normalized) weights.
            num += w * xs[i+j-1]
            den += w
        end
        wm[i] = num / den                              # Since weights are not normalized, must divide by their sum.
    end

    return wm
end


end # module Finance

