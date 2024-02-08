"""
    ContFrac(a::AbstractArray)

A continued fraction in (0,1]. The member `a` is a vector of the continued
fraction coefficients so that ω = 1/(a[1] + 1/(a[2] + ...))
"""
struct ContFrac
    a::AbstractArray
end

function cont_frac_eval(a::AbstractArray)
    Na = length(a);

    x = 0;
    for ii = Na:-1:1
        x = 1/(x + a[ii]);
    end

    return x
end

function big_cont_frac_eval(a::AbstractArray)
    Na = length(a);

    x = BigFloat(0.);
    for ii = Na:-1:1
        x = 1/(x + a[ii]);
    end

    return x
end

"""
    big_cont_frac_eval(c::ContFrac)

Find a high precision representation of the continued fraction as a `BigFloat`. To
set the precision, use `precision`.
"""
function big_cont_frac_eval(c::ContFrac)
    big_cont_frac_eval(c.a)
end

"""
    evaluate(c::ContFrac)

Find a floating point representation of the continued fraction.
"""
function evaluate(c::ContFrac)
    a = c.a;
    return cont_frac_eval(a)
end

"""
    ContFrac(x::Number; tol = 1e-15)

Find a continued fraction representation of `x` to tolerance `tol`.
"""
function ContFrac(x::Number; tol = 1e-15)
    N = 5000;
    a = ones(Int, N);
    x = mod(x,1);
    x0 = x;
    for ii = 1:N
        if x == zero(x)
            a[ii] = typemax(Int)
            x = zero(x)
        else
            x = 1/x;
            a[ii] = floor(Int, x);
            x = mod(x,1);
        end

        approx = (typeof(x0) <: Float64) ? cont_frac_eval(a[1:ii]) : big_cont_frac_eval(a[1:ii])
        if abs(x0 - approx) < tol
            return ContFrac(a[1:ii])
        end
    end

    return nothing
end

"""
    partial_frac(c::ContFrac, n::Integer)

Find the `n`th convergent of the continued fraction `c` as a rational number.
"""
function partial_frac(c::ContFrac, n::Integer)
    x = 0 // 1;
    a = c.a;

    for ii = n:-1:1
        x = 1 // (x + a[ii]);
    end

    return x;
end

"""
    denoms(c::ContFrac)

Return a list of the denominators of the continued fractions of `c`
"""
function denoms(c::ContFrac)
    a = c.a;
    N = length(a);

    D = zeros(Int, N);
    for ii = 1:N
        x = partial_frac(c, ii);
        D[ii] = denominator(x);
    end

    D
end

"""
    big_partial_frac(c::ContFrac, n::Integer)

Find the `n`th convergent of the continued fraction `c` as a rational number as
 `BigInt`s. Set the precision with `precision`.
"""
function big_partial_frac(c::ContFrac, n::Integer)
    x = BigInt(0) // BigInt(1);
    a = c.a;

    for ii = n:-1:1
        x = BigInt(1) // (x + BigInt(a[ii]));
    end

    return x;
end

"""
    big_denoms(c::ContFrac)

Return a list of the denominators of the continued fractions of `c` as
`BigInt`s. Set the precision with `precision`.
"""
function big_denoms(c::ContFrac)
    a = c.a;
    N = length(a);

    D = zeros(BigInt, N);
    for ii = 1:N
        x = big_partial_frac(c, ii);
        D[ii] = denominator(x);
    end

    D
end
