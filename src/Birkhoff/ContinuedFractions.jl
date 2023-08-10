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

function big_cont_frac_eval(c::ContFrac)
    big_cont_frac_eval(c.a)
end


function eval(c::ContFrac)
    a = c.a;
    return cont_frac_eval(a)
end

function ContFrac(x::Number; tol = 1e-15)
    N = 5000;
    a = ones(Int64, N);
    x = mod(x,1);
    x0 = x;
    for ii = 1:N
        x = 1/x;
        a[ii] = floor(Int64, x);
        x = mod(x,1);
        
        # println(a[ii])
        approx = (typeof(x0) <: Float64) ? cont_frac_eval(a[1:ii]) : big_cont_frac_eval(a[1:ii])
        if abs(x0 - approx) < tol
            return ContFrac(a[1:ii])
        end
        # println("typeof(x0) = $(typeof(x0)), typeof(approx) = $(typeof(approx)), err = $(abs(x0 - approx))")
    end
    
    return nothing
end

function partial_frac(c::ContFrac, n::Integer)
    x = 0 // 1;
    a = c.a;
    
    for ii = n:-1:1
        x = 1 // (x + a[ii]);
    end
    
    return x;
end

function denoms(c::ContFrac)
    a = c.a;
    N = length(a);
    
    D = zeros(Int64, N);
    for ii = 1:N
        x = partial_frac(c, ii);
        D[ii] = denominator(x);
    end
    
    D
end

function big_partial_frac(c::ContFrac, n::Integer)
    x = BigInt(0) // BigInt(1);
    a = c.a;
    
    for ii = n:-1:1
        x = BigInt(1) // (x + BigInt(a[ii]));
    end
    
    return x;
end

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