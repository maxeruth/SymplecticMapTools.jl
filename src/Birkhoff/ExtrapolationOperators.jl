## Hankel stuff
# Create Block Hankel matrix with first column and final row given by ac and ar, i.e.
#      _                                              _
#     |ac[:,1]            | ac[:,2] | ... |___________|
#     |___________________|_________| ... |___________|
#     |ac[:,2]            | ac[:,3] | ... |           |
#     |___________________|_________| ... |___________|
# A = |...
#     |___________________|_________| ... |___________|
#     |ac[:,end-1]
#     |___________________|_________| ... |___________|
#     |ac[:, end]=ar[:,1] | ar[:,2] | ... | ar[:,end] |
#     |_                                             _|
function BlockHankelMatrix(ac::AbstractArray, ar::AbstractArray)
    @assert ac[:,end] == ar[:,1]

    a = vcat(vec(ac), vec(ar[:,2:end]));

    d, N1 = size(ac);
    d, N2 = size(ar);

    A = zeros(d*N1, N2);
    for ii = 1:N2
        A[:, ii] = a[(1:d*N1).+(d*(ii-1))];
    end
    A
end

struct BlockHankelPlan
    # Problem size
    d::Integer
    N1::Integer
    N2::Integer

    # Forward multiplication plans
    F
    Finv

    # Transpose multiplication plans
    FT
    FTinv

    # Hankel coefficients
    ahat::AbstractArray

    # Buffer arrays
    xbuf::AbstractVector
    ybuf::AbstractArray

    function BlockHankelPlan(ac, ar)
        d, N1 = size(ac);
        d, N2 = size(ar);
        N = nextprod((2,3), N1+N2-1)
        xbuf = zeros(N)
        ybuf = zeros(d, N);

        a = hcat(zeros(d, N-N1-N2+1), ar[:, end:-1:1], ac[:, end-1:-1:1]);

        F = plan_rfft(xbuf)
        ahat = rfft(a, 2);
        Finv = plan_irfft(ahat, N, 2)
        FT = plan_rfft(ybuf, 2);
        FTinv = plan_irfft(rfft(xbuf), N)
        new(d, N1, N2, F, Finv, FT, FTinv, ahat, xbuf, ybuf)
    end
end

function fast_block_hankel_multiply(a_plan::BlockHankelPlan, x; transpose = false)
    N1 = a_plan.N1;
    N2 = a_plan.N2;
    d = a_plan.d;
    F = a_plan.F;
    Finv = a_plan.Finv;
    FT = a_plan.FT;
    FTinv = a_plan.FTinv;
    ahat = a_plan.ahat;
    xbuf = a_plan.xbuf;
    ybuf = a_plan.ybuf;
    N = length(xbuf);

    if transpose
        ybuf[:, end:-1:end-N1+1] = reshape(x, d, N1);
        yhat = FT*ybuf
        # ayhat = [ahat[:,ii]'*yhat[:,ii] for ii = 1:size(yhat,2)];
        # ahatconj = conj(ahat);
        # println("ahatconj = $(size(ahatconj)), yhat = $(size(yhat))")
        ayhat = vec(ones(d)'*(conj.(ahat) .* yhat));

        x = FTinv*ayhat
        return x[1:N2]
    else
        xbuf[1:N2] = x;
        y = Finv*(ahat * Diagonal(F*xbuf))
        y = y[end:-1:1, :];
        return y[end:-1:end-d*N1+1]
    end
end



# Make a block Hankel linear operator with fast multiplication. Inspired by
# ToeplitzMatrices.jl. ac is an array of the coefficients for the first column,
# ar is a vector of coefficients for the final row. We require that
# ac[end] == ar[1].
#
# The operator can (mostly) be used as a normal matrix, i.e. the operations A*b
# and A'*b, will work as expected. For performing linear solves, it is recommended
# to use lsqr from IterativeSolvers.jl, but it is not implemented right now.
#
# See BlockHankelMatrix for a schematic
function block_hankel_linear_operator(ac, ar)
    plan = BlockHankelPlan(ac, ar);

    function block_hankel_prod!(res, v)
        res[:] = fast_block_hankel_multiply(plan, v; transpose = false);
    end

    function block_hankel_tprod!(res, v)
        res[:] = fast_block_hankel_multiply(plan, v; transpose = true);
    end

    LinearOperator(Float64, plan.d*plan.N1, plan.N2, false, false, block_hankel_prod!,block_hankel_tprod!)
end

##  The "P" matrix (imposes time reversal symmetry)
function mpe_p(K::Integer)
    M = K+1;
    N = 2K+1;
    function P_mul!(res, v, α, β::T) where T
        if β == zero(T)
            res[1] = v[M]
            for ii = 2:M
                res[ii] = (α)*(v[M-ii+1]+v[M+ii-1])
            end
        else
            res[1] = β*res[1] + v[M]
            for ii = 2:M
                res[ii] = β*res[ii] + (α)*(v[M-ii+1]+v[M+ii-1])
            end
        end
    end

    function P_tmul!(v, res, α, β::T) where T
        if β == zero(T)
            v[M] = α*res[1]
            for ii = 2:M
                v[M-ii+1] = (α)*res[ii]
                v[M+ii-1] = v[M-ii+1]
            end
        else
            v[M] = β*v[M] + α*res[1]
            for ii = 2:M
                x = (α)*res[ii]
                v[M-ii+1] = β*v[M-ii+1] + x
                v[M+ii-1] = β*v[M+ii-1] + x
            end
        end
    end

    LinearOperator(Float64, M, N, false, false, P_mul!, P_tmul!)
end

# Used for when c_0 is constrained rather than c_(±K-1)
function mpe_p_2(K::Integer)
    M = K;
    N = 2K+1;
    mid = K+1
    function P_mul!(res, v, α, β::T) where T
        if β == zero(T)
            for ii = 1:M
                res[ii] = (α)*(v[mid-ii]+v[mid+ii])
            end
        else
            for ii = 1:M
                res[ii] = β*res[ii] + (α)*(v[mid-ii]+v[mid+ii])
            end
        end
    end

    function P_tmul!(v, res, α, β::T) where T
        if β == zero(T)
            for ii = 1:M
                v[mid-ii] = (α)*res[ii]
                v[mid+ii] = v[mid-ii]
            end
        else
            for ii = 1:M
                x = (α)*res[ii]
                v[mid-ii] = β*v[mid-ii] + x
                v[mid+ii] = β*v[mid+ii] + x
            end
        end
    end

    LinearOperator(Float64, M, N, false, false, P_mul!, P_tmul!)
end

# """
#     rre_p(K::Integer)
#
# Function used for iterative solutions of the RRE problem. The matrix, when
# applied to a vector of
# """
function rre_p(K::Integer)
    M = K;
    N = 2K+1;
    mid = K+1
    function P_mul!(res, v, α, β::T) where T
        if β == zero(T)
            for ii = 1:M
                res[ii] = α*(v[mid-ii]+v[mid+ii]-v[mid-ii+1]-v[mid+ii-1])
            end
        else
            for ii = 1:M
                res[ii] = β*res[ii] + α*(v[mid-ii]+v[mid+ii]-v[mid-ii+1]-v[mid+ii-1])
            end
        end
    end

    function P_tmul!(v, res, α, β::T) where T
        if β == zero(T)
            v[mid] = 0.
            for ii = 1:M
                αr = α*res[ii]
                v[mid - ii + 1] -= αr
                v[mid + ii - 1] -= αr
                v[mid + ii    ] = αr
                v[mid - ii    ] = αr
            end
        else
            for ii = 1:M
                αr = α*res[ii]
                v[mid - ii + 1] -= αr
                v[mid + ii - 1] -= αr
                v[mid + ii    ] = β*v[mid + ii] + αr
                v[mid - ii    ] = β*v[mid - ii] + αr
            end
        end
    end

    LinearOperator(Float64, M, N, false, false, P_mul!, P_tmul!)
end
