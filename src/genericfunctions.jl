@doc raw"""
    designmatrix(s::Any, q::Any, f::Function)::AbstractMatrix

- s is an array of setpoints
- q is an array of query points
- f is a function that maps the functor [Domain](@ref) to a solution y

The function to creates an array of domain nodes using the standard basis eᵢ of 
the vector space. The setpoint or query point can be any abstract notion of the 
input. For examples, a numerical value corresponding to a setting, a string label 
(["bin A", "bin B", ...), or a list of pixel coordinate [(1,1), (1,2), ...].  The
function f must accept a single argument of type [Domain](@ref) and provide a 
numerical mapping between input ``x`` and output ``y`` at query point ``q``. The
function designmatrix then returns a design matrix ``\mathbf{A}`` such that 

``
y = \mathbf{A}x
``

where x is an array of numerical input values. In the case that [q] = [s], the 
shortcut ```designmatrix(s, f)``` can be used.
"""
function designmatrix(s::Any, q::Any, f::Function)::AbstractMatrix
    n = length(s)
    x(i) = @_ map(_ == i ? 1.0 : 0.0, 1:n)
    nodes = [Domain(s, x(i), q[j]) for i ∈ 1:n, j ∈ 1:length(q)]
    return @> map(f, nodes) transpose copy
end

designmatrix(s::Any, f::Function) = designmatrix(s, s, f)

@doc raw"""
    forwardmodel(
        s::AbstractVector{T1},
        x::AbstractVector{T2},
        f::Function,
    )::AbstractArray where {T1<:Any,T2<:Number}

Forward model that maps numerical values [x] corresponding to 
setpoints [s] to output [y], given the function f. The function f must 
accept a single argument of type [Domain](@ref) and provide a 
numerical mapping between input ``x`` and output ``y`` at query point ``q``.

Note that the [designmatrix](@ref) and forward model are related

```
A = designmatrix(s, q, f)
y1 = A*x

y2 = forwardmodel(s, x, q, f)
y1 == y2
```
"""
function forwardmodel(
    s::AbstractVector{T1},
    x::AbstractVector{T2},
    q::AbstractVector{T3},
    f::Function,
)::AbstractArray where {T1<:Any,T2<:Number,T3<:Any}
    return @>> (@_ map(Domain(s, x, _), q)) map(f)
end
