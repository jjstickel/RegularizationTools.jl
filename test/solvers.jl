r = mdopen("shaw", 100, false)
A, x = r.A, r.x
rn = CSV.read("random_numbers.csv", DataFrame)
y = A  * x
b = y + 0.1y .* rn.x
x₀ = 0.5x

Ψ = setupRegularizationProblem(A,1)
λopt = @> solve(Ψ, b, alg=:gcv_svd) getfield(:λ)
@test round(λopt, digits = 3) ≈ 2.396
λopt = @> solve(Ψ, b, alg=:gcv_tr) getfield(:λ)
@test round(λopt, digits = 3) ≈ 2.396
λopt = @> solve(Ψ, b, alg =:L_curve) getfield(:λ)
@test round(λopt, digits = 3) ≈ 1.681
λopt = @> solve(Ψ, b, x₀, alg =:L_curve) getfield(:λ)
@test round(λopt, digits = 3) ≈ 1.725
