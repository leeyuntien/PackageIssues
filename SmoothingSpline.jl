abstract type SmoothingSpline end

function estParams(data)
    fn = generateSmoothingSplineFunction(data, lambda=1000)
    generateSplinePoints(fn, data)
end

function generateSmoothingSplineFunction(data, lambda=1000)
    betas = solveForBetas(data, lambda)
    return x => 
    begin
        basis = createBasisArray(x, data)
        return betas.reduce((acc, beta, i) => (acc += beta * basis[i]), 0)
    end
end

function createBasisArray(x, data)
    [
  1,
  x,
  x ** 2,
  x ** 3,
  ...data.x.map((x_k) => pos(x - x_k) ** 3),
]
end

function generateSplinePoints(fn, data)
    step = (max(data) - min(data)) / 1000
    xs = [min(data) + i * step for i in 1:1000]
    ys = fn.(xs)
end

function solveForBetas(data, lambda)
    x = createBasisArray(x, data.x)
    y = data.y

    U, sigma, V = x.SVD

    innerInverse = sigma * sigma + lambda == 0 ? 0 : 1 / (sigma * sigma + lambda)

    betas = V * innerInverse * diag(sigma) * U' * y
end 