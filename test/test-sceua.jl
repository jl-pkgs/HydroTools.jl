using HydroTools
using Test

kw = (; maxn=1e4, kstop=10, pcento=0.1, peps=0.001, iseed=1, iniflg=0)

@testset "functn1" begin
  function functn1(x)
    # This is the Goldstein - Price Function
    # Bound X1 = [-2, 2], X2 = [-2, 2]
    # Global Optimum:3.0, (0.0, -1.0)
    x1 = x[1]
    x2 = x[2]
    u1 = (x1 + x2 + 1.0)^2
    u2 = 19 .- 14 .* x1 + 3 .* x1^2 - 14 .* x2 + 6 .* x1 * x2 + 3 .* x2^2
    u3 = (2 .* x1 - 3 .* x2)^2
    u4 = 18 .- 32 .* x1 + 12 .* x1^2 + 48 .* x2 - 36 .* x1 * x2 + 27 .* x2^2
    u5 = u1 * u2
    u6 = u3 * u4
    (1 .+ u5) * (30 .+ u6)
  end

  bl = [-2, -2]
  bu = [2, 2]
  x0 = [1, 1]

  x, feval, exitflag = sceua(functn1, x0, bl, bu; maxn=2000)
  @test abs(feval - 3) <= 1e-6
end

@testset "functn2" begin
  # !存在错误
  function functn2(x)
    # This is the Rosenbrock Function
    # Bound:X1 = [-5, 5], X2 = [-2, 8]
    # Global Optimum:0, (1, 1)
    x1 = x[1]
    x2 = x[2]
    a = 100
    a * (x2 - x1^2)^2 + (1 - x1)^2
  end

  bl = [-5, -5]
  bu = [5, 5]
  x0 = [-1, 1]
  x, feval, exitflag = sceua(functn2, x0, bl, bu; kw...)
  # @show x, feval
  @test abs(feval - 0) <= 1e-6
end

@testset "functn3" begin
  function functn3(x)
    # This is the Six - hump Camelback Function.
    # Bound:X1 = [-5, 5], X2 = [-5, 5]
    # True Optima:-1.031628453489877, (-0.08983, 0.7126), (0.08983, -0.7126)
    x1 = x[1]
    x2 = x[2]
    (4 - 2.1 * x1^2 + x1^4 / 3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2
  end
  bl = [-5, -2]
  bu = [5, 8]
  x0 = [-0.08983, 0.7126]
  x, feval, exitflag = sceua(functn3, x0, bl, bu)

  @test abs(feval - (-1.031628453489877)) <= 1e-6
end

@testset "functn4" begin
  # !存在错误
  function functn4(x)
    # This is the Rastrigin Function
    # Bound:X1 = [-1, 1], X2 = [-1, 1]
    # Global Optimum:-2, (0, 0)
    x1 = x[1]
    x2 = x[2]
    x1^2 + x2^2 - cos(18.0 * x1) - cos(18.0 * x2)
  end

  bl = [-1, -1]
  bu = [1, 1]
  x0 = [-1, -1]

  x, feval, exitflag = sceua(functn4, x0, bl, bu; kw...)
  @test abs(feval - (-2)) <= 1e-3
end

@testset "functn5" begin
  # !存在错误
  function functn5(x)
    # This is the Griewank Function [2 - D | 10 - D]
    # Bound:X[i] = [-600, 600], for i = 1, 2, ... , 10
    # Global Optimum:0; at origin
    nopt = length(x)
    d = nopt == 2 ? 200 : 4000

    u1 = 0.0
    u2 = 1.0

    for j = 1:nopt
      u1 = u1 + x[j]^2 / d
      u2 = u2 * cos(x[j] / sqrt(j))
    end
    u1 - u2 + 1
  end
  bl = -600 * ones(10)
  bu = 600 * ones(10)
  x0 = ones(10) * -1.0
  x, feval, exitflag = sceua(functn5, x0, bl, bu; kw...)
  @test abs(feval - 0) <= 1e-6
end


@testset "functn6" begin
  function functn6(x)
    # This is the Shekel Function
    # Bound:X[i] = [0, 10], j = 1, 2, 3, 4
    # Global Optimum:-10.5364098252, (4, 4, 4, 4)
    # Data for Skekel function coefficients (n = 4, m = 10)
    a1 = [4.0 1.0 8.0 6.0 3.0 2.0 5.0 8.0 6.0 7.0
      4.0 1.0 8.0 6.0 7.0 9.0 5.0 1.0 2.0 3.6
      4.0 1.0 8.0 6.0 3.0 2.0 3.0 8.0 6.0 7.0
      4.0 1.0 8.0 6.0 7.0 9.0 3.0 1.0 2.0 3.6]
    c1 = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5]
    nopt = length(x)

    f = 0.0
    for i = 1:10
      u = 0.0
      for j = 1:nopt
        u = u + (x[j] - a1[j, i])^2
      end

      u = 1.0 / (u + c1[i])
      f = f - u
    end
    f
  end

  bl = zeros(4)
  bu = 10 * ones(4)
  x0 = [4, 4, 4, 3.0]

  x, feval, exitflag = sceua(functn6, x0, bl, bu; maxn=1e4)
  @test abs(feval - -10.5364098252) <= 1e-5
end

