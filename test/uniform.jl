using GeneralizedSampling
using Base.Test

M = 10

# 1D
begin
	xi = float(collect(-M:M))
	@test isuniform(xi)

	@test !isuniform(2*rand(M)-1)
end


# 2D
begin
	xi = grid(M)
	@test isuniform(xi)

	@test !isuniform(2*rand(M,2)-1)
end

