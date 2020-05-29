#=
Julia program to compute applying Gaussian filter to 2D function
=#

"""
Module for parameters & variables
"""
module ParamVar
    """
    Global parameters
    """
    struct Parameters
        σ::Float64  # Kernel width
        x_lim::Float64  # Boundary of region [-x_lim, x_lim]
        dx::Float64  # Spatial resolution
        x::Array{Float64}  # coordinate
        N::Int64  # Number of points in one direction
    end

    """
    Property of function
    """
    mutable struct Function
        int_G::Array{Float64}  # Integral of Gaussian kernel
        f::Array{Float64}  # Function
        f_G::Array{Float64}  # Gaussian filtered function
    end
end


"""
Module for computation
"""
module Computation

    """
    Compute relative distance of two points
    """
    function compute_relative_distance_square(point_1, point_2)
        return distance_square = (point_1 - point_2) ^ 2
    end


    """
    Compute Gaussian weight
    1D: G(x,x',σ²) = 1/[(2π)^1/2 σ] exp{ 1/(2σ²) (x'-x)² }
    3D: G(x,x',σ²) = 1/[(2π)^3/2 σ³] exp{ 1/(2σ²) [(x'-x)² + (y'-y)² + (z'-z)²] }
    """
    function compute_gaussian_weight(σ, r²)
        return 1.0 / (sqrt(2.0*π) * σ) * exp(r² / (2.0*σ^2))
    end


    """
    Compute integral of Gaussian kernel
    """
    function compute_gaussian_kernel_integral(param, kernel_integral)

        for itr_x_out = 1:param.N

            tmp = 0.0
            for itr_x_in = 1:param.N
                # Compute square of relative distance
                r² = compute_relative_distance_square(param.x[itr_x_out], param.x[itr_x_in])

                # Sum up Gaussian weight
                tmp += compute_gaussian_weight(param.σ, r²)
            end

            kernel_integral[itr_x_out] = tmp
        end

    end


end


# ========================================
# Main function
# ========================================

## Declare modules
using Printf
using Plots

using .ParamVar
using .Computation:
    compute_gaussian_kernel_integral


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
σ = π/5.0
x_lim = π
dx = 0.01
x = -x_lim:dx:x_lim
N = length(x)
param = ParamVar.Parameters(
    σ, x_lim, dx, x, N
)

int_G = Array{Float64}(undef, param.N)
f = Array{Float64}(undef, param.N)
f_G = Array{Float64}(undef, param.N)

func = ParamVar.Function(
    int_G, f, f_G
)


# ----------------------------------------
## Compute integral of Gaussian kernel
# ----------------------------------------
compute_gaussian_kernel_integral(param, func.int_G)

# ----------------------------------------
## Set function
# ----------------------------------------
set_function(param, func.f)

# ----------------------------------------
## Apply Gaussian filter to the function
# ----------------------------------------
apply_Gaussian(param, func)

# ----------------------------------------
## Output result
# ----------------------------------------
out_functions(param, func)
