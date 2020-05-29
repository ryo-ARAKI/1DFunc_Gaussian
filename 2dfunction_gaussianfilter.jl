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
        x_lim::Float64  # Boundary of region [-x_lim, x_lim]
        dx::Float64  # Spatial resolution
        x::Array{Float64}  # coordinate
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



# ========================================
# Main function
# ========================================

## Declare modules
using Printf
using Plots

using .ParamVar


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
x_lim = π
dx = 0.01
x = -x_lim:dx:x_lim
param = ParamVar.Parameters(
    x_lim, dx, x
)

int_G = Array{Float64}(undef, length(x))
f = Array{Float64}(undef, length(x))
f_G = Array{Float64}(undef, length(x))

func = ParamVar.Function(
    int_G, f, f_G
)


# ----------------------------------------
## Compute integral of Gaussian kernel
# ----------------------------------------
compute_Gaussian_kernel(param, func.int_G)

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
