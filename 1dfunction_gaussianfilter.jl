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
        return exp(-r² / (2.0*σ^2)) / (sqrt(2.0*π) * σ)
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


    """
    Apply Gaussian filtering
    """
    function apply_gaussian(param, func)

        for itr_x_out = 1:param.N

            tmp = 0.0
            for itr_x_in = 1:param.N
                # Compute square of relative distance
                r² = compute_relative_distance_square(param.x[itr_x_out], param.x[itr_x_in])

                # Sum up contribution of Gaussian weight in the nearby of original point
                tmp += compute_gaussian_weight(param.σ, r²) * func.f[itr_x_in]
            end

            func.f_G[itr_x_out] = tmp / func.int_G[itr_x_out]
        end
    end


    """
    Set test function
    """
    function set_function(param, test_function)
        test_function .= sin.(param.x) + cos.(3.0 * param.x) + sin.(5.0 * param.x)
    end
end


"""
Module for output figures
"""
module PlotFigures
    using Printf
    using PyPlot
    using PyCall
    sns = pyimport("seaborn")
    sns.set(
        context="talk",
        style="white",
        palette="plasma",
        font="sans-serif",
        font_scale=1,
        color_codes=false,
    )
    rc("text", usetex ="true")


    """
    Plot function
    """
    function out_functions(param, func)
        # Figure setting
        pygui(true)
        fig = figure()
        ax = gca(
            xlabel=L"$x$", ylabel=L"$f$",
            xlim=[-param.x_lim, param.x_lim],
            xticks=[-π, -π/2.0, 0.0, π/2.0, π],
            xticklabels=[L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"],
        )

        # Plot original/filtered function
        ax.plot(
            param.x, func.f,  # original
            param.x, func.f_G,  # filtered
        )

        ax.legend(["original", "filtered"])

        # Save figure
        filename = @sprintf("./fig/func_σ=%.3f.png", param.σ)
        savefig(
            filename,
            bbox_inches="tight", pad_inches=0.1
        )
    end


    """
    Plot Gaussian function
    """
    function out_gaussian(param)
        # Figure setting
        pygui(true)
        fig = figure()
        ax = gca(
            xlabel=L"$x$", ylabel=L"$f$",
            xlim=[-param.x_lim, param.x_lim],
            xticks=[-π, -π/2.0, 0.0, π/2.0, π],
            xticklabels=[L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"],
        )

        # Set Gaussian function
        gaussian = Array{Float64}(undef, param.N)
        gaussian .= exp.(-param.x.^2 / (2.0*param.σ^2)) / (sqrt(2.0*π) * param.σ)

        # Plot Gaussian function
        ax.plot(param.x, gaussian)

        ax.legend(["Gaussian"])

        # Save figure
        filename = @sprintf("./fig/gaussian_σ=%.3f.png", param.σ)
        savefig(
            filename,
            bbox_inches="tight", pad_inches=0.1
        )
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
    compute_gaussian_kernel_integral,
    apply_gaussian,
    set_function
using .PlotFigures:
    out_functions,
    out_gaussian


# ----------------------------------------
## Declare parameters & mutable structs
# ----------------------------------------
σ = π/3.0
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
out_gaussian(param)

# ----------------------------------------
## Set function
# ----------------------------------------
set_function(param, func.f)

# ----------------------------------------
## Apply Gaussian filter to the function
# ----------------------------------------
apply_gaussian(param, func)

# ----------------------------------------
## Output result
# ----------------------------------------
out_functions(param, func)
