import Pkg
Pkg.activate(".")
Pkg.instantiate()


using DelimitedFiles
using LinearAlgebra
using Printf
using SparseArrays

using BandedMatrices
using DataStructures: top
using Interpolations
using JSON
using OrdinaryDiffEq, DiffEqCallbacks
using Parameters
using RecursiveArrayTools
using Roots
using StaticArrays
using SummationByPartsOperators

import FFTW; FFTW.set_num_threads(1)

using PyCall, LaTeXStrings; import PyPlot; plt=PyPlot
inset_locator = pyimport("mpl_toolkits.axes_grid.inset_locator")

cycler = pyimport("cycler").cycler
line_cycler   = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                 cycler(marker=["4", "2", "3", "1", "+", "x", "."]))

plt.rc("axes", prop_cycle=line_cycler)
plt.rc("text", usetex=true)
plt.rc("text.latex", preamble="\\usepackage{newpxtext}\\usepackage{newpxmath}\\usepackage{commath}\\usepackage{mathtools}")
plt.rc("font", family="serif", size=18.)
plt.rc("savefig", dpi=100)
plt.rc("legend", loc="best", fontsize="medium", fancybox=true, framealpha=0.5)
plt.rc("lines", linewidth=2.5, markersize=10, markeredgewidth=2.5)


function convergence_orders(Ns, errors)
    @assert length(Ns) == length(errors)
    
    orders = similar(errors, length(errors)-1)
    orders[1] = NaN
    for i in 2:length(errors)
        orders[i-1] = - log(errors[i-1] / errors[i]) / log(Ns[i-1] / Ns[i])
    end
    orders
end

function linear_regression(x, y)
    A = hcat(one.(x), x)
    param = A \ y
    param
end

function even_values(val_N)
    res = Vector{eltype(val_N)}()
    for N in val_N
        N1 = iseven(N) ? N : N-1
        push!(res, N1)
    end
    res
end

function odd_values(val_N)
    res = Vector{eltype(val_N)}()
    for N in val_N
        N1 = iseven(N) ? N-1 : N
        push!(res, N1)
    end
    res
end

function evenodd_values(val_N)
    res = Vector{eltype(val_N)}()
    for N in val_N
        N1, N2 = iseven(N) ? (N, N+1) : (N-1, N)
        push!(res, N1)
        push!(res, N2)
    end
    res
end
