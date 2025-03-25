module mean_intensity

# Линкер найти надо
# using Base.Libdl

using DataFrames, StatsBase

function integrate_mean_intensity(I::Vector{T}, L::T, N::Int) where T<: Real
    dx = L / N;
    X = [-L/2 + dx * i for j in 1:N for i in 1:N]
    Y = [-L/2 + dx * j for j in 1:N for i in 1:N]
    R = @. sqrt(X^2 + Y^2)
    
    data = DataFrame(x = X, y = Y, r = R, I = I)

    center = (0.0, 0.0)
    N_bins = N / 2
    dr = L / N_bins
    r_edges = Vector(0:dr:maximum(data.r))
    maximum(data.r) > last(r_edges) && push!(r_edges, last(r_edges) + dr)
    data.r_bin = @. floor(Int, data.r / dr) * dr + dr/2
    data_grouped = groupby(data, :r_bin)
    result = combine(data_grouped, :I => mean => :I_mean)

    I0 = similar(I)
    for i in eachindex(I)
        I0[i] = result.I_mean[isapprox.(result.r_bin, R[i]; atol = dr)] |> first
    end

    return I0
end

const integrate_mean_intensity = Base.@cfunction(
    integrate_mean_intensity, 
    Ptr{Cdouble}, 
    (Ptr{Cdouble}, Ref{Cdouble}, Ref{Cint})
)

end