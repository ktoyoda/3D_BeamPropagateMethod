using JLD2
using Plots
plotting_t = 1:10
z = 1400
struct data
    E::Array{ComplexF64,3}
    N::Array{Float64,3}
end

function data_arra(t)
    strt = string(t)
    jldopen("xyz1105um_1s_"*strt*"_E.jld2","r") do file1 
    jldopen("xyz1105um_1s_"*strt*"_N.jld2","r") do file2 
        E1 = file1[strt*"/E"]   
        N2 = file2[strt*"/N"]
        return data(E1,N2)
    end
    end
end

for t in plotting_t
    strt = string(t)

    Exy = view(data_arra(t).E,:,:,z)
    Ixy = abs.(Exy)
    Nxy = view(data_arra(t).N,:,:,z)
    x = range(-100,100,length = 100)
    y = range(-100,100,length = 100)
    p1 = heatmap(x, y, Ixy,levels = 200,xlabel = "Z [mm]", ylabel = "X [μm]" )
    strt = string(t)
    display(p1)
    p2 = heatmap(x*10^6, y*10^6, Nxy,levels = 200,xlabel = "Z [mm]", ylabel = "X [μm]" )
    strt = string(t)
    display(p1)
end