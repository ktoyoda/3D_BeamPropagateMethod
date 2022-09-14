using JLD2
using Plots
plotting_t = [4 ]
@show plotting_t
plotting_z = 1:10
struct data
    E::Array{ComplexF64,3}
    N::Array{Float64,3}
end

function data_arra(t)
    strt = string(t)
    jldopen("1105_532_/xyz1105um_1s_"*strt*"_E.jld2","r") do file1 
    jldopen("1105_532_/xyz1105um_1s_"*strt*"_N.jld2","r") do file2 
        E1 = file1[strt*"/E"]   
        N2 = file2[strt*"/N"]
        return data(E1,N2)
    end
    end
end



for t in plotting_t
    for z in plotting_z
        Exy = view(data_arra(t).E,:,:,(z*100-99))
        Ixy = abs.(Exy)
        Nxy = view(data_arra(t).N,:,:,(z*100-99))
        
        x = range(-100,100,length = size(Exy,1))
        y = range(-100,100,length = size(Exy,1))

        
        strt = string(t)
        p1 = heatmap(x, y, Ixy,levels = 200,xlabel = "Z [mm]", ylabel = "X [μm]", aspect_ratio = 1)
        strt = string(t)
        display(p1)
    end
end


for t in plotting_t
    for z in plotting_z
        Exy = view(data_arra(t).E,:,:,(z*100-99))
        Ixy = abs.(Exy)
        Nxy = view(data_arra(t).N,:,:,(z*100-99))
        
        x = range(-100,100,length = size(Exy,1))
        y = range(-100,100,length = size(Exy,1))

        strt = string(t)
        p2 = heatmap(x*10^6, y*10^6, Nxy,levels = 200,xlabel = "Z [mm]", ylabel = "X [μm]" ,clim = (1.47,1.5), aspect_ratio = 1)
        strt = string(t)
        display(p2)
    end
end
