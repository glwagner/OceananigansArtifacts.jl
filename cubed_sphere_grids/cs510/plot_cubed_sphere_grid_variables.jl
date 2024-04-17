using JLD2
using PyPlot

function quick_and_dirty_plot(var, var_name)
    PyPlot.pcolormesh(var)
    PyPlot.title(var_name)
    PyPlot.colorbar()
    PyPlot.savefig("$var_name.png", dpi=300)
    PyPlot.close("all")
    @info "Saved $var_name.png"
end

face = 1 # just plot one face for now

jldopen("cubed_sphere_510_grid.jld2", "r") do file
    for key in keys(file["face$face"])
        quick_and_dirty_plot(file["face$face/$key"], key)
    end
end
