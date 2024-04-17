using JLD2

# Each file (big-endian, 32 bits) contains 12 records of 128 x 60 values

Nφ = 128
Nλ = 60 

function convert_lat_lon_bin(filepath, Nmonths)

    lat_lon_data = Array{Float32}(undef, Nλ*Nφ*Nmonths)

    read!(filepath, lat_lon_data)

    lat_lon_data = bswap.(lat_lon_data) |> Array{Float64}
    if Nmonths == 1 
        lat_lon_var  = reshape(lat_lon_data, Nφ, Nλ)
    else
        lat_lon_var  = reshape(lat_lon_data, Nφ, Nλ, Nmonths)
    end

    return lat_lon_var 
end


files = [("bathymetry_lat_lon_128x60",       "bathymetry",  1),
         ("sea_surface_temperature_25_128x60x12", "sst25", 12),
         ("tau_x_128x60x12",                      "tau_x", 12),
         ("tau_y_128x60x12",                      "tau_y", 12)]


for (filename, varname, Nmonths) in files

    var_data = convert_lat_lon_bin("$(filename)_FP32.bin", Nmonths)
    
    jldopen("$filename.jld2", "w") do file       
        if Nmonths == 1 
            file[varname] = var_data
        else
            for month in 1:Nmonths
                file["$varname/$month"] = var_data[:,:,month] 
            end
        end
    end
end
