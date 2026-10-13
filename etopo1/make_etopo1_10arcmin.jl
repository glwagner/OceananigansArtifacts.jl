using NCDatasets

stride = 10 # arc-minutes between samples
url = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO1_Ice_g_gmt4.nc"

z, λ, φ = NCDataset(url) do dataset
    i = 1 + stride÷2 : stride : dataset.dim["lon"] - stride÷2
    j = 1 + stride÷2 : stride : dataset.dim["lat"] - stride÷2
    nomissing(dataset["z"][i, j]), dataset["lon"][i], dataset["lat"][j]
end

NCDataset("ETOPO1_Ice_g_10arcmin.nc", "c") do dataset
    dataset.attrib["title"] = "ETOPO1 ice-surface elevation sampled every 10 arc-minutes at cell centers"
    dataset.attrib["source"] = url
    defVar(dataset, "lon", λ, ("lon",), attrib = ["units" => "degrees_east"])
    defVar(dataset, "lat", φ, ("lat",), attrib = ["units" => "degrees_north"])
    defVar(dataset, "z", z, ("lon", "lat"), attrib = ["units" => "m"], deflatelevel = 9)
end
