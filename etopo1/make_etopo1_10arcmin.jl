using NCDatasets
using Downloads

stride = 10 # arc-minutes between samples
offset = stride ÷ 2

thredds_url = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO1_Ice_g_gmt4.nc"
erddap_url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude" *
             "[$offset:$stride:$(10800 - offset)][$offset:$stride:$(21600 - offset)]"

read_thredds() = NCDataset(thredds_url) do dataset
    i = 1 + offset : stride : dataset.dim["lon"] - offset
    j = 1 + offset : stride : dataset.dim["lat"] - offset
    nomissing(dataset["z"][i, j]), dataset["lon"][i], dataset["lat"][j], thredds_url
end

read_erddap() = NCDataset(Downloads.download(erddap_url)) do dataset
    Float32.(nomissing(dataset["altitude"][:, :])), dataset["longitude"][:], dataset["latitude"][:], erddap_url
end

z, λ, φ, url = try
    read_thredds()
catch error
    @warn "Reading ETOPO1 from NCEI failed; trying CoastWatch ERDDAP" error
    read_erddap()
end

NCDataset("ETOPO1_Ice_g_10arcmin.nc", "c") do dataset
    dataset.attrib["title"] = "ETOPO1 ice-surface elevation sampled every 10 arc-minutes at cell centers"
    dataset.attrib["source"] = url
    defVar(dataset, "lon", λ, ("lon",), attrib = ["units" => "degrees_east"])
    defVar(dataset, "lat", φ, ("lat",), attrib = ["units" => "degrees_north"])
    defVar(dataset, "z", z, ("lon", "lat"), attrib = ["units" => "m"], deflatelevel = 9)
end
