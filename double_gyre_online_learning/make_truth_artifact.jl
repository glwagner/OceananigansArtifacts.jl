# Turn the raw truth output (truth_run_nh.jl) into the example's artifact:
#   truth_hourly_monthNN.jld2  hourly coarse u, v for month NN of the first recorded year, as FieldTimeSeries
#   truth_statistics.jld2      KE, MKE, EKE over all recorded years (daily samples), in yearly and 6-year blocks,
#                              and the time-mean coarse velocities
#   julia --project=twin twin/bench/make_truth_artifact.jl <truth_directory> <artifact_directory> [months]

using Oceananigans
using Oceananigans.Units
using JLD2
using Printf
using Statistics

truth_directory, artifact_directory = ARGS[1], ARGS[2]
months = parse(Int, get(ARGS, 3, "3"))
mkpath(artifact_directory)

const Lx = 4000kilometers
const Nx = 64
const spinup = 4 * 360days
const steps_per_hour = 15

grid = RectilinearGrid(Float32; size=(Nx, Nx), x=(0, Lx), y=(0, Lx), topology=(Bounded, Bounded, Flat))

# Snapshots written by a JLD2Writer: (time, u, v), sorted, keeping only times on the schedule
function read_snapshots(paths; on_schedule)
    snapshots = Dict{Float64, Tuple{Matrix{Float32}, Matrix{Float32}}}()
    for path in paths
        jldopen(path) do file
            for key in keys(file["timeseries/t"])
                key == "serialized" && continue
                t = Float64(file["timeseries/t/$key"])
                on_schedule(parse(Int, key), t) || continue
                snapshots[t] = (file["timeseries/u/$key"], file["timeseries/v/$key"])
            end
        end
    end
    times = sort(collect(keys(snapshots)))
    return times, [snapshots[t] for t in times]
end

##### Hourly training data, one file per month

parts = [joinpath(truth_directory, "coarse_velocities_part$p.jld2") for p in 1:12]
hourly_times, hourly = read_snapshots(parts; on_schedule = (iteration, t) -> iteration % steps_per_hour == 0 && t < spinup + 360days)
@info "Read $(length(hourly_times)) hourly snapshots spanning $(prettytime(hourly_times[end] - hourly_times[1]))"

u, v = XFaceField(grid), YFaceField(grid)
for month in 1:months
    hours = findall(t -> spinup + (month - 1) * 30days < t <= spinup + month * 30days, hourly_times)
    path = joinpath(artifact_directory, @sprintf("truth_hourly_month%02d.jld2", month))
    rm(path; force=true)
    times = hourly_times[hours] .- hourly_times[hours[1]]
    u_series = FieldTimeSeries{Face, Center, Center}(grid, times; backend=OnDisk(), path, name="u")
    v_series = FieldTimeSeries{Center, Face, Center}(grid, times; backend=OnDisk(), path, name="v")
    for (n, h) in enumerate(hours)
        set!(u, hourly[h][1])
        set!(v, hourly[h][2])
        set!(u_series, u, n)
        set!(v_series, v, n)
    end
    @info @sprintf("Wrote %s: %d hours, %.1f MB", path, length(hours), filesize(path) / 1e6)
end

##### Statistics from daily samples of all recorded years

daily_paths = filter(f -> startswith(basename(f), "coarse_velocities_daily"), readdir(truth_directory; join=true))
daily_times, daily = read_snapshots(vcat(parts, daily_paths); on_schedule = (iteration, t) -> iteration % (24steps_per_hour) == 0)
years = (daily_times[end] - daily_times[1]) / 360days
@info @sprintf("Statistics from %d daily snapshots, %.1f years", length(daily_times), years)

kinetic_energy(u, v) = (mean(abs2, u) + mean(abs2, v)) / 2

function block_statistics(snapshots)
    ū = mean(first.(snapshots))
    v̄ = mean(last.(snapshots))
    total = mean(kinetic_energy(s...) for s in snapshots)
    mean_flow = kinetic_energy(ū, v̄)
    return (; kinetic_energy = total, mean_kinetic_energy = mean_flow, eddy_kinetic_energy = total - mean_flow, ū, v̄)
end

all_years = block_statistics(daily)
blocks(length_days) = [block_statistics(daily[findall(t -> b <= (t - daily_times[1]) / length_days < b + 1, daily_times)])
                       for b in 0:floor(Int, years * 360days / length_days) - 1]
yearly = blocks(360days)
six_yearly = blocks(6 * 360days)

for (name, statistics) in (("yearly", yearly), ("6-year", six_yearly), ("all", [all_years]))
    @info "$name KE, MKE, EKE (1e-3 m² s⁻²): " * join((@sprintf("(%.2f, %.2f, %.2f)", 1e3 .* (s.kinetic_energy, s.mean_kinetic_energy, s.eddy_kinetic_energy)...) for s in statistics), " ")
end

path = joinpath(artifact_directory, "truth_statistics.jld2")
jldsave(path; kinetic_energy = all_years.kinetic_energy, mean_kinetic_energy = all_years.mean_kinetic_energy,
        eddy_kinetic_energy = all_years.eddy_kinetic_energy, mean_u = all_years.ū, mean_v = all_years.v̄,
        yearly = [(s.kinetic_energy, s.mean_kinetic_energy, s.eddy_kinetic_energy) for s in yearly],
        years, sampling = "daily")
@info "Wrote $path"
