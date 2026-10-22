# Truth for the online-learning double gyre example: a rigid-lid 2D NonhydrostaticModel (the coarse model's equations)
# at 1024² (3.9 km) with WENO(order=9) and no closure. After a 4-year spin-up, saves coarse velocities on the
# 64² C-grid, Gaussian-filtered with σ = 62.5 km as in Maddison (2026): hourly for one year (training data),
# and daily for twelve years (statistics).
#
# The filter acts on the streamfunction ψ = -∫₀ʸ u dy, which vanishes on the walls: ψ is filtered with odd reflection
# at the walls and sampled at the coarse cell corners, and the coarse face velocities are differences of ψ. They are
# then exactly divergence-free with no flow through the walls. Without the filter, this is the box average of the
# face velocities.
#
#   julia --project=twin twin/bench/truth_run_nh.jl <output_directory>
#   TRUTH_SELF_TEST=true julia --project=twin twin/bench/truth_run_nh.jl

using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: on_architecture
using CUDA
using JLD2
using Printf

const L = 2000kilometers
const Nf = 1024
const Nc = 64
const coarsening = Nf ÷ Nc
const Δc = 2L / Nc
const σ = 62.5kilometers / (2L / Nf)  # filter width in fine cells
const spinup = 4 * 360days
const hourly_recording = 360days      # training data
const recording = 12 * 360days        # statistics, as Maddison's 12-year averages

parameters = (τ₀=0.1, ρ₀=1000.0, H=500.0, r=1e-7, L)
@inline wind_and_drag(x, y, t, u, p) = - p.τ₀ * cos(π * y / p.L) / (p.ρ₀ * p.H) - p.r * u
@inline drag(x, y, t, v, p) = - p.r * v

function truth_model()
    grid = RectilinearGrid(GPU(); size=(Nf, Nf), x=(0, 2L), y=(0, 2L), halo=(6, 6),
                           topology=(Bounded, Bounded, Flat))
    forcing = (u = Forcing(wind_and_drag; field_dependencies=:u, parameters),
               v = Forcing(drag; field_dependencies=:v, parameters))
    return NonhydrostaticModel(grid; advection=WENO(order=9), coriolis=BetaPlane(f₀=7e-5, β=2e-11),
                               forcing, timestepper=:RungeKutta3)
end

# Gaussian weights (truncated at 4σ) for sampling at the coarse corners, with the odd reflection about the wall
# corners folded in, so that filtering and coarse-graining is ψc = G ψ Gᵀ with G of size (Nc+1) × (Nf+1).
function coarse_graining_matrix(; filter=true)
    n = Nf + 1
    G = zeros(Nc + 1, n)
    radius = filter ? ceil(Int, 4σ) : 0
    weights = [exp(-k^2 / 2σ^2) for k in -radius:radius]
    weights ./= sum(weights)
    for I in 1:Nc+1, (w, k) in zip(weights, -radius:radius)
        j = 1 + coarsening * (I - 1) + k
        j, sign = j < 1 ? (2 - j, -1) : j > n ? (2n - j, -1) : (j, 1)
        G[I, j] += sign * w
    end
    return on_architecture(GPU(), G)
end

struct CoarseVelocities{P, G}
    ψ :: P
    G :: G
end

CoarseVelocities(model; filter=true) =
    CoarseVelocities(Field(CumulativeIntegral(-model.velocities.u, dims=2)), coarse_graining_matrix(; filter))

function (coarse::CoarseVelocities)(model)
    compute!(coarse.ψ)
    ψc = on_architecture(CPU(), coarse.G * interior(coarse.ψ, :, :, 1) * coarse.G')
    return (u = Float32.(-diff(ψc, dims=2) ./ Δc), v = Float32.(diff(ψc, dims=1) ./ Δc))
end

function progress(simulation)
    u, v = simulation.model.velocities
    kinetic_energy = (sum(u^2) / length(u) + sum(v^2) / length(v)) / 2  # whole-domain mean
    @printf("TRUTH day %5.0f  KE %.4e  max|u| %.3f  wall %s\n", simulation.model.clock.time / day, kinetic_energy,
            maximum(abs, u), prettytime(simulation.run_wall_time))
    flush(stdout)
    return nothing
end

function main(directory)
    mkpath(directory)
    model = truth_model()
    set!(model, u=(x, y) -> 1e-3 * sin(π * x / 2L) * sin(2π * y / 2L))
    simulation = Simulation(model; Δt=4minutes, stop_time=spinup)
    add_callback!(simulation, progress, TimeInterval(30days))
    simulation.output_writers[:checkpointer] = Checkpointer(model; dir=directory, prefix="checkpoint",
                                                            schedule=TimeInterval(30days), cleanup=true)
    run!(simulation, pickup=true)  # the spin-up; picking up from a checkpoint past it returns immediately

    u, v, w = model.velocities
    coarse = CoarseVelocities(model)
    coarse_velocities = (u = model -> coarse(model).u, v = model -> coarse(model).v)
    simulation.output_writers[:fine] = JLD2Writer(model, (; u, v); dir=directory, filename="fine_velocities",
                                                  schedule=TimeInterval(30days), overwrite_files=false)

    # Hourly coarse velocities for the first recorded year (training data), then daily ones (statistics)
    if model.clock.time < spinup + hourly_recording
        simulation.output_writers[:hourly] = JLD2Writer(model, coarse_velocities; dir=directory, filename="coarse_velocities",
                                                        schedule=TimeInterval(1hour), file_splitting=TimeInterval(30days),
                                                        overwrite_files=false)
        simulation.stop_time = spinup + hourly_recording
        run!(simulation)
        pop!(simulation.output_writers, :hourly)
    end

    simulation.output_writers[:daily] = JLD2Writer(model, coarse_velocities; dir=directory, filename="coarse_velocities_daily",
                                                   schedule=TimeInterval(1day), file_splitting=TimeInterval(360days),
                                                   overwrite_files=false)
    simulation.stop_time = spinup + recording
    run!(simulation)
    return nothing
end

function self_test()
    model = truth_model()
    set!(model, u=(x, y) -> 0.3 * sin(π * x / 2L) * sin(2π * y / 2L) + 0.1 * sin(9π * x / 2L) * sin(13π * y / 2L))
    for _ in 1:20
        time_step!(model, 4minutes)
    end
    box, gaussian = CoarseVelocities(model; filter=false)(model), CoarseVelocities(model)(model)

    u, v = (on_architecture(CPU(), interior(f, :, :, 1)) for f in model.velocities[(:u, :v)])
    box_u = dropdims(sum(reshape(u[1:coarsening:end, :], Nc + 1, coarsening, Nc), dims=2), dims=2) ./ coarsening
    divergence(c) = maximum(abs, diff(c.u, dims=1) .+ diff(c.v, dims=2)) / maximum(abs, c.u)
    rms_vorticity(c) = sqrt(sum(abs2, diff(c.v, dims=1)[:, 2:end-1] .- diff(c.u, dims=2)[2:end-1, :]) / (Nc - 1)^2) / Δc
    coarse = CoarseVelocities(model)
    CUDA.@sync coarse(model)
    elapsed = @elapsed for _ in 1:10; CUDA.@sync coarse(model); end
    @printf("SELFTEST unfiltered u vs box average: %.2e relative\n", maximum(abs, box.u .- box_u) / maximum(abs, box_u))
    @printf("SELFTEST relative divergence: box %.2e, Gaussian %.2e; max |normal flow| on walls %.2e\n",
            divergence(box), divergence(gaussian), max(maximum(abs, gaussian.u[[1, end], :]), maximum(abs, gaussian.v[:, [1, end]])))
    @printf("SELFTEST rms vorticity ratio Gaussian/box %.2f; %.1f ms per coarse-graining\n",
            rms_vorticity(gaussian) / rms_vorticity(box), 100elapsed)
end

get(ENV, "TRUTH_SELF_TEST", "false") == "true" ? self_test() : main(ARGS[1])
