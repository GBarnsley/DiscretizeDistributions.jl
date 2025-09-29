using Random, Distributions, .DiscretizeDistributions, Plots, IntervalArithmetic, StatsPlots
Random.seed!(1234)
N = 100000

## 1 Single censoring
continuous_waiting_time = truncated(Gamma(3.1, 2.1), lower = 1, upper = 20)
time_step = 1//1
discrete_waiting_time = discretize(continuous_waiting_time, time_step, method = :left_aligned)
discrete_waiting_time.support ./= time_step

# simulate
c_sims = rand(continuous_waiting_time, N)
d_sims = rand(discrete_waiting_time, N)

# summaries
centred_discrete = discretize(continuous_waiting_time, time_step, method = :centred)
println(round(mean(c_sims); digits = 2))
println(round(mean(centred_discrete); digits = 2))

println(round(median(c_sims); digits = 2))
println(round(median(centred_discrete); digits = 2))

println(round(var(c_sims); digits = 2))
println(round(var(centred_discrete); digits = 2))

# visual comparison
time_steps = range(minimum(continuous_waiting_time), maximum(continuous_waiting_time))
plot(time_steps, pdf.(continuous_waiting_time, time_steps),
    label = "Continuous PDF", lw = 2, xlabel = "Time to event B, from A")
values = discrete_waiting_time.support
counts = [sum(v .== d_sims) for v in values] ./ N
bar!(values .* time_step .+ (time_step / 2), counts,
    label = "Discrete approximation (bar height adjusted for centering)", alpha = 0.5, bar_width = 1.0)
savefig("docs/src/assets/single_censoring.png")

## 2 Double censoring
continuous_waiting_time_two = truncated(Gamma(3.1, 2.1), lower = 1, upper = 100)

discrete_waiting_time_two = discretize(continuous_waiting_time_two, time_step, method = :left_aligned)
discrete_waiting_time_two.support ./= time_step

# simulate

c_sims_two = c_sims .+ rand(continuous_waiting_time_two, N)

d_sims_two = d_sims .+ rand(discrete_waiting_time_two, N)

#summaries
d_sims_two_adjusted = ((d_sims_two .+ 1) .* time_step)

println(round(mean(c_sims_two); digits = 2))
println(round(mean(d_sims_two_adjusted); digits = 2))

println(round(median(c_sims_two); digits = 2))
println(round(median(d_sims_two_adjusted); digits = 2))

println(round(var(c_sims_two); digits = 2))
println(round(var(d_sims_two_adjusted); digits = 2))

#visual
density(c_sims_two, label = "Continuous PDF", xlabel = "Time to event C, from A")
values = [discrete_waiting_time.support..., discrete_waiting_time_two.support .+ maximum(discrete_waiting_time.support)...]
values = values[values .< 40]
counts = [sum(v .== d_sims_two) for v in values] ./ N
bar!(values .* time_step .+ (time_step / 2), counts,
    label = "Discrete approximation\n(bar height adjusted for centering)", alpha = 0.5, bar_width = 1.0)
savefig("docs/src/assets/double_censoring.png")
