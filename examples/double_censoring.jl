using Random, Distributions, .DiscretizeDistributions, Plots, IntervalArithmetic
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
histogram!(d_sims .* time_step .+ (time_step / 2), bins = time_steps, normalize = true,
    label = "Discrete approximation (adjusted for centering)", alpha = 0.5)
savefig("docs/src/assets/single_censoring.png")

## 2 Double censoring
continuous_waiting_time_two = truncated(Gamma(3.1, 2.1), lower = 1, upper = 100)

discrete_waiting_time_two = discretize(continuous_waiting_time_two, time_step, method = :left_aligned)
discrete_waiting_time_two.support ./= time_step

# simulate

c_sims_two = c_sims .+ rand(continuous_waiting_time_two, N)

d_sims_two = d_sims .+ rand(discrete_waiting_time_two, N) .+ 1

#summaries
d_sims_two_adjusted = (d_sims_two .* time_step)

println(round(mean(c_sims_two); digits = 2))
println(round(mean(d_sims_two_adjusted); digits = 2))

println(round(median(c_sims_two); digits = 2))
println(round(median(d_sims_two_adjusted); digits = 2))

println(round(var(c_sims_two); digits = 2))
println(round(var(d_sims_two_adjusted); digits = 2))

#visual
histogram(c_sims_two, bins = 0:40, label = "Continuous PDF",
    normalize = true, alpha = 0.5, xlabel = "Time to event C, from A")
histogram!(d_sims_two_adjusted, bins = 0:40, normalize = true,
    label = "Discrete approximation", alpha = 0.5)
savefig("docs/src/assets/double_censoring.png")
