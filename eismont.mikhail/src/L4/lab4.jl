using ControlSystems
using Plots

plot_path_pre = "./lab4_plots/"

coefs_nonclosed = [[10], [1, 3, 2, 2, 1, 1]]
coefs_closed = [[10], [1, 3, 2, 2, 1, 11]]
W_nonclosed = tf(coefs_nonclosed[1], coefs_nonclosed[2])
W_closed = tf(coefs_closed[1], coefs_closed[2])

# Годограф Найквиста
p = nyquistplot(W_nonclosed, xlabel="Re A(w)", ylabel="Im A(w)")
display(p)
savefig(p, string(plot_path_pre, "nyquist.svg"))

# Годограф Михайлова
D_closed = tf(coefs_closed[2], coefs_closed[1])

omega_range = 0:0.0001:1.25
D_closed_values = [D_closed(im * omega) for omega in omega_range]
U_closed = [v[1] for v in real.(D_closed_values)]
V_closed = [v[1] for v in imag.(D_closed_values)]

p = plot(U_closed, V_closed,
    xlabel="Re D(w)", ylabel="Im D(w)")
scatter!([U_closed[1]], [V_closed[1]], label="w = 0", marker=:circle)

display(p)
savefig(p, string(plot_path_pre, "mikhaylov.svg"))

# Влияние переходных процессов
for (W_cur, fnsuffix, secs) in ((W_nonclosed, "nonclosed", 30), (W_closed, "closed", 21.6))
    step_timerange = range(0, secs, length=1000)
    local step_W = step(W_cur, step_timerange)
    global p = plot(step_W)
    display(p)
    savefig(p, string(plot_path_pre, "step_", fnsuffix, ".svg"))
end
