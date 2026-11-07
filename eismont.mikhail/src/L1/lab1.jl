# для первого запуска (на последующие не влияет) - добавление Giac и Plots (Giac при импорте должен перекомпилироваться при первом запуске, для этого он должен быть добавлен в манифест-файл проекта Engee)
using Pkg
Pkg.add("Giac")
Pkg.add("ControlSystems")
Pkg.add("Plots")

# импорт библиотек, объявления
using Giac
using Giac.Commands: laplace, ilaplace, simplify, evalf, Dirac, Heaviside
using ControlSystems
using Plots

savefig_path_pre = "./lab1_plots/"

@giac_var t s K T T1 T1a T2

# параметры согласно варианту задания
K = 1.4
T = 0.5
T1 = 0.7
T2 = 0.9
# для апериодического звена 2-го порядка (T1 >= 2*T2)
T1a = 2.0

# передаточные функции
Wset = [
    K / s,
    K * s,
    K,
    K / (T * s + 1),
    K / (T2^2 * s^2 + T1a * s + 1),
    K / (T2^2 * s^2 + T1 * s + 1),
    K / (T * s^2 + 1),
    K / (s * (T * s + 1)),
    K * s / (T * s + 1),
    K * (T * s + 1),
    K * (T * s + 1) / s
]
# передаточные функции для численного решения (ControlSystems)
s_cs = tf("s")
Wset_cs = [
    K / s_cs,
    K * s_cs,
    K,
    K / (T * s_cs + 1),
    K / (T2^2 * s_cs^2 + T1a * s_cs + 1),
    K / (T2^2 * s_cs^2 + T1 * s_cs + 1),
    K / (T * s_cs^2 + 1),
    K / (s_cs * (T * s_cs + 1)),
    K * s_cs / (T * s_cs + 1),
    K * (T * s_cs + 1),
    K * (T * s_cs + 1) / s_cs
]

# входные и выходные сигналы
utStep = Heaviside(t)
utImpulse = Dirac(t)

UsStep = laplace(utStep(t), t, s)
UsImpulse = laplace(utImpulse(t), t, s)

YsStepSet = [UsStep * Ws for Ws in Wset]
YsImpulseSet = [UsImpulse * Ws for Ws in Wset]

ytStepSet = [ilaplace(Ys, s, t) for Ys in YsStepSet]
ytImpulseSet = [ilaplace(Ys, s, t) for Ys in YsImpulseSet]

# обход ограничения символьного вычисления ilaplace для случаев "s в числителе"
for i in range(1, size(Wset)[1])
    for ytSet in [ytStepSet, ytImpulseSet]
        local ytStr = string(ytSet[i])
        local should_be_dirac_derivative = occursin("ilaplace", ytStr)
        if should_be_dirac_derivative == true
            # выделение подстроки с выражением, переданным в ilaplace
            local pattern = r"ilaplace\(\[?([^)]*), *s, *t\]?\)"
            local matches = match(pattern, ytStr)
            if matches !== nothing
                # подготовка выражения для производной дельта-функции
                local innerStr = string(matches.captures[1])
                local innerExprNew = substitute(giac_eval(innerStr), s => 1)
                ytSet[i] = evalf(giac_eval(replace(ytStr, pattern => string("Dirac(", innerExprNew, ", 1)"))))
            end
        end
    end
end

# вычисление значений функций на интервале и построение графиков
t_range = 0.0:0.01:16.0
fytStepSet = [build_function(yt, t) for yt in ytStepSet]
fytImpulseSet = [build_function(yt, t) for yt in ytImpulseSet]

iFunction = [1, 1]
for ytSet in [fytStepSet, fytImpulseSet]
    global iFunction[2] = 1
    for iyt in range(1, size(ytSet)[1])
        yt = ytSet[iyt]
        local ySeq = map(
            y -> (isinf(y) == true || isnan(y) == true) ? NaN : Float64(y),
            yt.(t_range))
        # инициализация графика, показ символьного решения
        local p = plot(t_range, ySeq,
            xlabel = "t",
            ylabel = "y(t)",
            linewidth = 2,
            legend = false,
            color = :blue)
        # показ численного решения (кроме выражений W с s в числителе, не поддерживаемых step и impulse)
        if !(iFunction[2] in [2, 10])
            local W_cs = Wset_cs[iyt]
            local numeric_solve = nothing
            if iFunction[1] == 1
                numeric_solve = step(W_cs, t_range) # переходная характеристика
            elseif iFunction[1] == 2
                numeric_solve = impulse(W_cs, t_range) # импульсная характеристика
            end
            plot!(numeric_solve, line=:dash, linewidth = 2, legend = false, color = :violet)
        end
        # наложение асимптот, где были бесконечности (при вычислении помечены как NaN)
        for ty in collect(zip(t_range, ySeq))
            if isnan(ty[2])
                vline!([ty[1]],
                    line=:dot,
                    linewidth = 2,
                    legend = false,
                    color = :darkgreen)
            end
        end
        # принуждение к включению y=0 в отображаемую область графика с помощью маркерной метки в (0,0)
        scatter!([0], [0], markersize=0, primary=false)
        # сохранение
        savefig(p, string(savefig_path_pre, iFunction[1], "_", iFunction[2], ".svg"))
        global iFunction[2] += 1
    end
    global iFunction[1] += 1
end

# Графики переходного процесса апериодического звена при изменении T
T_alts = [T, T/2, T/5, T*2, T*5]

p = plot([], [])
for T_alt in T_alts
    W4cur = K / (T_alt * s + 1)
    YsStep = UsStep * W4cur
    ytStep = ilaplace(YsStep, s, t)
    plot!(t_range, build_function(ytStep, t).(t_range),
        xlabel = "t",
        ylabel = "y(t)",
        linewidth = 2,
        label = string("T = ", T_alt),
        legend = :bottomright)
end
savefig(p, string(savefig_path_pre, "0_step_4_variable_T.svg"))

# Графики переходного процесса звена №6 при изменении T1
T_alts = [2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
line_colors = ["#AAAA00", "#AA9900", "#AA8800", "#AA7700", "#AA6600", "#AA5500", "#AA4400", "#AA3300", "#AA2200", "#AA1100", "#AA0000"]

p = plot([], [])
for iT1_alt in range(1, size(T_alts)[1])
    T1_alt = T_alts[iT1_alt]
    W6cur = K / (T2^2 * s^2 + T1_alt * s + 1)
    YsStep = UsStep * W6cur
    ytStep = ilaplace(YsStep, s, t)
    plot!(t_range, build_function(ytStep, t).(t_range),
        xlabel = "t",
        ylabel = "y(t)",
        linewidth = 2,
        label = string("T1 = ", T1_alt),
        legend = :bottomright,
        color = line_colors[iT1_alt])
end
savefig(p, string(savefig_path_pre, "0_step_6_decreasing_T_to_type7.svg"))

# Графики переходного процесса в зависимости от коэффициента K на примере звена №6
K_alts = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
line_colors = ["#AAAA00", "#AA9900", "#AA8800", "#AA7700", "#AA6600", "#AA5500", "#AA4400", "#AA3300", "#AA2200", "#AA1100"]

p = plot([], [])
for iK_alt in range(1, size(K_alts)[1])
    K_alt = K_alts[iK_alt]
    W6cur = K_alt / (T2^2 * s^2 + T1 * s + 1)
    YsStep = UsStep * W6cur
    ytStep = ilaplace(YsStep, s, t)
    plot!(t_range, build_function(ytStep, t).(t_range),
        xlabel = "t",
        ylabel = "y(t)",
        linewidth = 2,
        label = string("K = ", K_alt),
        legend = :bottomright,
        color = line_colors[iK_alt])
end
savefig(p, string(savefig_path_pre, "0_step_6_decreasing_K.svg"))
