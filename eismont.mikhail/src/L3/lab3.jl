using ControlSystems
using Plots

savefig_path_pre = "./lab3_plots/"

s = tf('s')
Wset = [
    K -> K / s,
    K -> K * s,
    K -> K + 0 * s,
    (K, T) -> K / (T * s + 1),
    (K, T1a, T2) -> K / (T2^2 * s^2 + T1a * s + 1),
    (K, T1, T2) -> K / (T2^2 * s^2 + T1 * s + 1),
    (K, T) -> K / (T * s^2 + 1),
    (K, T) -> K / (s * (T * s + 1)),
    (K, T) -> K * s / (T * s + 1),
    (K, T) -> K * (T * s + 1),
    (K, T) -> K * (T * s + 1) / s
]

# наборы коэффициентов: первый элемент - согласно варианту, остальные - альтернативные
Kset   = [1.4,  0.5, 1.0, 2.0, 3.0, 5.0]
Tset   = [0.5,  0.4, 0.7, 1.0, 2.0, 5.0]
T1set  = [0.7,  0.4, 0.5, 0.8, 1.0, 1.3]
T1aset = [2.5,  2.4, 2.7, 3.5, 4.2, 5.0]
T2set  = [0.9,  0.7, 0.8, 1.0, 1.1, 1.2]

Kcolors = ["#00FF00", "#770000", "#993300", "#BB6600", "#DD9900", "#FFCC00"]
T1colors = ["#00FF00", "#550088", "#773399", "#9966AA", "#BB99BB", "#DDCCCC"]
T2colors = ["#009966", "#008855", "#007744", "#006633", "#005522"]

for i in 1:1:length(Wset)
    local W = nothing
    local paramSetsK = nothing
    local paramSetsT = []
    if i in [1, 2, 3]
        paramSetsK = collect(zip(Kcolors, Kset))
    elseif i in [5, 6]
        local T12pairs = collect(zip(
            [
                i == 5 ? T1aset : T1set;
                fill(i == 5 ? T1aset[1] : T1set[1], length(T2set) - 1)
            ], [
                fill(T2set[1], length(i == 5 ? T1aset : T1set));
                T2set[2:length(T2set)]
            ]
        ))
        paramSetsK = collect(zip(
            Kcolors,
            Kset,
            fill((i == 5 ? T1aset[1] : T1set[1], T2set[1]), length(Kset))
        ))
        paramSetsT = collect(zip(
            [T1colors; T2colors],
            fill(Kset[1], length(T12pairs)),
            T12pairs
        ))
    elseif i in [4, 7, 8, 9, 10, 11]
        paramSetsK = collect(zip(
            Kcolors,
            Kset,
            fill(Tset[1], length(Kset))
        ))
        paramSetsT = collect(zip(
            T1colors[1:length(Tset)],
            fill(Kset[1], length(Tset)),
            Tset
        ))
    end
    if paramSetsK != nothing
        for (plotFunc, plotOverFunc, fnPrefix) in [(nyquistplot, nyquistplot!, "nqst"), (bodeplot, bodeplot!, "bode")]
            for (paramSets, varyingVar) in [(paramSetsK, "K"), (paramSetsT, "T")]
                if length(paramSets) > 0
                    local p = nothing
                    local plotCreated = false
                    for j in length(paramSets):-1:1
                        local paramSet = paramSets[j]
                        local (color, K, T1, T2) = fill(nothing, 4)
                        local W = nothing
                        if i in [1, 2, 3]
                            (color, K) = paramSet
                            W = Wset[i](K)
                        elseif i in [5, 6]
                            (color, K, (T1, T2)) = paramSet
                            W = Wset[i](K, T1, T2)
                        elseif i in [4, 7, 8, 9, 10, 11]
                            (color, K, T1) = paramSet
                            W = Wset[i](K, T1)
                        end
                        if plotCreated
                            plotOverFunc(W,
                                label=string("S", j),
                                title=string("Chain#", i, ", <", varyingVar, ">"),
                                color=color,
                                line=j == 1 ? :solid : :dot,
                                linewidth=(j == 1 ? 2 : 3),
                                legend=:outertopright,
                                xlims=:auto, ylims=:auto
                            )
                        else
                            p = plotFunc(W,
                                label=string("S", j),
                                title=string("Chain#", i, ", <", varyingVar, ">"),
                                color=color,
                                line=j == 1 ? :solid : :dot,
                                linewidth=(j == 1 ? 2 : 3),
                                legend=:outertopright,
                                xlims=:auto, ylims=:auto
                            )
                            plotCreated = true
                        end
                    end
                    if p != nothing
                        display(p)
                        savefig(p, string(savefig_path_pre, fnPrefix, i, "_", varyingVar, ".svg"))
                    end
                end
            end
        end
    end
end
