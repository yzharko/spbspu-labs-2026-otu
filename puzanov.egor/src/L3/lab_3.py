import control.matlab as cm
import matplotlib.pyplot as plt

s = cm.tf("s")

LINKS = {
     1: ("Интегрирующее (идеальное)",     lambda p: p["K"] / s),
     2: ("Дифференцирующее (идеальное)",  lambda p: p["K"] * s),
     3: ("Усилительное (безынерционное)", lambda p: cm.tf([float(p["K"])], [1.0])),
     4: ("Апериодическое 1-го порядка",   lambda p: p["K"] / (p["T"] * s + 1)),
     5: ("Апериодическое 2-го порядка",   lambda p: p["K"] / (p["T2"]**2 * s**2 + p["T1"] * s + 1)),
     6: ("Колебательное",                 lambda p: p["K"] / (p["T2"]**2 * s**2 + p["T1"] * s + 1)),
     7: ("Консервативное",                lambda p: p["K"] / (p["T"]**2 * s**2 + 1)),
     8: ("Инерционное интегрирующее",     lambda p: p["K"] / (s * (p["T"] * s + 1))),
     9: ("Реальное дифференцирующее",     lambda p: p["K"] * s / (p["T"] * s + 1)),
    10: ("Форсирующее 1-го порядка",      lambda p: p["K"] * (p["T"] * s + 1)),
    11: ("Изодромное",                    lambda p: p["K"] * (p["T"] * s + 1) / s),
}

PARAMS = {i: {"K": 2.2, "T": 0.2} for i in LINKS}
PARAMS[5].update({"T1": 0.7, "T2": 0.3})
PARAMS[6].update({"T1": 0.4, "T2": 0.3})

# ЛАЧХ и ЛФЧХ
fig, axes = plt.subplots(6, 4, figsize=(24, 32))
for i in LINKS:
    name = LINKS[i][0]
    W = LINKS[i][1](PARAMS[i])
    # чтобы два графика одного звена были друг под другом
    r, c = divmod(i-1, 4)
    mag_ax, phase_ax = axes[2*r, c], axes[2*r+1, c]

    cm.bode(W, dB=True, ax=[mag_ax, phase_ax], grid=True, color=f"C{i}")
    mag_ax.set_title(f"№ {i}, {name}: ЛАЧХ", fontsize=9)
    mag_ax.set_ylabel("аплитуда, дБ", fontsize=9)
    mag_ax.set_xlabel("частота, рад/с", fontsize=9)
    phase_ax.set_title(f"№ {i}, {name}: ЛФЧХ", fontsize=9)
    phase_ax.set_ylabel("фаза, град", fontsize=9)
    phase_ax.set_xlabel("частота, рад/с", fontsize=9)
    
for ax in axes.flat:
    if not ax.lines:
        ax.axis("off")

fig.savefig("bode.pdf", bbox_inches="tight")
fig.subplots_adjust(left=0.036, bottom=0.046, right=0.994, top=0.944, wspace=0.162, hspace=0.81)
plt.show()

# АФЧХ
fig, axes = plt.subplots(4, 3, figsize=(24, 28))
for i in LINKS:
    ax = axes.flat[i-1]
    name = LINKS[i][0]
    W = LINKS[i][1](PARAMS[i])
    
    cm.nyquist(W, ax=ax, color=f"C{i}")
    ax.grid(True)
    ax.set_title(f"№ {i}, {name}: АФЧХ", fontsize=9)
    ax.set_xlabel("Реальная часть", fontsize=9)
    ax.set_ylabel("Мнимая часть", fontsize=9)

for ax in axes.flat:
    if not ax.lines:
        ax.axis("off")
    
fig.savefig("nyquist.pdf", bbox_inches="tight")
fig.subplots_adjust(left=0.032, bottom=0.035, right=0.994, top=0.968, wspace=0.12, hspace=0.45)
plt.show()

# Изучение влияния параметров
K0, T0 = 1.0, 1.0
T_vals = [T0 *0.2, T0 *0.5, T0, T0*1.5, T0*2]
K_vals = [K0 *0.2, K0 *0.5, K0, K0*1.5, K0*2]

# отрисовывает графики с влиянием параметров
def varying_bode(
        data_to_plot: list[tuple[int, dict[str, float], str]], # список троек "(номер звена, параметры, какой параметр менять на графике)"
        varying_alues: list[float], # список значений что будет принимать меняющийся параметр
        save_image_name: str
        ):
    fig, axes = plt.subplots(6, 4, figsize=(24, 32))
    for idx, (link_index, param, key) in enumerate(data_to_plot):
        name = LINKS[link_index][0]
        r, c = divmod(idx, 4)
        mag_ax, phase_ax = axes[2*r, c], axes[2*r+1, c]

        for v in varying_alues:
            param[key] = v
            W = LINKS[link_index][1](param)
            cm.bode(W, dB=True, ax=[mag_ax, phase_ax], grid=True, label=f"{key} = {v:g}")

        mag_ax.set_title(f"№ {link_index}, {name}: ЛАЧХ", fontsize=9)
        mag_ax.set_ylabel("аплитуда, дБ", fontsize=9)
        mag_ax.set_xlabel("частота, рад/с", fontsize=9)
        mag_ax.legend()
        phase_ax.set_title(f"№ {link_index}, {name}: ЛФЧХ", fontsize=9)
        phase_ax.set_ylabel("фаза, град", fontsize=9)
        phase_ax.set_xlabel("частота, рад/с", fontsize=9)
        phase_ax.legend()

    for ax in axes.flat:
        if not ax.lines:
            ax.axis("off")

    fig.savefig(save_image_name, bbox_inches="tight")
    fig.subplots_adjust(left=0.036, bottom=0.046, right=0.994, top=0.944, wspace=0.162, hspace=0.58)
    plt.show()

# Влияние Т
data_to_plot = []
for i in range(4, 12):
    if i not in (5, 6):
        data_to_plot.append((i, PARAMS[i].copy(), "T"))
    elif i == 5:
        p = PARAMS[i].copy()
        p["T2"] = T_vals[0] / 2
        data_to_plot.append((i, p.copy(), "T1"))
        p["T1"] = T_vals[-1] * 2
        data_to_plot.append((i, p.copy(), "T2"))
    else:
        p = PARAMS[i].copy()
        p["T2"] = T_vals[-1] / 2
        data_to_plot.append((i, p.copy(), "T1"))
        p["T1"] = T_vals[0] * 2
        data_to_plot.append((i, p.copy(), "T2"))

varying_bode(data_to_plot, T_vals, "variyng_T_bode.pdf")

# влияние K
data_to_plot = [(i, PARAMS[i].copy(), "K") for i in range(1, 12)]
varying_bode(data_to_plot, K_vals, "variyng_K_bode.pdf")

