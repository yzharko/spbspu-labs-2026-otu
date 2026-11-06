import control.matlab as cm
import sympy as sp

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = (10, 5)
plt.rcParams["font.size"] = 11
plt.rcParams["axes.grid"] = True

s, t = sp.symbols('s t')
K, T, T1, T2 = sp.symbols('K T T1 T2', positive=True)

LINKS = {
    1:  ("Интегрирующее (идеальное)",     K/s),
    2:  ("Дифференцирующее (идеальное)",  K*s),
    3:  ("Усилительное (безынерционное)", K),
    4:  ("Апериодическое 1-го порядка",   K/(T*s + 1)),
    5:  ("Апериодическое 2-го порядка",   K/(T2**2*s**2 + T1*s + 1)),
    6:  ("Колебательное",                 K/(T2**2*s**2 + T1*s + 1)),
    7:  ("Консервативное",                K/(T**2*s**2 + 1)),
    8:  ("Инерционное интегрирующее",     K/(s*(T*s + 1))),
    9:  ("Реальное дифференцирующее",     K*s/(T*s + 1)),
    10: ("Форсирующее 1-го порядка",      K*(T*s + 1)),
    11: ("Изодромное",                    K*(T*s + 1)/s),
}

PARAMS = {i: {K: 2.2, T: 0.2} for i in LINKS}
PARAMS[5].update({T1: 0.7, T2: 0.3})   # T1 >= 2*T2  => T1 отличается от варианта
PARAMS[6].update({T1: 0.4, T2: 0.3})   

results_sympy = {}
time          = np.linspace(0, 10, 1001)

U_step    = sp.laplace_transform(sp.Heaviside(t),  t, s, noconds=True)
U_impulse = sp.laplace_transform(sp.DiracDelta(t), t, s, noconds=True)

for i in LINKS:
    name, W = LINKS[i]
    W_num = sp.simplify(W.subs(PARAMS[i]))

    Y_step    = sp.simplify(W_num * U_step)
    h         = sp.simplify(sp.inverse_laplace_transform(Y_step, s, t))

    Y_impulse = sp.simplify(W_num * U_impulse)
    w         = sp.simplify(sp.inverse_laplace_transform(Y_impulse, s, t))
    print(f"№{i}, h: {h}")
    print(f"№{i}, w: {w}")
    print()
    results_sympy[i] = (name, h, w)

# Функция для перевода из объекта sympy в массив числел для графика
def to_numpy(expr, tt):
    e = sp.simplify(expr.subs(sp.Heaviside(t), sp.Integer(1)))
    has_delta = bool(e.has(sp.DiracDelta))
    if has_delta:
        e = e.replace(sp.DiracDelta, lambda *a: sp.Integer(0))
    f = sp.lambdify(t, e, modules='numpy')
    y = np.asarray(f(tt))

    return np.broadcast_to(y, tt.shape).astype(float), has_delta

# Графики h(t) (sympy)
fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat

for i in LINKS:
    ax = axes[i - 1]
    name, h_sym, _ = results_sympy[i]
    ax.set_title(f"№{i}. {name}\nh(t)", fontsize=9)
    ax.set_xlabel('t, c'); ax.set_ylabel('h(t)')
    y, has_delta = to_numpy(h_sym, time)
    ax.plot(time, y, lw=1.8, color=('r' if has_delta else 'b'))

axes[-1].axis('off')
fig.suptitle('Переходные функции h(t)', fontsize=13)
plt.tight_layout()
plt.show()

# Графики w(t) (sympy)
fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat

for i in LINKS:
    ax = axes[i - 1]
    name, _, w_sym = results_sympy[i]
    ax.set_title(f"№{i}. {name}\nw(t)", fontsize=9)
    ax.set_xlabel('t, c'); ax.set_ylabel('w(t)')
    y, has_delta = to_numpy(w_sym, time)
    ax.plot(time, y, lw=1.8, color=('r' if has_delta else 'g'))

axes[-1].axis('off')
fig.suptitle('Импульсные переходные функции w(t)', fontsize=13)
plt.tight_layout()
plt.show()

# Проверка при помощи пакета contol.matlab

results_ct = {}
s_tf = cm.tf("s")

LINKS_CT = {
    1: lambda p: p[K] / s_tf,
    2: lambda p: p[K] * s_tf,
    3: lambda p: cm.tf([float(p[K])], [1.0]),
    4: lambda p: p[K] / (p[T] * s_tf + 1),
    5: lambda p: p[K] / (p[T2] ** 2 * s_tf**2 + p[T1] * s_tf + 1),
    6: lambda p: p[K] / (p[T2] ** 2 * s_tf**2 + p[T1] * s_tf + 1),
    7: lambda p: p[K] / (p[T] ** 2 * s_tf**2 + 1),
    8: lambda p: p[K] / (s_tf * (p[T] * s_tf + 1)),
    9: lambda p: p[K] * s_tf / (p[T] * s_tf + 1),
    10: lambda  p: p[K] * (p[T] * s_tf + 1),
    11: lambda  p: p[K] * (p[T] * s_tf + 1) / s_tf,
}

for i in LINKS:
    name = LINKS[i][0]
    W_check = LINKS_CT[i](PARAMS[i])

    try:
        y_h, _ = cm.step(W_check,    T=time)
        y_w, _ = cm.impulse(W_check, T=time)
        y_h    = np.squeeze(y_h)
        y_w    = np.squeeze(y_w)
        is_improper = False
        
    except Exception as e:
        print(e)
        y_h = None
        y_w = None
        is_improper = True

    results_ct[i] = (name, y_h, y_w, is_improper)

# Графики h(t) по control.matlab
fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat

for i in LINKS:
    ax = axes[i - 1]
    name, y_h, _, is_imp = results_ct[i]
    ax.set_title(f"№{i}. {name}\nh(t) — control", fontsize=9)
    ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
    if is_imp:
        ax.text(0.1, 0.5, "не удалось построить график для проверки")
    else:
        ax.plot(time, y_h, color="b", lw=1.8)

axes[-1].axis('off')
fig.suptitle('Переходные функции h(t) (control.matlab)', fontsize=13)
plt.tight_layout()
plt.show()

# Графики w(t) по control.matlab
fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat

for i in LINKS:
    ax = axes[i - 1]
    name, _, y_w, is_imp = results_ct[i]
    ax.set_title(f"№{i}. {name}\nw(t) — control", fontsize=9)
    ax.set_xlabel('t, c'); ax.set_ylabel('w(t)'); ax.grid(True)
    if is_imp:
        ax.text(0.1, 0.5, "не удалось построить график для проверки")
    else:
        ax.plot(time, y_w, color="g", lw=1.8)

axes[-1].axis('off')
fig.suptitle('Импульсные переходные функции w(t) — (control.matlab)', fontsize=13)
plt.tight_layout()
plt.show()

# Исследование влияния T
K0, T0 = 1.0, 1.0
T_vals = [T0 *0.2, T0 *0.5, T0, T0*1.5, T0*2]
tmax   = 10
time_short = np.linspace(0, tmax, 4000)
time_long = np.linspace(0, tmax*10, 4000*10)

U_step = sp.laplace_transform(sp.Heaviside(t), t, s, noconds=True)

fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat
link_index = 4
graph_index = 0
#  кол-во графиков равно количеству звеньев в которых включен параметр T (все 11 - 3 первых) + 2 тк у звеньев 5 и 6 нужны дополнительные графики для второго параметра
number_of_graphs = len(LINKS) - (link_index-1) + 2
while graph_index < number_of_graphs:

    if link_index != 5 and link_index != 6:
        ax = axes[graph_index]
        name, W    = LINKS[link_index]
        for Tv in T_vals:
            W_num  = sp.simplify(W.subs({K: K0, T: Tv}))
            Y_step = sp.simplify(W_num * U_step)
            h_T    = sp.inverse_laplace_transform(Y_step, s, t)
            y, _   = to_numpy(h_T, time_short)
            ax.plot(time_short, y, label=f'T = {Tv:g}')
        ax.set_title(f"№{link_index}. {name}\nh(t)", fontsize=9)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.legend(loc=("upper right" if link_index == 7 else "best"))
        graph_index += 1
        link_index += 1

    elif link_index == 5:
        ax = axes[graph_index]
        name, W    = LINKS[link_index]
        for Tv in T_vals:
            W_num  = sp.simplify(W.subs({K: K0, T1: Tv, T2: T_vals[0]/2}))
            Y_step = sp.simplify(W_num * U_step)
            h_T    = sp.inverse_laplace_transform(Y_step, s, t)
            y, _   = to_numpy(h_T, time_short)
            ax.plot(time_short, y, label=f'T1 = {Tv:g}')
        ax.set_title(f"№{link_index}. {name}\nh(t)", fontsize=9)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.legend()
        graph_index += 1

        ax = axes[graph_index]
        name, W    = LINKS[link_index]
        for Tv in T_vals:
            W_num  = sp.simplify(W.subs({K: K0, T1: T_vals[-1]*2, T2: Tv}))
            Y_step = sp.simplify(W_num * U_step)
            h_T    = sp.inverse_laplace_transform(Y_step, s, t)
            y, _   = to_numpy(h_T, time_short)
            ax.plot(time_short, y, label=f'T2 = {Tv:g}')
        ax.set_title(f"№{link_index}. {name}\nh(t)", fontsize=9)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.legend()
        graph_index += 1
        link_index += 1

    else:
        ax = axes[graph_index]
        name, W    = LINKS[link_index]
        for Tv in T_vals:
            W_num  = sp.simplify(W.subs({K: K0, T1: Tv, T2: T_vals[-1]*2}))
            Y_step = sp.simplify(W_num * U_step)
            h_T    = sp.inverse_laplace_transform(Y_step, s, t)
            y, _   = to_numpy(h_T, time_long)
            ax.plot(time_long, y, label=f'T1 = {Tv:g}')
        ax.set_title(f"№{link_index}. {name}\nh(t)", fontsize=9)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.legend(loc='upper right')
        graph_index += 1

        ax = axes[graph_index]
        name, W    = LINKS[link_index]
        for Tv in T_vals:
            W_num  = sp.simplify(W.subs({K: K0, T1: T_vals[0]/2, T2: Tv}))
            Y_step = sp.simplify(W_num * U_step)
            h_T    = sp.inverse_laplace_transform(Y_step, s, t)
            y, _   = to_numpy(h_T, time_long)
            ax.plot(time_long, y, label=f'T2 = {Tv:g}')
        ax.set_title(f"№{link_index}. {name}\nh(t)", fontsize=9)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.legend(loc='upper right')
        graph_index += 1
        link_index += 1

axes[-1].axis('off')
axes[-2].axis('off')

fig.suptitle('Влияние коэффициента T', fontsize=12)
plt.tight_layout()
plt.show()

# Исследование влияния K
K0, T0 = 1.0, 1.0
K_vals = [K0 *0.2, K0 *0.5, K0, K0*1.5, K0*2]
tmax   = 12
time_short = np.linspace(0, tmax, 4000)
time_long = np.linspace(0, tmax*4, 4000*4)

U_step = sp.laplace_transform(sp.Heaviside(t), t, s, noconds=True)

fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flat
for i in LINKS:
    ax = axes[i - 1]
    name, W    = LINKS[i]
    for Kv in K_vals:
        if i != 5 and i != 6:
            W_num  = sp.simplify(W.subs({K: Kv, T: T0}))
            use_long = False
        elif i == 5:
            W_num  = sp.simplify(W.subs({K: Kv, T1: 2.1*T0, T2: T0}))
            use_long = True
        else:
            W_num  = sp.simplify(W.subs({K: Kv, T1: T0, T2: 2.1*T0}))
            use_long = True

        time = time_long if use_long else time_short

        Y_step = sp.simplify(W_num * U_step)
        h_T    = sp.inverse_laplace_transform(Y_step, s, t)
        y, _   = to_numpy(h_T, time)
        ax.plot(time, y, label=f'K = {Kv:g}')
    ax.set_title(f"№{i}. {name}\nh(t)", fontsize=9)
    ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
    ax.legend(loc=("upper right" if i <7 and i > 3 else "best"))

axes[-1].axis('off')

fig.suptitle('Влияние коэффициента K', fontsize=13)
plt.tight_layout()
plt.show()

# апериодическое звено (№4) — изменение T
K0, T0 = 1.0, 1.0
T_vals = [T0 / 5, T0 / 2, T0, 2 * T0, 5 * T0]
tmax   = 8 * max(T_vals)
time_T = np.linspace(0, tmax, 4000)

U_step = sp.laplace_transform(sp.Heaviside(t), t, s, noconds=True)
W_4    = LINKS[4][1]

fig, ax = plt.subplots()
ax.grid(True)
for Tv in T_vals:
    W_num  = sp.simplify(W_4.subs({K: K0, T: Tv}))
    Y_step = sp.simplify(W_num * U_step)
    h_T    = sp.inverse_laplace_transform(Y_step, s, t)
    y, _   = to_numpy(h_T, time_T)
    ax.plot(time_T, y, label=f'T = {Tv:g}')

ax.set_title('Апериодическое звено: влияние постоянной времени T')
ax.set_xlabel('t'); ax.set_ylabel('h(t)')
ax.legend()
plt.show()

# колебательное звено (№6) — изменение T1
T2_fix  = 1.0
T1_vals = [1.99, 1.5, 1.0, 0.5, 0.2, 0.05]
time_O  = np.linspace(0, 60, 6001)

U_step = sp.laplace_transform(sp.Heaviside(t), t, s, noconds=True)
W_6    = LINKS[6][1]

fig, ax = plt.subplots()
ax.grid(True)
for T1v in T1_vals:
    W_num  = sp.simplify(W_6.subs({K: 1.0, T1: T1v, T2: T2_fix}))
    Y_step = sp.simplify(W_num * U_step)
    h_T    = sp.simplify(sp.inverse_laplace_transform(Y_step, s, t))
    y, _   = to_numpy(h_T, time_O)
    ax.plot(time_O, y, label=f'T1 = {T1v:g}')

ax.axhline(1.0, ls='--', c='k', lw=0.8)
ax.set_title('Колебательное звено: T2 = 1, уменьшение T1')
ax.set_xlabel('t, c'); ax.set_ylabel('h(t)')
ax.legend()
plt.show()
