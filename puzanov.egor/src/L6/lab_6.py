import control.matlab as cm
import matplotlib.pyplot as plt
import numpy as np
from pidtune import pidtune

s = cm.tf("s")
sigma_req = 5
t_p_req = 15

W_open = 1 / (((s + 0.1) ** 2) * (s + 1))
W_closed = cm.feedback(W_open, 1)

# 1. Оценка исходной замкнутой системы
print("Полюса замкнутой системы")
for pole in cm.pole(W_closed):
    print(pole)

[[[a0, a1, a2, a3]]] = W_closed.den
matrix = np.array(
    [[a1, a3, 0],
     [a0, a2, 0],
     [0, a1, a3]]
)
def getMinor(matrix, number):
    return np.linalg.det(matrix[:number, :number])

print(f"Определители Гурвица: D1={getMinor(matrix, 1)} D2={getMinor(matrix, 2)} D3={getMinor(matrix, 3)}\n")

stepinfo = cm.stepinfo(W_closed)
print(f"t_п = {stepinfo['SettlingTime']}, sigma = {stepinfo['Overshoot']}, ошибка = {cm.dcgain(cm.feedback(1, W_open))}\n")

# 2. Настройка ПИД методом Зиглера-Николса
time = np.linspace(0, 100, 10001)
y_h, T = cm.step(W_open, T=time)

# Нахождение точки перегиба
dy_dx = np.gradient(y_h, T)
d2y_dx2 = np.gradient(dy_dx, T)

zero_idx = np.where(np.diff(np.sign(d2y_dx2)))[0][0]
t_flip = np.interp(0, [d2y_dx2[zero_idx], d2y_dx2[zero_idx+1]], [T[zero_idx], T[zero_idx+1]])
h_flip = np.interp(t_flip, T, y_h)
point_flip = (t_flip, h_flip)

# kx+b
k = np.interp(t_flip, T, dy_dx)
b = h_flip - k*t_flip

stepinfo = cm.stepinfo(W_open)
L = -b/k
T = (stepinfo["SteadyStateValue"]-b)/k - L

fig, axes = plt.subplots(1, 2, figsize=(19, 9))
for ax, limits in zip(axes.flat, [{"x": (-2, 20), "y": (-2, 58)}, None]):
    if limits is not None:
        ax.set_xlim(limits["x"])
        ax.set_ylim(limits["y"])

    ax.set_title("Реакция системы")
    ax.set_xlabel("t, c"); ax.set_ylabel("h(t)"); ax.grid(True)
    ax.plot(time, y_h, label="W_откр")
    ax.axline(xy1=point_flip, slope=k, c="g", linestyle="--")
    ax.axvline(T+L, linestyle="--",c="r")
    ax.axhline(stepinfo["SteadyStateValue"], linestyle="--",c="gray")
    ax.axhline(0, linestyle="--",c="gray")
    ax.plot(*point_flip, 'go', label=f'перегиб: ({point_flip[0]:.4f}, {point_flip[1]:.4f})')
    ax.plot(L, 0, "C1o", label=f'L = {L}')
    ax.plot(T+L, stepinfo["SteadyStateValue"], "ro", label=f'T = {T}')
    ax.legend()

fig.savefig("task_2.pdf", bbox_inches="tight")

Kp = 1.2*T/L 
Ti = 2*L 
Td = 0.5*L
print(f'Метод 1: L = {L}, T = {T} c')
print(f'Kp = {Kp}, Ti = {Ti} c, Td = {Td} c\n')

#  3. Оценка системы с ПИД регулятором
k1 = Kp; k2 = Kp/Ti; k3 = Kp*Td
C = k1 + k2/s + k3*s
W_open_fixed = C*W_open
W_closed_fixed = cm.feedback(W_open_fixed, 1)
stepinfo = cm.stepinfo(W_closed_fixed)
print(f'C(s) = {k1} + {k2}/s + {k3}*s')
print(
    f"t_п = {stepinfo['SettlingTime']} (<= {t_p_req}), "
    f"sigma = {stepinfo['Overshoot']} (<= {sigma_req}), "
    f"ошибка = {cm.dcgain(cm.feedback(1, W_open_fixed))}\n"
)

# Переходные характеристики до и после
time = np.linspace(0, 100, 1001)

fig, axes = plt.subplots(1, 2, figsize=(19, 9))
for ax, limits in zip(axes.flat, [{"x": (-1, 30), "y": (-2, 3)}, None]):
    if limits is not None:
        ax.set_xlim(limits["x"])
        ax.set_ylim(limits["y"])

    for name, W in zip(
        ["Исходная W", "Скорректированная W"], [W_closed, W_closed_fixed]
    ):
        ax.set_title("Реакция системы")
        ax.set_xlabel("t, c"); ax.set_ylabel("h(t)"); ax.grid(True)
        y_h, _ = cm.step(W, T=time)
        ax.plot(time, y_h, label=name)
        ax.legend()

fig.savefig("task_3.pdf", bbox_inches="tight")

# 4. Автоматическая настройка с помощью pidtune
C_auto, stats = pidtune(W_open)
W_open_auto = C_auto*W_open
W_closed_auto = cm.feedback(W_open_auto, 1)
stepinfo = cm.stepinfo(W_closed_auto)
print(
    "Автоматическая настройка: "
    f"Kp={stats["Kp"]} Ki={stats["Ki"]} Kd={stats["Kd"]} | "
    f"t_п = {stepinfo['SettlingTime']}, "
    f"sigma = {stepinfo['Overshoot']}, "
    f"ошибка = {cm.dcgain(cm.feedback(1, W_open_auto))}\n"
)

time = np.linspace(0, 100, 1001)
fig, ax= plt.subplots(1, 1, figsize=(8, 9))
color = "C1"
for name, W in zip(
    ["Метод Зиглера—Николса", "Автоматичекая настройка"], [W_closed_fixed, W_closed_auto]
):
    ax.set_title("Реакция системы"); ax.set_xlabel("t, c"); ax.set_ylabel("h(t)"); ax.grid(True)
    y_h, _ = cm.step(W, T=time)
    ax.plot(time, y_h, label=name, c=color)
    ax.legend()
    color = "C2"

fig.savefig("task_4.pdf", bbox_inches="tight")

[Gm_manual, Pm_manual, _, _] = cm.margin(W_open_fixed)
Gm_manual = 20*np.log10(Gm_manual)
[Gm_auto, Pm_auto, _, _] = cm.margin(W_open_auto)
Gm_auto = 20*np.log10(Gm_auto)
print("запасы устойчивости:")
print(f"Метод Зиглера—Николса - запас по амплитуде: {Gm_manual} дБ, запас по фазе: {Pm_manual} град.")
print(f"Автонастройка - запас по амплитуде: {Gm_auto} дБ, запас по фазе: {Pm_auto} град.")

# чтобы понять почему запасы по амплитуде бесконечны
omega = np.logspace(-3, 4, 10000)
fig, (mag_ax, phase_ax) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
for W, label, i in [(W_open_fixed, "Метод Зиглера—Николса", 1), (W_open_auto, "Автонастройка", 2)]:
    cm.bode([W], omega=omega, dB=True, ax=[mag_ax, phase_ax], grid=True, color=f"C{i}", label=label)
    mag_ax.set_title("ЛАЧХ", fontsize=9)
    mag_ax.set_ylabel("аплитуда, дБ", fontsize=9)
    mag_ax.set_xlabel("частота, рад/с", fontsize=9)
    phase_ax.set_title("ЛФЧХ", fontsize=9)
    phase_ax.set_ylabel("фаза, град", fontsize=9)
    phase_ax.set_xlabel("частота, рад/с", fontsize=9)
fig.savefig("task_additional.pdf", bbox_inches="tight")
