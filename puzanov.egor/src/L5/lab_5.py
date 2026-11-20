import control.matlab as cm
import matplotlib.pyplot as plt
import numpy as np

s = cm.tf("s")

W_open = 100 / (s * (0.1 * s + 1) * (0.02 * s + 1) * (0.01 * s + 1))
W_closed = cm.feedback(W_open, 1)

# 1. Проверка устойчивости по Гурвицу
[[[a0, a1, a2, a3, a4]]] = W_closed.den
matrix = np.array(
    [
        [a1, a3, 0, 0],
        [a0, a2, a4, 0],
        [0, a1, a3, 0],
        [0, a0, a2, a4],
    ]
)

def getMinor(matrix, number):
    return np.linalg.det(matrix[:number, :number])

print(
    "Определители Гурвица:  "
    f"D1={getMinor(matrix, 1)}  "
    f"D2={getMinor(matrix, 2)}  "
    f"D3={getMinor(matrix, 3)}  "
    f"D4={getMinor(matrix, 4)}"
)

# Исходные запасы устойчивости
[Am, Pm, Wcg, Wcp] = cm.margin(W_open)
print(f"start - Am: {Am}, Pm:{Pm}, Wcg: {Wcg}, Wcp: {Wcp}")

time = np.linspace(0, 20, 2001)
y_h, _ = cm.step(W_closed, T=time)

fig, axes = plt.subplots(3, 1, figsize=(16, 18))

ax1 = axes.flat[0]
ax1.set_title("реакция закрытой системы")
ax1.set_xlabel("t, c")
ax1.set_ylabel("h(t)")
ax1.grid(True)
ax1.plot(time, y_h)

ax2 = axes.flat[1]
ax3 = axes.flat[2]
cm.bode(W_open, dB=True, ax=[ax2, ax3], grid=True)
ax2.set_title("ЛАЧХ открытой системы")
ax2.set_ylabel("аплитуда, дБ")
ax2.set_xlabel("частота, рад/с")
ax3.set_title("ЛФЧХ открытой системы")
ax3.set_ylabel("фаза, град")
ax3.set_xlabel("частота, рад/с")

fig.savefig("task_1.pdf", bbox_inches="tight")

# 2. Расчёт параметров желаемой ЛАЧХ
Tpp_req = 0.8
M_req = 1.3
sigma_req = 25

gamma_req = np.degrees(np.arcsin(1 / M_req))
wc = 2 * np.pi / Tpp_req
mag_i, ph_i, _ = cm.bode(W_open, [wc], plot=False)
mag_i = mag_i[0]
ph_i = np.degrees(ph_i[0])
gamma_base = 180 + ph_i

print(mag_i, ph_i, gamma_base, gamma_req)

# 3. Рассчёт ПКУ
# W_ku(s) = K * (T1*s + 1) / (T2*s + 1)
delta_phi = gamma_req - gamma_base
alpha = (1 + np.sin(np.radians(delta_phi))) / (1 - np.sin(np.radians(delta_phi)))
# w_m = 1/sqrt(T1*T2) = wc   =>   T1*T2 = 1/wc^2,  T1 = alpha*T2
T2 = 1 / (wc * np.sqrt(alpha))
T1 = alpha * T2

# Определяем K из условия L_желаемое(w) = 0
L_wc = 20 * np.log10(mag_i)
# Амплитудный вклад звена (T1 s+1)/(T2 s+1) на wc
L_zv_wc = 20 * np.log10(np.sqrt(1 + (T1 * wc) ** 2) / np.sqrt(1 + (T2 * wc) ** 2))
L_K = 0 - L_wc - L_zv_wc
K = 10 ** (L_K / 20)

# Передаточная функция КУ
W_ku = K * cm.tf([T1, 1], [T2, 1])

# 4. Скорректированная САУ
W_final = W_ku * W_open
W_closed_final = cm.feedback(W_final, 1)

[Am2, Pm2, Wcg2, Wcp2] = cm.margin(W_final)
print(f"fin - Am: {Am}, Pm:{Pm}, Wcg: {Wcg}, Wcp: {Wcp}")
M_fact = 1 / np.degrees(np.sin(Pm2))
print(f"Фактическое M = 1/sin(Pm) = {M_fact} (требуется <= {M_req})")

# 5. Показатели качества переходного процесса
info = cm.stepinfo(W_closed_final)
Tpp_fact = info["SettlingTime"]
sigma_fact = info["Overshoot"]

print(f"Соответствие: Tпп {Tpp_fact} <= {Tpp_req}, sigma {sigma_fact} <= {sigma_req}")

# 6. ЛАЧХ/ЛФЧХ исходной и скорректированной САУ
fig, axes = plt.subplots(2, 1, figsize=(26, 18))
for name, W in zip(["Исходная W", "W_ку", "Скорректированная W"], [W_open, W_ku, W_final]):
    [ax1, ax2] = axes.flat
    cm.bode(W, dB=True, ax=[ax1, ax2], grid=True, label=name)
    ax1.set_title(f"{name}: ЛАЧХ")
    ax1.set_ylabel("аплитуда, дБ")
    ax1.set_xlabel("частота, рад/с")
    ax1.legend()
    ax2.set_title(f"{name}: ЛФЧХ")
    ax2.set_ylabel("фаза, град")
    ax2.set_xlabel("частота, рад/с")
    ax2.legend()
    
fig.savefig("task_6.1.pdf", bbox_inches="tight")

#  Переходные характеристики до и после
time = np.linspace(0, 10, 1001)
fig, axes = plt.subplots(2, 1, figsize=(16, 18))
for ax, name, W in zip(axes.flat, ["Исходная W", "Скорректированная W"], [W_closed, W_closed_final]):
    ax.set_title(f"{name}: реакция системы")
    ax.set_xlabel("t, c")
    ax.set_ylabel("h(t)")
    ax.grid(True)
    y_h, _ = cm.step(W, T=time)
    ax.plot(time, y_h, color="b")

fig.savefig("task_6.2.pdf", bbox_inches="tight")
