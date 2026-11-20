import control.matlab as cm
import matplotlib.pyplot as plt
import numpy as np

s = cm.tf("s")

W_open = 1 / (s**3 + 0.1 * s**2 + 0.1 * s + 1)
W_closed = 1 / (s**3 + 0.1 * s**2 + 0.1 * s + 2)

# годографы Найквиста
fig, axes = plt.subplots(1, 2, figsize=(20, 7))

                            # полный граф и регион около точки (-1, 0)
for ax, limits in zip(axes, [{"x": (-1.1, 1.1), "y": (-1.1, 1.1)}, {"x": (-1.1, -0.9), "y": (-0.1, 0.1)}]):
    cm.nyquist(W_open, ax=ax)
    ax.grid(True)
    ax.set_xlim(limits["x"])
    ax.set_ylim(limits["y"])
    ax.set_title("Годограф разомкнутой системы")
    ax.set_xlabel("Вещественная часть")
    ax.set_ylabel("Мнимая часть")

fig.savefig("nyquist.pdf", bbox_inches="tight")
plt.tight_layout()
plt.show()

# Кривая михайлова
omega = np.logspace(-3, 2, 10000)
omega = np.insert(omega, 0, 0)

_, den = cm.tfdata(W_closed)
den = np.asarray(den).ravel()

D_jw = np.polyval(den, 1j * omega)
X = D_jw.real
Y = D_jw.imag

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

                            # разные масштабы
for ax, limits in zip(axes, [{"x": (1.90, 2.1), "y": (-0.1, 0.1)}, {"x": (-5, 5), "y": (-150, 50)}, None]):
    if limits is not None:
        ax.set_xlim(limits["x"])
        ax.set_ylim(limits["y"])

    ax.plot(X, Y, 'b-', label='Кривая Михайлова')
    ax.plot(X[0], Y[0], 'ro', label=rf'Начало ($\omega=0$): ({X[0]:.1f}, {Y[0]:.1f})')

    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_title('Критерий Михайлова')
    ax.set_xlabel("Вещественная часть")
    ax.set_ylabel("Мнимая часть")
    # границы четвертей плоскости
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.legend()

fig.savefig("michael.pdf", bbox_inches="tight")
plt.tight_layout()
plt.show()

# прямая проверка устойчивости
print("Полюса разомкнутой системы")
for pole in cm.pole(W_open):
    print(pole)
print("Полюса замкнутой системы")
for pole in cm.pole(W_closed):
    print(pole)

time = np.linspace(0, 20, 1501)

fig, axes = plt.subplots(2, 2, figsize=(16, 18))

for row, W, title in zip(axes, [W_open, W_closed], ["реакция разомкнутой системы", "реакция замкнутой системы"]):
    y_h, _ = cm.step(W, T=time)
    # чтобы показать колебания до t=10
    row[0].set_xlim(-1, 11)
    row[0].set_ylim(-60, 60)
    for ax in row:
        ax.set_title(title)
        ax.set_xlabel('t, c'); ax.set_ylabel('h(t)'); ax.grid(True)
        ax.plot(time, y_h, color="b")

fig.savefig("direct_check.pdf", bbox_inches="tight")
plt.tight_layout()
plt.show()
