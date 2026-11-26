from pprint import pprint

import control.matlab as cm
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

#  1. Исходные данные
n = 2
A = np.array([[0.5, -10], [1, -2]])
B = np.array([0, 1])
C = np.array([3, 1])
M = np.array([[4, 0], [-2, 0.5]])

# 2. Замена базиса
M_inv = np.linalg.inv(M)
A_new = M_inv @ A @ M
B_new = M_inv @ B
C_new = C @ M

print("M_inv:")
pprint(M_inv)
print("A_new:")
pprint(A_new)
print("B_new:")
pprint(B_new)
print("C_new:")
pprint(C_new)

# 3. Переход к передаточной функции
s = sp.symbols("s")
A_s = sp.Matrix(A)
B_s = sp.Matrix(B).reshape(n, 1)
C_s = sp.Matrix(C).reshape(1, n)

sI_A = s * sp.eye(n) - A_s
W_sym = sp.simplify(C_s * sI_A.inv() * B_s)[0, 0]
print(f"\n\nW(s) аналитически: {W_sym}")

W_tf = cm.ss2tf(A_new, B_new, C_new, 0)
print(f"W(s) ss2tf: {W_tf}\n\n")

# 4. Моделирование переходного процесса
time = np.linspace(0, 10, 1001)
sys_old = cm.ss(A, B, C, 0)
sys_new = cm.ss(A_new, B_new, C_new, 0)
fig, axes = plt.subplots(3, 1, figsize=(16, 18))

for ax, W, title in zip(
    axes.flat,
    [sys_old, sys_new, W_tf],
    ["Исходная модель", "Модель в новом базисе", "Передаточная функция W(s)"],
):
    y_h, _ = cm.step(W, T=time)
    ax.set_title(title)
    ax.set_xlabel("t, c")
    ax.set_ylabel("h(t)")
    ax.grid(True)
    ax.plot(time, y_h)

fig.savefig("comparison.pdf", bbox_inches="tight")

# 5. Управляемость и Наблюдаемость
U = cm.ctrb(A, B)
rank_U = np.linalg.matrix_rank(U)
print("\nМатрица управляемости:")
pprint(U)
print(f"ранг: {rank_U}")

N = cm.obsv(A, C)
rank_N = np.linalg.matrix_rank(N)
print("Матрица наблюдаемости:")
pprint(N)
print(f"ранг: {rank_N}")

# 6. Канонические формы
[[[_, a1, a0]]] = W_tf.den
[[[b1, b0]]] = W_tf.num

A_c = np.array([[0, -a0], [1, -a1]])
B_c = np.array([1, 0])
C_c = np.array([b1, b0 - a1 * b1])
sys_c = cm.ss(A_c, B_c, C_c, 0)

A_o = np.array([[0, 1], [-a0, -a1]])
B_o = np.array([b1, b0 - a1 * b1])
C_o = np.array([1, 0])
sys_o = cm.ss(A_o, B_o, C_o, 0)

time = np.linspace(0, 10, 1001)
fig, axes = plt.subplots(2, 1, figsize=(16, 18))

for ax, W, title in zip(
    axes.flat,
    [sys_c, sys_o],
    ["Каноническая управляемая форма", "Каноническая наблюдаемая форма"],
):
    y_h, _ = cm.step(W, T=time)
    ax.set_title(title)
    ax.set_xlabel("t, c")
    ax.set_ylabel("h(t)")
    ax.grid(True)
    ax.plot(time, y_h)

fig.savefig("canonical.pdf", bbox_inches="tight")
