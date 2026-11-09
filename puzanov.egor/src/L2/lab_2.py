import control.matlab as cm
from control import interconnect
import matplotlib.pyplot as plt
import numpy as np

s = cm.tf("s")

W1 = 1/(s+1)
W2 = cm.tf([2], [1]) # создание объекта TransferFunction в виде константы 2
W3 = 5/(s**2+3*s+2)

# Аналитические предобразования из отчета
W4 = W3/(1+W3)
W5 = W2*W4
W6 = W5/(1+W5)
W7 = W1+1/W2
W8 = W7*W6
W9 = W8/(1+W8)

# проверка через функции матлаба
"""
Выразить через функции series, parallel и feedback у меня получалось только выполняя те же преобразования,
что и выше, поэтому выполнять одну и ту же работу дважды мне показалось бесполезно.
Я нашел способ создание схемы через задание звеньев с наименованием их связей - чтобы было более
понятно тут нарисована схема какие связи как называются, а также какими переменными в коде обозначены
сами звенья и сумматоры:

!                                                 a  (продолжение)
!                                                 ↓
!                     ┌───────────────────────────────────────────────────┐
!        SM1          |                   SM2      с              d       | SM3
! u       v           |                    v       ↓              ↓      +| ↙                   y
! ────────o───────────┴──>|W1|────────────>o─────────────>|W2|───────────>o─────────────|W3|────┬───>
!        -|   ↑            ^        ↑     -|               ^             -|    ↑         ^      |
!         |   a        W1_check     b      |            W2_check          |    e      W3_check  |
!         |                                |                              |                     |
!         |              y (продолжение) → └──────────────────────────────┴─────────────────────┤
!         └─────────────────────────────────────────────────────────────────────────────────────┘
!                                                ↑
!                                                y (продолжение)
"""
SM1 = cm.summing_junction(inputs=['u', '-y'], output='a')
W1_check = cm.tf(W1, inputs='a', outputs='b')
SM2 = cm.summing_junction(inputs=['b', '-y'], output='c')
W2_check = cm.tf(W2, inputs='c', outputs='d')
SM3 = cm.summing_junction(inputs=['d', 'a', '-y'], output='e')
W3_check = cm.tf(W3, inputs='e', outputs='y')

W_final = interconnect(
    [W1_check, W2_check, W3_check, SM1, SM2, SM3],
    inputs='u',
    outputs='y'
)
print(W_final)

# сравнение моделей
time = np.linspace(0, 10, 1001)

y_h, _ = cm.step(W9, T=time)
y_h_check, _ = cm.step(W_final, T=time)
y_w, _ = cm.impulse(W9, T=time)
y_w_check, _ = cm.impulse(W_final, T=time)
y_s, _, _ = cm.lsim(W9, U=np.sin(time), T=time)
y_s_check, _, _ = cm.lsim(W_final, U=np.sin(time), T=time)

fig, axes = plt.subplots(3, 2, figsize=(18, 12))
axes = axes.flat

GRAPHS = [
    ("скачок (аналитические преобразования)", y_h),
    ("скачок (проверка схемой)", y_h_check),
    ("импульс (аналитические преобразования)", y_w),
    ("импульс (проверка схемой)", y_w_check),
    ("синусоида (аналитические преобразования)", y_s),
    ("синусоида (проверка схемой)", y_s_check),
]

for i in range(len(GRAPHS)):
    ax = axes[i]
    ax.set_title(GRAPHS[i][0])
    ax.set_xlabel('t, c'); ax.set_ylabel('y(t)'); ax.grid(True)
    ax.plot(time, GRAPHS[i][1], color=("b" if i%2==0 else "g"), lw=1.8)

fig.suptitle('Проверка аналитических преоброзований')
plt.tight_layout()
plt.show()

