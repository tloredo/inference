
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#x = np.linspace(98.42, 99.21, 100)
x = np.linspace(147.63, 148.31, 100)
y = np.random.random((len(x)))

def func(x, pos):
    if pos is None:
        return ''
    else:
        return "%.2f" % (x)

ax = plt.subplot(111)

formatter_x = FuncFormatter(func)
ax.xaxis.set_major_formatter(formatter_x)

plt.plot(x, y)
plt.xlim(np.min(x), np.max(x))
plt.show()


