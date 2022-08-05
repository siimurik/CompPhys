import pandas as pd
import matplotlib.pyplot as plt

# Setting background for graph
plt.style.use('dark_background')

# Reading in the general results
df_euler    = pd.read_table('result1_euler.dat', skiprows=4, delim_whitespace=True,
                                names=['t', 'x(t)'] )
df_pred_cor = pd.read_table('result1_pred_cor.dat', skiprows=4, delim_whitespace=True,
                                names=['t', 'x(t)'] )
df_rk       = pd.read_table('result1_rk.dat', skiprows=4, delim_whitespace=True,
                                names=['t', 'x(t)'] )

# Printing the values from dataframe
#print(df.head())

# Assigning names to datacolumns
t_e =  df_euler['t']
x_e =  df_euler['x(t)']
t_pc = df_pred_cor['t']
x_pc = df_pred_cor['x(t)']
t_rk = df_rk['t']
x_rk = df_rk['x(t)']

# Plotting x(t), v(t) and E
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t_e, x_e,   label="Euler")
ax.plot(t_pc, x_pc, label="Predictor-Corrector")
ax.plot(t_rk, x_rk, label="Runge-Kutta")
ax.set_xlabel('Time')
ax.set_ylabel('Position')
plt.title("The Carrying Capacity Problem")
plt.grid()
plt.legend()
plt.show()