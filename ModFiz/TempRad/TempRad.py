# %%
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from IPython.display import display, Math, Latex
import math


# %%
t, Tp, T, U = np.loadtxt('adatok.txt', unpack=True, skiprows=1, dtype=float)

for i in range(len(T)):
    T[i]+=273.15
    Tp[i]+=273.15
    #U[i]*=0.001
    
    

d = 42 # Mérések száma

t_dat = t.reshape(-1, d)
Tp_dat = Tp.reshape(-1, d)
T_dat = T.reshape(-1, d)
U_dat = U.reshape(-1, d)


print(Tp_dat[0], '\n\n', T_dat[0], '\n\n', U_dat[0])

# %%
print(U_dat[0], '\n\n', U_dat[1]) 

# %%
def ead(t,A,b,U):
    return A*-1*np.exp(-t*b)+U

popt_arr = np.array([])

for k in range (len(t_dat)):
    fig, ax = plt.subplots() #a plotolási paramétereim kiszögezése
    fig.set_figheight(8) 
    fig.set_figwidth(10)

    popt, pcov=curve_fit(ead, t_dat[k], U_dat[k],p0=[1, 0.01, U_dat[k][0]])#1000*

    #rect = patches.Rectangle((5, 1.75), 10, 0.5, linewidth=1, edgecolor='r', facecolor='none')
    #ax.add_patch(rect)
    

    ax.scatter(t_dat[k], U_dat[k],s=10,c='b')#*1000
    ax.plot(t_dat[k], ead(t_dat[k], popt[0],popt[1], popt[2]), c='r', linewidth=5, label=r'$U(t) \,=\, - A \cdot e^{-t \cdot b} + U_0 $')#*1000
    
    ax.set_xlabel(r'$t \,\, [s]$', fontsize=20)
    ax.set_ylabel(r'$U \,\, [mV]$', fontsize=20)
    ax.legend(fontsize=12, loc='lower right')
    #ax.set_title(r'Első rugó statikus terhelése')

    ax.set_xticks(np.arange(0, 420.1, 420/7))
    
    ax.grid( which='major', color='black', linestyle=':')
    ax.minorticks_on()
    ax.grid(which='minor', color='black', linestyle=':')
    
    ax.margins(0.05, 0.05)


    perr=np.sqrt(np.diag(pcov))
    
    popt_arr = np.append(popt_arr, popt) # Paraméter-vektor készítés
    #print(k+1, ".  A, b, U:", popt)
    fig.savefig('U-t-' + str(k+1) + '.pdf', bbox_inches = "tight")

print('\n\n', popt_arr)
popt_mtrx = popt_arr.reshape(-1,3) # Paraméter-mátrix készítés

print('\n\n', popt_mtrx)


# %%
T_avrg = np.array([])
T4_avrg = np.array([])
T_err = np.array([])

for i in range(len(T_dat)):
    T_avrg = np.append(T_avrg, np.mean(T_dat[i]))
    T4_avrg = np.append(T4_avrg, T_avrg[i]**4)
    T_err = np.append(T_err, np.std(T_dat[i]))
    
print('\n\n', T_avrg) # Átlaghőmérsékletek Kelvinben
print('\n\n', T4_avrg) # Átlaghőmérsékletek negyedi hatványa Kelvinben
print('\n\n', T_err) # Hőmérsékleti hibák



# %%
(T_avrg[0])**4
print(67.15**4)

# %%
S = 6
alfa = 40
A = 0.0001

# Adott, irodalmi mennyiségek

sigma = []
U_inf = []

for i in range(len(popt_mtrx)):
    U_inf.append(popt_mtrx[i][2])
    sigma.append((alfa*U_inf[i])/(S*A*T4_avrg[i]))
    
sigma = np.array(sigma)
    
#print('\n\n', sigma) # Stefan-Boltzmann állandók

sigma_avrg = np.average(sigma)
sigma_err = np.std(sigma)

print(U_inf, '\n\n', sigma, '\n\n', sigma_avrg, '\n\n', sigma_err) # Átlag és szórás

# %%
sigma4 = np.array([])
sigma1 = np.array([])

for i in range(len(U_dat)):
    sigma4 = np.append(sigma4, ((S*A)/alfa)*((U_dat[1][i]-U_dat[0][i])/((T_dat[1][i])**4-(T_dat[0][i])**4)))
    sigma1 = np.append(sigma1, ((S*A)/alfa)*((U_dat[1][i]-U_dat[0][i])/((T_dat[1][i])-(T_dat[0][i]))))
    
print('\n\n', sigma4, '\n\n', sigma1) # Stefan-Boltzmann állandók
print('\n\n', np.average(sigma4), '\n\n', np.std(sigma4), '\n\n', np.average(sigma1), '\n\n', np.std(sigma1)) # Átlag és szórás 

# %%
P=[]

for i in range(len(sigma)):
    P.append(sigma[i]*T4_avrg[i])

P=np.array(P)

print('\n\n', P) # Sugárzási teljesítmények

# %%
def exponencialis(x,a,c):
    return(a * np.power(x,c))
def exponencialis2(x,c):
    return(c * np.power(x,4))

popt3,pcov3=curve_fit(exponencialis,T_avrg,P)
perr3 = np.sqrt(np.diag(pcov3))

popt4,pcov4=curve_fit(exponencialis2,T_avrg,P)
perr4 = np.sqrt(np.diag(pcov4))



fig, ax = plt.subplots() #a plotolási paramétereim kiszögezése
fig.set_figheight(8) 
fig.set_figwidth(10)


ax.plot(T_avrg, P, "o", color="blue")
ax.plot(T_avrg, exponencialis(T_avrg,*popt3),label='Legjobban illeszkedő egyenes', color='limegreen', linewidth=2)
ax.plot(T_avrg, exponencialis2(T_avrg,*popt4),label='Illesztett elméleti egyenes', color='red', linewidth=2)


ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlabel(r"$T \,\, [K]$", fontsize=20)
ax.set_ylabel(r"$P \,\, [\frac{W}{m^2}]$", fontsize=20)
#title(r"A $T-P$ grafikon")
ax.legend(fontsize=12, loc='lower right')

ax.set_xticks(np.arange(0, 420.1, 420/7))
    
ax.grid( which='major', color='black', linestyle='-')
ax.minorticks_on()
ax.grid(which='minor', color='black', linestyle=':')
    
ax.margins(0.05, 0.05)

print (['a, c'])
print (*popt3)
print (*perr3)
print('\n')
print(np.log(popt3[0]), np.log(popt3[1]))

print(1/popt3[0]*perr3[0], 1/popt3[1]*perr3[1])
print('\n')

print (['c'])
print (*popt4)
print (*perr4)
print('\n')
print(np.log(popt4[0]))
print(1/popt4[0]*perr4[0])



fig.savefig('touching-lines.pdf', bbox_inches = "tight")

# %%
#stefan-boltzmann állandó

sig = np.exp(np.log(popt4[0]))
dsig = 1/popt4[0]*perr4[0]*sig

print('\n\n', sig, dsig)

# %%



