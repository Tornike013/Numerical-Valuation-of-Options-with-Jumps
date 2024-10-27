# -*- coding: utf-8 -*-
"""
Created on Fri May 10 22:15:31 2024

@author: User
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:12:34 2024

@author: User
"""


import math
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

i = 5
lambda_ = 0.1
sigma = 0.2
sigma_j = 0.5
K = 1
T = 1
t = 0.8
tau_= 0.2
r = 0  
n = 65
zeta = math.exp((sigma_j ** 2) / 2) - 1
lambda_der = lambda_ * (1 + zeta)   


def w(t, s, K, tau_, r, sigma, sigma_j, i):
    total_sum = 0
    for m in range(0, i-1):
        sigma_m = sigma ** 2 + (m * sigma_j ** 2) / tau_
        r_sub_m = r - lambda_ * zeta + (m * math.log(1 + zeta)) / tau_
        total = ((math.exp(-lambda_der * tau_) * (lambda_der * tau_) ** m) / math.factorial(m)) * black_scholes(s, K, tau_, r_sub_m, sigma_m)
        total_sum += total
    return total_sum

# Define the Black-Scholes function
def black_scholes(s, K, tau_, r, sigma_m, option_type='call'):
    d1 = (np.log(s / K) + (r + 0.5 * sigma_m ** 2) * tau_) / (np.sqrt(sigma_m) * np.sqrt(tau_))
    d2 = d1 - np.sqrt(sigma_m) * np.sqrt(tau_)

    if option_type == 'call':
        option_price = s * norm.cdf(d1) - K * np.exp(-r * tau_) * norm.cdf(d2)

    return option_price



x_star = 4
s = np.exp(np.linspace(-x_star, x_star, n))
analytical_solution = []

for price in s:
    each_price = w(t, price, K, tau_, r, sigma, sigma_j, i)
    analytical_solution.append(each_price)
    
print(analytical_solution)
    
    
    



def merton_numerical_scheme(T, x_star, n, q, sigma, r, lambda_, zeta, K, mu_J=0, sigma_j=0.5):
    # Discretize space and time
    h = 2 * x_star / (n - 1)
    k = 0.2  
    price = np.linspace(-x_star, x_star, n) 
    
    
    u = np.zeros((q, n))
    for i in range(n): 
        if math.exp(price[i]) - K >= 0:
            u[0, i] = math.exp(price[i]) - K
    
    
  
    def epsilon(tau_m, x_i, x_star):
        term1 = np.exp(x_i + 0.5 * sigma_j**2) * norm.cdf((x_i -x_star + sigma_j**2) / sigma_j)
        term2 = K * np.exp(-r * tau_m) * norm.cdf((x_i - x_star) /sigma_j)
        return term1 - term2
    
    
    
    def f(x):
        return np.exp(-(x - mu_J)**2 / (2 * sigma_j**2)) / np.sqrt(2 * np.pi * sigma_j**2)
    
    
    
    for m in range(1, q):
        tau_m = (m-1)*k
        C = np.zeros((n, n))
        D = np.zeros((n, n))
        b = np.zeros(n)
        omega0 = 3/2 if m > 1 else 1
        omega1 = 2 if m > 1 else 1
        omega2 = -1/2 if m > 1 else 0
        
        for i in range(1, n - 1):
            C[i, i - 1] = -k * sigma**2 / (2 * h**2) + k * (r - 0.5 * sigma**2 - lambda_ * zeta) / (2 * h)
            C[i, i] = k * sigma**2 / h**2 + k * (r + lambda_)
            C[i, i + 1] = -k * sigma**2 / (2 * h**2) - k * (r - 0.5 * sigma**2 - lambda_ * zeta) / (2 * h)
            D[i, 0] = -k * h * lambda_ * f(price[0] - price[i]) / 2 
            D[i, -1] = -k * h * lambda_ * f(price[-1] - price[i]) / 2
            
            for j in range(1, n - 1):
                D[i, j] = -k * h * lambda_ * f(price[j] - price[i])
            
            b[i] = k * lambda_ * epsilon(tau_m, price[i], x_star) + omega1 * u[m - 1, i] + omega2 * u[m - 2, i]
        
        C += omega0 * np.eye(n) + D
        C[0, 0], C[-1, -1] = 1, 1 
        b[0], b[-1] = 0, omega0 * (np.exp(x_star) - K * np.exp(-r *tau_m))
        u[m, :] = np.linalg.solve(C, b)
        
    return u



T=1                  
x_star = 4 
k = 0.2
q = int(0.2/k) + 1
sigma = 0.2
r = 0                
lambda_ = 0.1           
zeta = math.exp(sigma_j ** 2) / 2 - 1             
K = 1

u = merton_numerical_scheme(T, x_star, n, q, sigma, r, lambda_, zeta, K)

print(u[-1])




differences = analytical_solution-u[-1]
plt.plot(np.linspace(-4, 4, 65), differences)


