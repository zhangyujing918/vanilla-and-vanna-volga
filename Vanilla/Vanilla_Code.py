import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

class BS:
    
    def __init__(self, S, Klist, sigma, r, T, option_type):
        self.S = S
        self.Klist = Klist
        self.sigma = sigma
        self.r = r
        self.T = T
        self.option_type = option_type
        
    
    def price(self):
        prices = []
        for K in self.Klist:
            d1 = (np.log(self.S / K) + (self.r + self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            d2 = d1 - self.sigma * np.sqrt(self.T)
            
            if self.option_type.lower() == 'call':
                prices.append(self.S * norm.cdf(d1) - K * np.exp(-self.r * self.T) * norm.cdf(d2))
                
            elif self.option_type.lower() == 'put':
                prices.append(K * np.exp(-self.r * self.T) * norm.cdf(-d2) - self.S * norm.cdf(-d1))
        
        return prices
        
    
    def delta(self):
        deltas = []
        for K in self.Klist:
            d1 = (np.log(self.S / K) + (self.r + self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            if self.option_type.lower() == 'call':
                deltas.append(norm.cdf(d1))
            elif self.option_type.lower() == 'put':
                deltas.append(-norm.cdf(-d1))
        return deltas
        
    
    def gamma(self):
        gammas = []
        for K in self.Klist:
            d1 = (np.log(self.S / K) + (self.r + self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            gammas.append(norm.pdf(d1) / (self.S * self.sigma * np.sqrt(self.T)))
        return gammas
    
    
    def vega(self):
        vegas = []
        for K in self.Klist:
            d1 = (np.log(self.S / K) + (self.r + self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            vegas.append(self.S * norm.pdf(d1) * np.sqrt(self.T))
        return vegas
    
    
    def theta(self):
        thetas = []
        for K in self.Klist:
            d1 = (np.log(self.S / K) + (self.r + self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            d2 = d1 - self.sigma * np.sqrt(self.T)
            if self.option_type.lower() == 'call':
                thetas.append(-(self.S * norm.pdf(d1) * self.sigma) / (2 * np.sqrt(self.T)) - K * np.exp(-self.r * self.T) * norm.cdf(d2))
            elif self.option_type.lower() == 'put':
                thetas.append(-(self.S * norm.pdf(d1) * self.sigma) / (2 * np.sqrt(self.T)) + K * np.exp(-self.r * self.T) * norm.cdf(-d2))
        return thetas
        
    
    def rho(self):
        rhos = []
        for K in self.Klist:
            d2 = (np.log(self.S / K) + (self.r - self.sigma ** 2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
            if self.option_type.lower() == 'call':
                rhos.append(K * self.T * np.exp(-self.r * self.T) * norm.cdf(d2))
            elif self.option_type.lower() == 'put':
                rhos.append(-K * self.T * np.exp(-self.r * self.T) * norm.cdf(-d2))
        return rhos


####################################################
# # 随strike变化

strike = np.linspace(80, 120, 100)
o = BS(100, strike, 0.15, 0.02, 1/12, 'call')
o2 = BS(100, strike, 0.15, 0.02, 3/12, 'call')
o3 = BS(100, strike, 0.15, 0.02, 6/12, 'call')
fig, axs = plt.subplots(2, 3, figsize = (18, 12), dpi = 180)

axs[0, 0].plot(strike, o.price(), color='b')
axs[0, 0].plot(strike, o2.price(), color='r')
axs[0, 0].plot(strike, o3.price(), color='green')
axs[0, 0].set_title('Option Price vs. Strike Price', fontsize = 12)
axs[0, 0].legend(['1M', '3M', '6M'])

axs[0, 1].plot(strike, o.delta(), color='b')
axs[0, 1].plot(strike, o2.delta(), color='r')
axs[0, 1].plot(strike, o3.delta(), color='green')
axs[0, 1].set_title('Delta vs. Strike Price',  fontsize = 12)
axs[0, 1].legend(['1M', '3M', '6M'])

axs[0, 2].plot(strike, o.gamma(), color='b')
axs[0, 2].plot(strike, o2.gamma(), color='r')
axs[0, 2].plot(strike, o3.gamma(), color='green')
axs[0, 2].set_title('Gamma vs. Strike Price',  fontsize = 12)
axs[0, 2].legend(['1M', '3M', '6M'])

axs[1, 0].plot(strike, o.theta(), color='b')
axs[1, 0].plot(strike, o2.theta(), color='r')
axs[1, 0].plot(strike, o3.theta(), color='green')
axs[1, 0].set_title('Theta vs. Strike Price', fontsize = 12)
axs[1, 0].legend(['1M', '3M', '6M'])

axs[1, 1].plot(strike, o.vega(), color='b')
axs[1, 1].plot(strike, o2.vega(), color='r')
axs[1, 1].plot(strike, o3.vega(), color='green')
axs[1, 1].set_title('Vega vs. Strike Price', fontsize = 12)
axs[1, 1].legend(['1M', '3M', '6M'])

axs[1, 2].plot(strike, o.rho(), color='b')
axs[1, 2].plot(strike, o2.rho(), color='r')
axs[1, 2].plot(strike, o3.rho(), color='green')
axs[1, 2].set_title('Rho vs. Strike Price', fontsize = 12)
axs[1, 2].legend(['1M', '3M', '6M'])

plt.suptitle('Option Greeks vs. Strike Price -- Call', fontsize=16)
plt.show()










