import numpy as np
import sympy as sp
from sympy import symbols,Eq,solve
from scipy.stats import norm
N = norm.cdf
n = norm.pdf


------------Mid1------------------
------------bid-------------------
def find_bid(probability):
    o = symbols("o")
    eq1 = o/(1+o)
    re = solve(Eq(eq1,probability),o)
    print(solve(Eq(eq1,probability),o))
    return re
find_bid(9/100)

______________bet problem_________________
def bet_cal(bid1,bid2,m,maxlose):
    x = symbols("x")
    eq1 = (-x+bid2*(m-x))
    eq2 = (bid1*x-(m-x))
    profit = (m*(bid1*bid2-1))/(2+bid1+bid2)
    n = m**2 / profit
    print("bond1:",solve(Eq((eq1),0)))
    print("bond2:",solve(Eq((eq2),0)))
    print("x:",solve(Eq(eq1,eq2),x))
    print("profit:",profit)
    print("How much need to get",m,":",n)
    print("lose bond1:",solve(Eq(eq1,maxlose)))
    print("lose bond2:",solve(Eq(eq2,maxlose)))
    return
bet_cal(5,1.4,100,-10)


________________Horse bet______________
def Horse_bet(first_bet, second_bet, third_bet, earn_in_totalbet_x_percent):
    o = symbols("o")
    eqa = -first_bet*o + second_bet + third_bet
    eqb = first_bet - second_bet*o + third_bet
    eqc = first_bet + second_bet - o*third_bet
    oa = solve(Eq(eqa, earn_in_totalbet_x_percent), o)
    ob = solve(Eq(eqb, earn_in_totalbet_x_percent), o)
    oc = solve(Eq(eqc, earn_in_totalbet_x_percent), o)
    print("O_a:", oa)
    print("O_b:", ob)
    print("O_c:", oc)
    return
Horse_bet(64000,27000,9000,10000)

____Solve for Interest Rate Problem______
print(np.log(79.5/110)/-4)


__________Mean_and_Variance______________
x = [(4/14),(7/14),(3/14)]
r = [7,10,14]
def Mean_Variance(x, r):
    Ex = []
    varx = []
    Re = []
    for i,j in enumerate(x):
        Ex.append(j * r[i])
    for i,j in enumerate(x):
        varx.append(j * r[i]**2)
    Varx = sum(varx) - sum(Ex)**2
    E_x = sum(Ex)
    Re.append(E_x)
    Re.append(Varx)
    print (Re)
    return Re
Mean_Variance(x,r)

___________Tylar_Expansion__________
x,y = sp.symbols("x,y")
#Change Equation Here
f = sp.exp(x)*sp.log(1+y)
def second_ty(f):
    fx = sp.diff(f, x)
    fy = sp.diff(f, y)
    fxx = sp.diff(f, x, x)
    fxy = sp.diff(f, x, y)
    fyy = sp.diff(f, y, y)
    x0 = 0
    y0 = 0
    f_val = f.subs([(x, x0), (y, y0)])
    fx_val = fx.subs([(x, x0), (y, y0)])
    fy_val = fy.subs([(x, x0), (y, y0)])
    fxx_val = fxx.subs([(x, x0), (y, y0)])
    fxy_val = fxy.subs([(x, x0), (y, y0)])
    fyy_val = fyy.subs([(x, x0), (y, y0)])
    print(f_val+fx_val*(x-x0)+fy_val*(y-y0)+(1/2)*fxx_val*(x-x0)**2+fxy_val*(x-x0)*(y-y0)+\
    (1/2)*fyy_val*(y-y0)**2)
# Print the result
    return
second_ty(f)

-------------------Mid2----------------------
# Change this part for Black_scholes Calculation
#------------------------------------------------------------
S = 90
K = 100
T = 1
r = 0.1
sigma = 0.3
#-------------------------------------------------------------

def BS_CALL(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S * N(d1) - K * np.exp(-r * T) * N(d2)

def BS_PUT(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return K * np.exp(-r * T) * N(-d2) - S * N(-d1)

def delta_call(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    dt = N(d1)
    return dt

def vega(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    ve = S * n(d1) * np.sqrt(T)
    return ve

def theta_call(S,K,T,r,sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    the = -((S * n(d1) * sigma)/(2 * np.sqrt(T))) - r * K * np.exp(-r * T) * N(d2)
    return the

def theta_put(S,K,T,r,sigma):
    d1 = (np.log(S / K) + (r + sigma ** 2 / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    the = -((S * n(d1) * sigma) / (2 * np.sqrt(T))) + r * K * np.exp(-r * T) * N(-d2)
    return the

print("Value of Call:",BS_CALL(S,K,T,r,sigma))
print("Value of Put:",BS_PUT(S,K,T,r,sigma))
print("Call Delta", delta_call(S,K,T,r,sigma))
print("Put Delta", delta_call(S,K,T,r,sigma)-1)
print("Vega", vega(S,K,T,r,sigma))
print("Call Theta:",theta_call(S,K,T,r,sigma))
print("Put Theta:",theta_put(S,K,T,r,sigma))

#---------------MODE Calculation-----------
# Change variable here!
# miu in float!!!
s = 100
miu = 0.13
sigma = 0.20
T = 1
#-------------------------------------------
def MODE(s,miu,sigma,T):
    S = s * np.exp((miu - (3/2) * sigma**2 )* T)
    return S

def Prob_mode(s,miu,sigma,T):
    prob = (1 / (np.sqrt(2*np.pi))) * (1/s) * (np.exp((sigma**2-miu)*T)/(sigma * np.sqrt(T)))
    return prob

print("MODE Value:", MODE(s,miu,sigma,T))
print("Prob of MODE:",Prob_mode(s,miu,sigma,T))
