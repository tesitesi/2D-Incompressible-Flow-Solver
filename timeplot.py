import matplotlib.pyplot as plt 

beta = [0.3,0.5,0.6,0.7,0.8,0.9,1,1.5]
t = [33553,29732,28073,28135,28662,28757,30401,32977]

err = [434,523,434,453,454,423,434,530]

plt.plot(beta,t)
plt.errorbar(beta,t,yerr=err,fmt='ro',ecolor='g')

plt.xlabel('beta')
plt.ylabel('Run Time [ms]')
plt.show()