
#Practice Julia advanced lessons


#%% Lesson 9: Working with packages
# using Pkg
# Pkg.activate("Testpackage1")
#%%
using Testpackage1
Testpackage1.greet()

Testpackage1.myspecialfunction(42)

# using Pkg
# Pkg.instantiate()



#%% Lesson 10: Plotting

# using Pkg
# Pkg.add("Plots")

using Plots

plotly()
# Pkg.add("PlotlyBase")
using Plots
using PlotlyBase
x= 1:0.01:10*pi #array from 1 to 10pi using stepsize 0.01
y = sin.(x)
plot(x,y,label="sinx")
plot!(xlab="x", ylab="f(x)")

y2=sin.(x).^2 #double dot to treat individually
plot!(x,y2,label="sin(x)^2", color=:red, line=:dash)

xaxis!(:log10)
plot!(legend=:bottomleft)
#savefig("img1c.png")
#! is the append function and thus alters the current plot

#plotly works for nice interactive plots
plotly()
x=1:0.1:3*pi
y=1:0.1:3*pi

xx = reshape([xi for xi in x for yj in y], length(y), length(x))
yy = reshape([yj for xi in x for yj in y], length(y), length(x))
zz = sin.(xx).*cos.(yy)

# plot3d(xx,yy,zz,label=:none, st=:surface)
# plot!(xlab="x", ylab="y", zlab="sin(x)*cos(y)")
#
# savefig("img2")

# Pkg.add("ORCA")
using ORCA
savefig("img2.png")

gr()

plot(x,y)



#%% Lesson 11: Numerical integration

# using Pkg
# Pkg.add("QuadGK")

using QuadGK
func1(x)=exp(-x^2)
res, err = quadgk(func1,-Inf, Inf)
#Quadgk integrates
abs(res-sqrt(pi))/sqrt(pi)
res, err = quadgk(func1,-Inf, Inf, rtol=1e-15)
abs(res-sqrt(pi))/sqrt(pi)

#the latter is relative error

func2(x,y,z) = x+y^3+sin(z)

x=5
z=3
arg(y)=func2(x,y,z)
quadgk(arg,1,3)

#%% Lesson 12: Units of measurement
# using Pkg
# Pkg.add("Unitful")
using Unitful

one_meter=1*u"m"

b = uconvert(u"km", one_meter)

one_meter

#%% Lesson 13: Interacting with python
# using Pkg
# Pkg.add("PyCall")
using PyCall

#install a python library
math =pyimport("math")
math.sin(3)

using Pkg
Pkg.add("PyPlot")
using PyPlot

x=1:0.1:2*pi
y =sin.(x)
plot(x,y,label="sin(x)")
#%% Lesson 14: Data storage
#Struct for dataframe

#%% Lesson 15:  Parallel computing
