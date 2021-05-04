# ---- Lesson 1: Variables & types ------#

#%%
my_name= "Magali"
my_favourite_number = 8
my_favourtie_pie = 3.1415
print(my_name)

a=2
b=3
som = a+b
diff = a-b
prod = a*b
quot = b/a
power =a^3
modulus = b%a

print(som,diff, prod, quot, power, modulus)


#conversion of data type
b = convert(Float64,a)
a
b

#%%

#------ Lesson 2: Functions --------
function plus_two(x)
    return x+2
end

plus_two(2)


plus_three(x) = x+3
plus_four = x -> x+4

using Pkg
Pkg.add("QuadGK")
using QuadGK

f(x,y,z) = (x^2 +2y)*z

#------ Data structures -----
a=[1,2,3,4,5]
b=[1.2,3,4,5]
c=["Hello", "Its me", "Magali"]

c[3]

append!(a, 6)
#Arrays must contain elements of the same type

typeof(a)
typeof(b)
typeof(c)

matrix = [1 2 3; 4 5 6]
#for matrix thus no commas in between elements, but only spaces

zeroes = zeros(2,3,4)
print(zeroes)

table= zeros(2,3,4)
for k in 1:4
    for j in 1:3
        for i in 1:2
            table[i,j,k]=i*j*k
        end
    end
end
#iterate first over columns and then over rows
print(table)

#slices
a = [1,2,3,4,5,6]
b = a[2:5]


mat = reshape([i for i in 1:16],4,4)
print([i for i in 1:16])
mat2 = mat[2:3,2:3]

#nested comprehension
[i+j for i in 1:10 for j in 1:5]

a= [1,2,3]
b=copy(a)

b[2]=42
print(a)
print(b)

#tuples
tuple = (1,2,3)
a,b,c = tuple
a

function return_multiple()
    return 42,43,44
end

a,b,c = return_multiple()
a

print("$a $b $c")

#splatting
function splat_me(a,b,c)
    return a*b*c
end

tuple = (1,2,3)
splat_me(tuple...)
#... causes unpacking

#named tuples
namedtuple = (a=1, b="hello")
namedtuple

namedtuple[:a]

#assign names to tuple
namedtuple2 = NamedTuple{(:a, :b)}((2, "hello2"))
namedtuple2[:b]

#dictionaries
person1=Dict("Name" => "Aurelio", "Phone" =>"123456789", "Shoe size"=>"40" )
person2=Dict("Name" => "Elena", "Phone" =>"123456789", "Shoe size"=>"36" )

addressBook = Dict("Aurelio"=> person1, "Elena"=> person2)

person3=Dict("Name" => "Vittorio", "Phone" =>"123456789", "Shoe size"=>"42" )

addressBook["vittorio"]=person3


#----- Lesson 4: Control flow ------

#absolute values
function absolute(x)
    if x>= 0
        return x
    else
        return -x
    end
end
absolute(-2)

#always write an "end" after an if statement

x=42

if x<1
    print("$x<1")
elseif x<3
    print("$x<3")
elseif x<100
    print("$x<100")
else
    print("$x is really big!")
end

#by using $x we tell julia that is must substitue the value of x
#%%
name1= "traveler"
name2 ="Magali"
print("Welcome $name1, this is $name2")

#%%
#for loops
for i in 1:10
    println(i^2)
end

for i =1:10
    println(i^2)
end

persons = ["Alice", "Bob", "Carla", "Daniel"]

for person in persons
    println("Hello $person, welcome to Julia")
end

#break of the forloop
for i in 1:100
    if i>10
        break
    else
         println(i^2)
    end
end

#continue will forcefully skip the current iteration
for i in 1:30
    if i%3 ==0
        continue
    else
        println(i)
    end
end

#% deelt door maar rond af

for i in 1:30
    if i==3
        continue
    else
        println(i)
    end
end

#While loop needs to continue until a certain condition is met

function whiletest()
    i=0
    while(i<30)
        println(i)
        i+=1
    end
 end

print(whiletest())

#enumerate
x  = ["a", "b", "c"]
for sign in enumerate(x)
    println(sign)
end

# returns value i, x(i)
#manual option below
enum_array =[(1,"a"), (2,"b"), (3,"c")]
for i in 1:length(x)
    println(enum_array[i])
end

#operating and storing with array values
my_array = collect(1:10)
my_array2 = zeros(10)
my_array3=zeros(10)
for (i, element) in enumerate(my_array)
    my_array2[i] = element^2
end
for i in my_array
    print(i), print(my_array[i])
    my_array3[i] = my_array[i]^2
end
print(my_array2)
print(my_array3)

#same results

#----Lesson 5: Broadcasting -----
#error operations
a=[1,2,3]
b=[4,5,6]
#the above are column vectors
#a*b
c=[4 5 6]
#row vector
a*c

c*a

d = reshape([1,2,3,4,5,6,7,8,9], 3, 3)
#make a matrix out of this

#broadcasting: operate element by element
a .*c

c .*a

a .*d

sin.(a) #this does not work without broadcasting, as it does not treat items individually

# ----- Lesson 6: Variable scope -----

#altering function scope
function example1()
    z = 42
    return
end

function example2()
    global z = 42
    return
end
example2()
z

function example3()
    z = 42
    return z
end

z = example3()

z

#let construct for new local scope
a = let
    i=3
    i+=5
    i
end

a

b = let i=5
    i+=42
    i
end

b

#only in this scope i is 5, afterwards nothing

d = begin
    i=42
    i+=1
    i
end

d
i

#begin blocks do no introduce a local scope
const C = 299792458

C = 3000000000

#C = 2.9998 *1e8
#cannot change the type of C, easier to stick to this value

# ------- Lesson 7: Modules ------
using SpecialFunctions

gamma(3)

sinint(5) #sine integral
#Only using particular functions
#using SpecialFucntions: gamma, sinint

function delta(x)
    println("I am another 'gamma' function")
    return x^2
end

delta(3)

module MyModule
    export func2
    a=42
    function func1(x)
        return x^2
    end
    function func2(x)
        return func1(x)+a
    end
end
#%%
using .MyModule
#. before because it is not an official package
func2(3)

MyModule.func1(3)

#include codes from other files

module MyBigModule
        include("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Practice Julia/bigmodule1.jl")
        include("/Users/magali/Documents/1. Master/1.4 Thesis/02 Execution/01 Model Sarah/Practice Julia/bigmodule2.jl")
    export func2big
end

using .MyBigModule
func2big(3)


#------- Lesson 8: Types --------
abstract type Person
end

abstract type Musician<: Person
end

#musician valt onder persoon

#mutable type and immutible type
mutable struct Rockstar<:Musician
    name::String
    instrument::String
    bandname::String
    headbandcolor::String
    instrumentsplayed::Int
end

struct Classicmusician<:Musician
    name::String
    instrument::String
end

mutable struct Physicist <:Person
    name::String
    sleephours::Float64
    favouritelanguage::String
end

aure = Physicist("Aureliuo", 6, "julia")

aure.sleephours

#works like a dataframe

aure_musician = Classicmusician("Aurelio", "Violin")

ricky = Rockstar("Riccardo", "Voice", "Black Lotus", "red", 2)

ricky.headbandcolor

function introduceme(person::Person)
    println("Hello, my name is $(person.name).")
end

introduceme(aure)

function introduceme2(person::Musician)
    println("Hello, my name is $(person.name) and I play $(person.instrument).")
end

introduceme2(aure_musician)

function introduceme3(person::Rockstar)
    if person.instrument =="Voice"
        println("Hello, my name is $(person.name) and I sing")
    else
        println("Hello my name is $(person.name) and I play $(person.instrument).")
    end

    println("My band name is $(person.bandname) and my favourite colour is $(person.headbandcolor)")
end

#:: indicates that Rockstar is a subtype of person, needs to be afforementions,
#this is placed lower in the type tree

introduceme3(ricky)

#type constructor
mutable struct mydata
    x::Float64
    x2::Float64
    y::Float64
    z::Float64
    function mydata(x::Float64, y::Float64)
        x2=x^2
        z=sin(x2+y)
        new(x,x2,y,z)
    end
end

mydata(2.0,3.0)

mutable struct mydata2{T<:Real}
    x::T
    x2::T
    y::T
    z::Float64
    function mydata2{T}(x::T, y::T) where {T<:Real}
        x2=x^2
        z=sin(x2+y)
        new(x,x2,y,z)
    end
end

#assign type later, for in real dimension

mydata2{Float64}(2.0,3.0)

mydata2{Int}(2,3)


#storing data that needs to be shared between funcitons inside a module

module Testmoduletypes
    export Circle, computeperimeter, compurearea, printcircleeq

    mutable struct Circle{T<:Real}
        radius::T
        perimeter::Float64
        area::Float64

        function Circle{T}(Radius::T) where T<:Real
            new(Radius,-1.0,-1.0)
        end
    end

    @doc raw"""
            computearea(circle::Circle)

     Compute the area of `circle` and store the value.
     """
    function computearea(circle::Circle)
        circle.area = pi*circle.radius^2
        return circle.area
    end

    function printcircleeq(xc::Real, yc::Real, circle::Circle)
        println("(x-$xc)^2 + (y-$yc)^2 = $(circle.radius^2)")
        return
    end
end

using .Testmoduletypes
circle1 = Circle{Float64}(5.0)
circle1.area
computearea(circle1)
