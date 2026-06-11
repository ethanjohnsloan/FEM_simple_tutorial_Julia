
using Plots
pythonplot()
using QuadGK

function hat(x, a, b)
    c = (a+b)/2

    if x < a || x > b
        return 0.0
    elseif a <= x < c
        return (x-a)/(c-a)
    else
        return (b-x)/(b-c)
    end
end
