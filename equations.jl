# Newton method
# f is function, x is initila value and err is the error

function Newton(f, x, err)
    # err = 1e-4   # default err is 1e-6
    fx = f(x)
    df = (fx - f(x-err))/err
    xNew = x - fx/df
    return abs(xNew - x) > err ? Newton(f, xNew, err) : xNew
end

