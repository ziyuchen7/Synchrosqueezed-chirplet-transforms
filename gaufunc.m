function g = gaufunc(x,y,sig)
g = exp(-(sum(x.^2,2)+sum(y.^2,2)'-2*x*y')/(2*sig^2));%/(2*pi*sig^2);
end