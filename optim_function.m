function f = optim_function(x,minx,maxx)
for i=1:length(x)
    if x(i)<minx
       %f(i)= abs(x(i)-minx); 
       f(i)= sqrt((x(i)-minx).^2+3^2)-3;
    elseif x(i)>maxx
       %f(i) = abs(x(i)-maxx);
       f(i)= sqrt((x(i)-maxx).^2+3^2)-3;
    else
       f(i) = 0;
    end
end