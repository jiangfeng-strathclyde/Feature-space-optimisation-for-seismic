
function [x]=erfenfa_new_use2(f,bot,top,err,init)

if f(bot)*f(top)<0


i=ceil((log(top-bot)- log(err))/log(2))-1; 
while abs(top-bot)>err
    x=(bot+top)/2;
    fx=f(x);
    if fx==0
        bot=x; top=x;
    elseif f(bot)*fx<0
        top=x;
    else
        bot=x;
    end
end

end

if f(bot)*f(top)>=0
    x=init;
    fx=0;
    i=1;
end


end
