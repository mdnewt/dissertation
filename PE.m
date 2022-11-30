%percent error function
function val = PE(obs,exp)
    val = 100 .* (obs - exp) ./ exp;
end