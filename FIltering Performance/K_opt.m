function output = K_opt(c_nx,n_x,norm_x)

    %optimal Epanechnikov kernel when our prior weights values are equal.
    if norm_x<1
        output = (1-norm_x^2)*(n_x+2)/(2*c_nx);
    else
        output = 0;
    end
end