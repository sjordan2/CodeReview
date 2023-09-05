% Uses Newton's Method to approximate roots for a function given its first
% derivative as input

function newtonsMethod_forReview()
    func_handle = @(x) (x^6 - x - 1); % f(x)
    dev_handle = @(x) (6*x^5 - 1); % f'(x)
    dev2_handle = @(x) (30*x^4); % f''(x)
     
    x_zero = 2.0;
    tolerance = 1.0e-3;
    iteration_max = 50;
     
    [x, iter] = newton_calc(func_handle, dev_handle, x_zero, tolerance, iteration_max);
    disp(x);
    desc = ["The approximated root is",x] ;
    disp(desc)
    disp(iter);
    desc2 = ["The iteration needed is ", iter];
    disp(desc2)
end

 
function [x, it] = newton_calc (fhandle, jhandle, xini, tol, itmax)
 
    alpha = 1.1347241384015195; % Known root of function
    table = array2table(zeros(0,7), 'VariableNames',{'n', 'x_n', 'f(x_n)',...
        'alpha - x_n', 'x_n+1 - x_n', 'R_1', 'R_2'});
    R1 = 0.0; % R1 and R2 are undefined for the zeroth iteration, 
                    % so set them equal to 0
    R2 = 0.0;
    x = xini; % Use x as a current root counter
    err = 1.0 + tol; % Set error to sufficiently high number to start off with
    iter = 0; % Set the current iteration equal to 0
     
    % Continue iterating until max iterations or tolerance has been reached.
    while((err > tol) && (iter < itmax)) 
        f = fhandle(x);
        F = jhandle(x);
        dx = - f / F;
    
        if(iter > 0)
            R1 = ratio_one(alpha, x, xold);
            R2 = ratio_two(alpha, x, xold);
        end
        newRow = {iter, x, f, alpha-x, dx, R1, R2};
        table = [table; newRow];
        xold = x;
        x = x + dx;
        err = abs(dx);
        iter = iter + 1;
    end
    
    f = fhandle(x);
    F = jhandle(x);
    dx = - f / F;
    newRow = {iter, x, f, alpha-x, dx, R1, R2};
    table = [table; newRow];
    disp(table);
    it = iter;
end
 
function R1 = ratio_one(alpha, x, xold)
    R1 = (alpha - x) / (alpha - xold);
end
 
function R2 = ratio_two(alpha, x, xold)
    R2 = (alpha - x) / (alpha - xold)^2;
end
