%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=poly_shift(A,p,x);

% This is the main routine that computes the shifted polynomial b from a
% running degree of polynomial, decreasing as we take more derivatives
pn=p;

% now create the new polynomial B by successfully computing derivatives and evaluating at the shift location
for der=0:p % loop over p derivatives
    % Horner's scheme for evaluating
    A_val=A(:,:,pn+1);
    for j=pn:-1:1
        A_val=A_val*x+A(:,:,j);
    end
    B(:,:,der+1)=A_val/factorial(der); %compute the polynomial coefficient from the value at shifted point
    % now take a derivative
    for j=0:pn
        A(:,:,j+1)=A(:,:,j+1)*j;
    end
    A=A(:,:,2:end); % drop lowest (zero) coefficient
    pn=pn-1;    % reduce degree after taking derivative
end

return
