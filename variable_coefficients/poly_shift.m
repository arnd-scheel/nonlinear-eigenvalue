%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=poly_shift(A,r,c,p,x);

% This is the main routine that computes the shifted polynomial b from a
% running degree of polynomial, decreasing as we take more derivatives
% coefficient of polynomial are rxc matrices, A is r*c*(p+1) matrix with coefficients stored in blocks concatenated horizontally, starting with lowest order
pn=p;

% now create the new polynomial B by successfully computing derivatives and evaluating at the shift location
for der=0:p % loop over p derivatives
    % Horner's scheme for evaluating
    A_val=A(1:r,(1:c)+c*(pn));
    for j=pn:-1:1
        A_val=A_val*x+A(1:r,(1:c)+c*(j-1));
    end
    B(1:r,(1:c)+c*(der))=A_val/factorial(der); %compute the polynomial coefficient from the value at shifted point
    % now take a derivative
    for j=0:pn
        A(1:r,(1:c)+c*(j))=A(1:r,(1:c)+c*(j))*j;
    end
    A=A(1:r,c+1:end); % drop lowest (zero) coefficient
    pn=pn-1;    % reduce degree after taking derivative
end

return
