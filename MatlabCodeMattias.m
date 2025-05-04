%% stirling series
%coeffficients
n = 0:100;
b = zeros(1,length(n));
b(1) = 1;
b(2) = 1;
b(3) = b(2)/3;
for i = 4:length(n)
    som = 0;
    for j = 2:i-2
        som = som + j*b(j+1)*b(i-j);
    end
    b(i)= (b(i-1)- som)/(i+1);
end
a = zeros(1,(length(n)-1)/2);
dubbelFact = 1;
for i = 1:(length(n)-1)/2
    dubbelFact = dubbelFact*(2*i - 1);
    a(i) = b(2*i)*dubbelFact;
end

%approximation for a specific x-value.
x = [1,2,3,4,30]; % not x = 0.
approx = zeros(length(x),length(a));
exact = factorial(x); % or equally gamma(x+1)
for k = 1:length(x)
    a_Div_n_k = zeros(1, length(a));
    a_Div_n_k(1) = a(1);
    for i = 2:length(a)
        a_Div_n_k(i) = a_Div_n_k(i-1) + a(i)/x(k)^(i-1);
    end
    approx(k,:) = (x(k)^(x(k))) * (exp(-x(k))) * (sqrt(2*pi*x(k))) * (a_Div_n_k);
end
OptTruncExact = zeros(2,length(x));
for i = 1:length(x)
    [OptTruncExact(1,i),I] = min(abs(exact(i)-approx(i,:)));
    OptTruncExact(2,i) = I;
end

N_Ster = x;
OptTrunc = zeros(2,length(x));
for i = 1:length(x)
    OptTrunc(1,i) = abs(exact(i)-approx(i,N_Ster(i)));
    OptTrunc(2,i) = N_Ster(i);
end

figure(1) %approximation in function of # terms
plot(1: length(a), abs(exact(1) - approx(1,:)),'.','MarkerSize',13)
hold on
plot(1: length(a), abs(exact(2) -approx(2,:)),'.','MarkerSize',13)
plot(1:length(a), abs(exact(3)-approx(3,:)),'.','MarkerSize',13)
%plot(0:length(a)-1, abs(exact(4)-approx(4,:)),'.','MarkerSize',13)
plot(OptTruncExact(2,1:3), OptTruncExact(1, 1:3) ,'r+','MarkerSize',10)
plot(OptTrunc(2,1:3), OptTrunc(1, 1:3) ,'go','MarkerSize',10)

yscale('log') % 1:length(a), abs(exact(4)-approx(4,:)),'.'
legend('x = 1','x = 2','x = 3', 'Minimal Error','Optimal Truncation') %, 'Absolute error (x = 30)'
xlabel('Number of Terms in the Approximation')
ylabel('Absolute Approximation Error (log-scale)')
title('The Stirling Series for the Gamma Function \Gamma(x+1)')
ylim([10^(-5),10^(25)])
xlim([0,25])
hold off


%% exp function
function som = borelExp(xi)
som = 0;
for n = 0:100
    som = som + xi.^n./(factorial(n).*factorial(n+1));
end
end

integrated_funcPos = [];
x = 0:0.001:10;
for k = 1:length(x)
fun = @(xi) borelExp(xi).*exp(-x(k)*xi);
    integrated_funcPos(k) = integral(fun,0,100);
end
integrated_funcNeg = [];
x2 = -10:0.001:0;
for k = 1:length(x2)
fun = @(xi) borelExp(xi).*exp(-x2(k)*xi);
    integrated_funcNeg(k) = integral(fun,0,-100);
end

figure(3)
plot(x,integrated_funcPos, -10:0.001:10, exp(1./(-10:0.001:10)),x2, integrated_funcNeg)
title('The Borel resummation of exp(\frac{1}{x})')
legend('resummation \theta = 0', 'exp(\frac{1}{x})','resummation \theta = \pi')
ylim([-1,10])
xlabel('x-value')