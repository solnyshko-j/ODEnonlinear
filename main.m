clear all, close all
Position1 = [108 540 560 420];
Position2 = [680 540 560 420];
Position3 = [1250 540 560 420];
figure('Position',Position1)
figure('Position',Position2)
figure('Position',Position3)
figure(2)
format
hold on
% y" = 2*xR^2*(1-theta+theta*pi/2*cos(pi/2*ksi))^2 .* 
%      (y'./(xR*(1-theta+theta*pi/2*cos(pi/2*ksi)))+ 1).*
%      (y + xR*((1-theta)*ksi +
%      theta*sin(pi/2*ksi)))-pi^2/4*y'*sin(pi/2*ksi)/(1-theta+theta*pi/2*cos(pi/2*ksi))
% x: [0, 1,56]
% y(0) = CL, y'(ksi_l)/(1-theta+theta*pi/2*cos(pi/2*ksi))/xR = CR
% y_a = tg(ksi(ksi)) - ksi(ksi)
xL    = 0;
xR    = 1.56;
ksiL = 0.0;
ksiR = 0.5;
YL    = 0;
YR    = -1.53276613369;
CL = 0;
CR = tan(xR)^2;
Eps   = 1e-4;
alpha = 0.0;
minmesh = 6;
maxmesh = 25;
mesh = minmesh;
maxiter = 30;
T = 10;
Pause = 1;

tic
M = 2 ^ (mesh + 2);
h = 1/M;
ksi = zeros(1, M + 1);
for m = 1:(M + 1)
  ksi(m) = (m - 1)*h;
end
theta = 0.9;
y_noise = sin(8 .* xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) + ...
    cos(8 .* xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi)));
yn = (1 - alpha) .* (tan(xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) -...
    xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) + alpha .* y_noise;
y_a = tan(xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) - xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi));
plot(ksi, yn, 'g-', ksi, y_a , 'r-', ksi, y_noise,'k--', 'LineWidth', 2),
title(sprintf('y(x): \\alpha=%0.2f',alpha)),drawnow
LEG2 = {'y_{ini}','y_a','y_{noise}'};
legend(LEG2,'Location','NorthWest')
xlabel('x'); ylabel('y')
fprintf(1,' alpha = %d theta = %d Eps = %d\n xL = %d xR = %d\n maxiter = %i maxmesh = %i\n',...
    alpha, theta,Eps, xL, xR, maxiter, maxmesh);
fprintf('Newtonian iterations have begun.\n');

for iter = 1:maxiter % Цикл по Ньютоновским итерациям
  
  figure(1)
  hold off
  plot([ksiL,ksiR],[0,0],'w.')
  tit1 = sprintf('Correction: iter = %i mesh = ',iter);
  title(tit1)
  pause(Pause)
  hold on
  LEG1 = {};
  
  for mesh = minmesh:maxmesh % Внутренний цикл
    
    M = 2 ^ (mesh + 2);
    h = 1/M;
    ksi = zeros(1, M + 1);
    for m = 1:(M + 1)
      ksi(m) = (m - 1)*h;
    end
  
 
    if (length(yn) > length(ksi)) % Интерполяция для y
      ynext = yn;
      k = 2 ^ (log2(length(yn) - 1) - log2(M));
      y = ksi;
      for i = 0:M
        y(i + 1) = yn(k * i + 1);
      end
    elseif (length(yn) == length(ksi))
      ksi_2 = (ksi(2:end) + ksi(1:end - 1)) / 2;
      
      P = ones(M,T);
      y_2 = zeros(M,1);
 
      for i = 1:M
        for j = 1:T
          for k = 1:T
            if k ~= j
              if(i <= T/2 - 1)
                P(i,j) =  P(i,j)*(ksi_2(i) - ksi(k)) / (ksi(j) - ksi(k));
                temp = j;
              elseif(i >= M - T/2 + 1)
                P(i,j) = P(i,j)*(ksi_2(i) - ksi(M - T + k))/(ksi(M - T + j) - ksi(M - T + k));
                temp = M - T + j;
              else
                P(i,j) = P(i,j)*(ksi_2(i) - ksi(i - T/2 + k))/(ksi(i - T/2 + j) - ksi(i - T/2 + k));
                temp = i - T/2 + j;
              end
            end
          end
          y_2(i) = y_2(i)+ yn(temp) * P(i,j);
          
        end 
      end
      
     ynext(1) = yn(1);
     for i = 1:M
       ynext(2*i) = y_2(i); 
       ynext(2*i + 1) = yn(i + 1);
     end
       y = yn;
    end % Конец интерполяции для y
    
    dydksi = 0*ksi;
    dydksi2 = 0*ksi;
    
%   Приближение первой производной
    dydksi(1) = (-3/2*y(1)+2*y(2)-1/2*y(3))/h;
     
    for m = 2:M
      dydksi(m) = (y(m+1) - y(m-1))/2/h;
    end
    dydksi(M+1) = (1/2*y(M-1)-2*y(M)+3/2*y(M+1))/h;
%   Приближение первой производной
    dydksi2(1) = (2*y(1)-5*y(2) + 4*y(3)-y(4))/h^2;
     
    for m = 2:M
      dydksi2(m) = (y(m-1) - 2*y(m)+y(m+1))/h^2;
    end
    dydksi2(M+1) = (-y(M-2)+4*y(M-1) - 5*y(M)+2*y(M+1))/h^2;
    
    q = 2*xR*(1-theta+theta*pi/2*cos(pi/2*ksi)).*...
        (y + xR*((1-theta).*ksi + theta .*sin(pi/2.*ksi))) +...
        1/4*pi^2*sin(pi/2*ksi)/(1-theta+theta*pi/2*cos(pi/2*ksi));
        
    p = 2*xR^2*(1-theta+theta*pi/2*cos(pi/2*ksi)).^2.*...
        (dydksi./(1-theta+theta*pi/2*cos(pi/2*ksi))/xR + 1);
    
    phi = 2*xR^2*(1-theta+theta*pi/2*cos(pi/2*ksi)).^2 .* ... 
      (dydksi./(xR*(1-theta+theta*pi/2*cos(pi/2*ksi)))+ 1).* ...
      (y + xR*((1-theta)*ksi + theta*sin(pi/2*ksi)))...
      -pi^2/4*dydksi.*sin(pi/2*ksi)./(1-theta+theta*pi/2*cos(pi/2*ksi))...
      - dydksi2;
   
    phi(M+1) = CR*xR*(1-theta+theta*pi/2*cos(pi/2*ksiR))-dydksi(M+1);
    clear dydksi;
    clear dydksi2;
    
    a = zeros(M+1, 1);
    b = zeros(M+1, 1);
    c = zeros(M+1, 1);
    b(1) = 1;
    phi(1) = YL - y(1);
    for m=2:M
      coefs1_1 = [-1/(2*h); 0; 1/(2*h)];
      coefs1_2 = [1/h^2; -2/h^2; 1/h^2];
      c(m) = coefs1_2(3) - q(m)*coefs1_1(3);
      b(m) = coefs1_2(2) - q(m)*coefs1_1(2)- p(m);
      a(m) = coefs1_2(1) - q(m)*coefs1_1(1);
    end
     
    coefs = [1/2; -2; 3/2]/h;
    a(M+1) = coefs(2) - b(M)/a(M)*coefs(1);
    b(M+1) = coefs(3)- c(M)/a(M)*coefs(1);
    phi(M+1) = phi(M+1) - phi(M)/a(M)*coefs(1);
    v_temp = Progonka(a,b,c, phi);  
    clear p;
    clear q;
    clear phi;
    clear a;
    clear b;
    clear c;
    delta(1) = 0;
    ratio(1) = 0;
    ratio(2) = 0;
    figure(1)
    c1 = 'rbkgmc';
    nc1 = mod(mesh - 1,6) + 1;
    t_ = round(log2(M))-1;
    plot(ksi(1:t_:end), v_temp(1:t_:end), c1(nc1))
    tit1 = sprintf('Correction: iter = %i mesh = %i M = %i \\delta_h = ', iter, ...
                    mesh, M);
    if mesh - minmesh >= 1
      tit1 = [tit1, sprintf('%0.2e R_h = ', delta(mesh - minmesh))];
      if mesh - minmesh >= 2
        tit1 = [tit1, sprintf('%0.3f', ratio(mesh - minmesh - 1))];
      end
    end
    title(tit1)
    LEG1(mesh - minmesh + 1) = {['M = ',int2str(M)]};
    legend(LEG1);
    xlabel('ksi'); ylabel('v')
    pause(Pause)
    
    if mesh >= (minmesh + 1)
      v_temp_avg = zeros(1,M / 2 + 1);
      v_temp_avg(1) = v_temp(1);
      for i = 2:(M / 2)
        v_temp_avg(i) = v_temp(2 * i - 1);
      end
      v_temp_avg(M / 2 + 1) = v_temp(M + 1);
      delta(mesh - minmesh + 1) = max(abs(v - v_temp_avg));
      if mesh >= (minmesh + 2)
        ratio(mesh - (minmesh + 1) + 2) = delta(mesh - (minmesh + 1) + 1) / delta(mesh - minmesh + 1);
      end
      if delta(mesh - minmesh + 1) < Eps / 2 % Выход из внутреннего цикла
        v = v_temp;
        fprintf(1,'iter = %2i mesh = %2i M = %7i d_h = %0.2e Eps/2 = %0.2e order = %0.3f\n',...
        iter, mesh, M, delta(mesh - minmesh + 1), Eps/2, log2(ratio(mesh - minmesh + 1)));
        break;
      end
    end
    
    fprintf(1,'iter = %2i mesh = %2i M = %7i d_h = %0.2e Eps/2 = %0.2e order = %0.3f\n',...
      iter, mesh, M, delta(mesh - minmesh + 1), Eps/2, log2(ratio(mesh - minmesh + 1)));
    
    v = v_temp;
    yn = ynext;
  end % Конец внутреннего цикла
  
  if length(v) < length(yn) % Интерполяция для v
    k = (log2(length(yn) - 1) - log2(length(v) - 1));
    for i = 1:k
      M = length(v) - 1;
      h = 1 / M;
      ksi = zeros(1, M + 1);
      for m = 1:(M + 1)
        ksi(m) = xL + (m - 1) * h;
      end

      ksi_2 = (ksi(2:end) + ksi(1:end - 1)) / 2;
      P = ones(M,T);
      v_2 = zeros(M,1);
      for i = 1:M
        for j = 1:T
          for k = 1:T
            if k ~= j
              if(i <= T/2 - 1)
                P(i,j) = P(i,j)*(ksi_2(i) - ksi(k))/(ksi(j) - ksi(k));
                temp = j;
              elseif(i >= M - T/2 + 1)
                P(i,j) = P(i,j)*(ksi_2(i) - ksi(M - T + k))/(ksi(M - T + j) - ksi(M - T + k));
                temp = M - T + j;
              else
                P(i,j) = P(i,j)*(ksi_2(i) - ksi(i - T/2 + k))/(ksi(i - T/2 + j) - ksi(i - T/2 + k));
                temp = i - T/2 + j;
              end
            end
          end
          v_2(i) =v_2(i) + v(temp) * P(i,j);
        end 
      end
      
      vnext = zeros(1,2*M + 1);
      vnext(1) = v(1);
      for i = 1:M
        vnext(2*i) = v_2(i); 
        vnext(2*i + 1) = v(i + 1);
      end
      v = vnext;
      h = 1/(2*M);
      ksi = zeros(1, 2*M + 1);
      for m = 1:(2*M + 1)
        ksi(m) =(m - 1)*h;
      end
    end
  end % Конец интерполяции для v
  
  yn =yn+v;
  
  fprintf(1,'For Newtonian iteration %i norm of residual ||v|| = %e\n',...
        iter,max(abs(v)));
  figure(2)
  nc2 = mod(iter - 1,5) + 1;
  c2 = 'bkgmc';
  t_ = round(log2(M))+1;
  plot(ksi(1:t_:end), yn(1:t_:end), [c2(nc2),'-'])
  LEG2(iter + 3) = {['i = ',int2str(iter)]};
  legend(LEG2,'Location','NorthWest')
  title(sprintf('y(x): \\alpha = %0.2f iter = %i ||v|| = %e',alpha,iter,max(abs(v))))
  pause(Pause)
  
  y_a_tmp = tan(xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) - xR.*((1-theta).*ksi + theta .*sin(pi/2 .* ksi));
  fprintf(1, 'max(|ya - yn|)=%f\n', max(abs(y_a_tmp - yn)));
  
  figure(3);
  plot(iter, log10(max(abs(v))),'b.', 'MarkerSize', 20)
  title(sprintf('Convergence: \\alpha = %0.2f \\epsilon/2 = %0.2e', alpha, Eps/2))
  xlabel('iteration'); ylabel('log_{10}(||v||)')
  hold on
  plot(0.1,log10(Eps/2),'w.',[0,maxiter],[log10(Eps/2),log10(Eps/2)],'k:')
  pause(Pause)

  if (max(abs(v)) < Eps / 2) % Выход из ньютоновского цикла
    
    break;
  end

end % Конец Ньютоновского цикла
fprintf(1,'%i iterations were done.\n',iter); 
if max(abs(v)) < Eps/2
    fprintf(1,'The Newton process has converged.\n');
else
    fprintf(1,'The Newton process has not converged.\n');
end
if max(abs(yn - (tan((1-theta).*ksi + theta .*sin(pi/2 .* ksi))) - (1-theta).*ksi + theta .*sin(pi/2 .* ksi))) < Eps
    fprintf(1,'The solution was found.\n');
else
    fprintf(1,'The solution was not found.\n');
end
fprintf(1,'alpha = %d eps/2 = %d\n', alpha, Eps/2);
toc