% Önerilen analitik çözümün grafiği
t = linspace(0, 1, 100);  % [0, 1] aralığında 100 nokta
y = tan(pi/4 * (1 - t));

% Grafiği çizdirme
figure;
plot(t, y, 'LineWidth', 2);
xlabel('t');
ylabel('y(t)');
title('Önerilen analitik çözüm: y''(t) + y(t)^2 = 0');
grid on;