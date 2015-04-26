% Plotting the log error

err = fliplr(log([2.3102 0.7667 0.2162 0.0562 0.0142 0.0036]));
h = fliplr(log([2.25 1.125 0.5625 0.281250 0.140625 0.070312]));
p = polyfit(h, err, 1); % Fit to a line
P = polyval(p, h);
disp(p(1))

set(gcf, 'color', 'w')
plot(h, err, h, P)
title('Calculating Convergence Rate')
xlabel('log(h)')
ylabel('log(err)')
legend('convergence', 'polyfit')


