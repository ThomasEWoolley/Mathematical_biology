function x = cobweb(fcnName, domain, initVal, nStart, nStop,c)

% Compute the iterates.
f    = inline(vectorize(fcnName));
x    = zeros(length(initVal),nStop);
x    = initVal;
for n = 1:nStop-1
    x(n+1) = f(x(n));
end


xMin = domain(1); xMax = domain(2);
xGrid = linspace(xMin, xMax);  
plot(xGrid,f(xGrid),'r-','linewidth',2);            % Plot the user-specified function.
hold on
plot(xGrid, xGrid,'k-','linewidth',2);        % Plot the diagonal.
axis([xMin xMax xMin xMax]);    % Scale the picture, and 
               % ensure the plot is square.

plot(initVal*[1 1],[0 x(2)],c)   % Draw vertical part of cobweb.
for n = nStart+1:nStop
        if (n > 1) 
            xplot = [x(n-1) x(n)];
            yplot = [x(n) x(n)];
            plot(xplot,yplot,c)  % Draw horizontal part of cobweb.
        end
        if (n < nStop)
            xplot = [x(n) x(n)];
            yplot = [x(n) x(n+1)];
            plot(xplot,yplot,c)   % Draw vertical part of cobweb.
        end
end 
axis('square');  
