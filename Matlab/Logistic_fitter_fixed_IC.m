ccc

Data1=[1790 	3.929
    1800 	5.308
    1810 	7.240
    1820 	9.638
    1830 	12.866
    1840 	17.069];

Data2=[1850 	23.192
    1860 	31.443
    1870 	38.558
    1880 	50.156
    1890 	62.948
    1900 	75.996
    1910 	91.972
    1920 	105.711
    1930 	122.775
    1940 	131.669];

Data3=[1950 	150.697
    1960 	179.323
    1970 	203.185
    1980 	226.546
    1990 	248.710];

[x,a,b,c,d,e,f] = lsqcurvefit(@logistic, [150 1e-2], [Data1(:,1);Data2(:,1)], [Data1(:,2);Data2(:,2)]);
logistic(x,[Data1(:,1);Data2(:,1);Data3(:,1)])

hold on
L=linspace(1790,1990);
plot(L,logistic(x,L))
plot(Data1(:,1),Data1(:,2),'b*')
plot(Data2(:,1),Data2(:,2),'k*')
plot(Data3(:,1),Data3(:,2),'r*')
axis tight
xlabel('Year','interpreter','latex')
ylabel('US Population in millions','interpreter','latex')
export_fig('../Pictures/US_logistic_population.png','-r300')
function F = logistic(x,xdata)
F=x(1)./(1+((x(1)-3.929)/3.929)*exp(-x(2)*(xdata-1790)));
end
