
%-------OECD high income data------%


clc;
clear all;

format longG 

%----Import data----%

[lhe]=xlsread("HEdata/HE.xlsx",'B52:U73'); % per capita health expenditure (dependent variable)
[lghe]=xlsread("HEdata/GHE.xlsx",'B52:U73'); % per capita gvt health expenditure 
[lgdp]=xlsread("HEdata/GDP.xlsx",'B52:U73'); % per capita GDP
[dr65]=xlsread("HEdata/DR65.xlsx",'B52:U73'); % Dependency rate: proportion of population aged above 65
[dr15]=xlsread("HEdata/DR15.xlsx",'B52:U73'); % Dependency rate: proportion of population aged below 15

[lheaus] = xlsread("HEdata/HE.xlsx","B14:U14");
[lhejap] = xlsread("HEdata/HE.xlsx","B10:U10");
[lhemex] = xlsread("HEdata/HE.xlsx","B2:U2");
[lhehun] = xlsread("HEdata/HE.xlsx","B30:U30");

[lgdpaus] = xlsread("HEdata/GDP.xlsx","B14:U14");
[lgdpjap] = xlsread("HEdata/GDP.xlsx","B10:U10");
[lgdpmex] = xlsread("HEdata/GDP.xlsx","B2:U2");
[lgdphun] = xlsread("HEdata/GDP.xlsx","B30:U30");

Z=zeros(5,20,22); % to store Zi
X=zeros(4,20,22); % to store Xi
Y=zeros(20,22); % to store Yi


%----Estimation----%

for i=1:22
%     lhei=zeros(1,20);
%     lghei=zeros(1,20);
%     lgdpi=zeros(1,20);
%     dr65i=zeros(1,20);
%     dr15i=zeros(1,20);

    lhei=lhe(i,:);
    lghei=lghe(i,:);
    lgdpi=lgdp(i,:);
    dr65i=dr65(i,:);
    dr15i=dr15(i,:);

    Xi=[lghei'  lgdpi'  dr65i'  dr15i'];
    Yi=lhei';
i;
    Zi=[Yi' ; Xi'];

    Z(:,:,i)=Zi;
    X(:,:,i)=Xi';
    Y(:,i)=Yi';
end

% Creating the weighted aggregate matrix Zw with wi=w=1/N where N=28

z=zeros(5,20);
for j=1:22

    z=z+Z(:,:,j);
    z=z;

end

Zw=(1/22)*z;

% Creating the projection matrix Mw
Hw=Zw'; % This matrix stores all the observed data
Mw=eye(20)-Hw*inv(Hw'*Hw)*Hw';

% Estimating the CCE estimator betai for each i=1:N

 betahat=zeros(4,1); % To store the four estimated coefficients for each i

for i=1:22
    xi=X(:,:,i)';
    yi=Y(:,i);
i;
    betaihat = inv(xi'*Mw*xi)*(xi'*Mw*yi);
    betahat =  [betahat betaihat];
    betaihat = zeros(4,1);
  
end
betahat;

% Estimating the CCEMG estimator for each four coefficients
b1 = betahat(1,:)';
b2 = betahat(2,:)';
b3 = betahat(3,:)';
b4 = betahat(4,:)';



% Plotting the kernal function for GDP coefficients
[f,xi,bw] = ksdensity(b2)
figure
plot(xi,f);
title("CCE estimation for the GDP coefficients for the high income group");
ylabel('Density');

weight = (1/22)*ones(1,23);
betaMG =[weight*b1; weight*b2; weight*b3; weight*b4];



%-------------------------------------------------------------------------------------------------------------------------------%

% Estimating the covariance matrix of CCEMG using the formula given in the paper for each iteration
    covMG=zeros(4,4);
  
    for i=1:22
     mat1i = betahat(:,i)-betaMG(:,1);
   
     covarianceMGi = covMG + (mat1i)*(mat1i)';
     covMG = covarianceMGi;
      
    end
    estcovMG = (1/21)*covMG;
    
% Estimating the covariance matrix of the CCEMG estimator under the new adjusted assumptions
   newcovMGniter = zeros(4,4);
   
   covMGnew = zeros(4,4);
   for i=1:22
     e_i = Mw*(Y(:,i)-(X(:,:,i)'*betahat(:,i)));
     x_i = Mw*X(:,:,i)';
     

     m = floor(20^(1/4));
%      m=T;
     sumi = 0;

     for t = 1:m
          h_it = e_i(t,:)'*x_i(t,:)';
          sumit = sumi + h_it*h_it';
          sumi = sumit;
   
     end
     sigmai0 = (1/20)*sumi;
     
     suminew = 0;
     
     for j=1:m
         weightedi = 0;
         for t=j+1:m
             h_it = e_i(t,:)'*x_i(t,:)';
             h_itj =  e_i(t-j,:)'*x_i(t-j,:)';
             sumitj = suminew + (h_it*h_itj'+h_itj*h_it');
            

         end
         weighteditj = (1-(j/(m+1)))*(1/20)*sumitj;
         weightedi = weightedi + weighteditj;

     end

     Siepsilon = sigmai0 + weightedi;
     

     sigmahatiMG = inv((X(:,:,i)*Mw*X(:,:,i)')/20)*Siepsilon*inv((X(:,:,i)*Mw*X(:,:,i)')/20);
     covarianceMGi2 = covMGnew + sigmahatiMG;      
     covMGnew = covarianceMGi2;
    
     
   end
   estcovMGnew =  (20/22)*estcovMG + (1/22)*covMGnew;

%--------------------------------------------------------------------------------------------------------------------------------%

%--Factor analysis using PCA--%

% u = Y -lghe'* b1(2:23) -lgdp'* b2(2:23) - dr65'*b3(2:23) - dr15'*b4(2:23);
% results=NbFactors(u,5);

betameanhigh = mean(betahat');
betasdhigh = std(betahat');
   
%----------------------------------------------------------------------------------------------------------------------------------%

%----Data Visualization----%

% ... (previous code)

% Scatter plot for lheaus against lgdpaus with regression line
figure;
subplot(1, 2, 1);
scatter(lgdpaus, lheaus, 'filled');
hold on;

% Fit and plot regression line
coeffs_aus = polyfit(lgdpaus, lheaus, 1);
fitline_aus = polyval(coeffs_aus, lgdpaus);
plot(lgdpaus, fitline_aus, 'r-', 'LineWidth', 1);

xlabel('log(GDP)');
ylabel('log(Health Care Expenditure)');
title('Australia');
grid on;
%legend('Data', 'Regression Line');

% Add annotation for slope
slope_text_aus = sprintf(' %.4f', coeffs_aus(1));
text(min(lgdpaus), max(lheaus),  ['\beta_{Aus} = ' slope_text_aus], 'Color', 'red', 'FontSize', 10, 'VerticalAlignment', 'top');

% Scatter plot for lhejap against lgdpjap with regression line
subplot(1, 2, 2);
scatter(lgdpjap, lhejap, 'filled');
hold on;

% Fit and plot regression line
coeffs_jap = polyfit(lgdpjap, lhejap, 1);
fitline_jap = polyval(coeffs_jap, lgdpjap);
plot(lgdpjap, fitline_jap, 'r-', 'LineWidth', 1);

xlabel('log(GDP)');
ylabel('log(Health Care Expenditure)');
title('Japan');
grid on;
%legend('Data', 'Regression Line');

 % Add annotation for slope
slope_text_jap = sprintf(' %.4f', coeffs_jap(1));
text(min(lgdpjap), max(lhejap),['\beta_{Jap} = ' slope_text_jap], 'Color', 'red', 'FontSize', 10, 'VerticalAlignment', 'top');


% Scatter plot for lhemex against lgdpmex with regression line
figure;
subplot(1, 2, 1);
scatter(lgdpmex, lhemex, 'filled');
hold on;

% Fit and plot regression line
coeffs_mex = polyfit(lgdpmex, lhemex, 1);
fitline_mex = polyval(coeffs_mex, lgdpmex);
plot(lgdpmex, fitline_mex, 'r-', 'LineWidth', 1);

xlabel('log(GDP)');
ylabel('log(Health Care Expenditure)');
title('Mexico');
grid on;
%legend('Data', 'Regression Line');

% Add annotation for slope
slope_text_mex = sprintf(' %.4f', coeffs_mex(1));
text(min(lgdpmex), max(lhemex),['\beta_{Mex} = ' slope_text_mex], 'Color', 'red', 'FontSize', 10, 'VerticalAlignment', 'top');

% Scatter plot for lhehun against lgdphun with regression line
subplot(1, 2, 2);
scatter(lgdphun, lhehun, 'filled');
hold on;

% Fit and plot regression line
coeffs_hun = polyfit(lgdphun, lhehun, 1);
fitline_hun = polyval(coeffs_hun, lgdphun);
plot(lgdphun, fitline_hun, 'r-', 'LineWidth', 1);

xlabel('log(GDP)');
ylabel('log(Health Care Expenditure)');
title('Hungary');
grid on;
%legend('Data', 'Regression Line');

% Add annotation for slope
slope_text_hun = sprintf(' %.4f', coeffs_hun(1));
text(min(lgdphun), max(lhehun), ['\beta_{Hun} = ' slope_text_hun], 'Color', 'red', 'FontSize', 10, 'VerticalAlignment', 'top');

% ... (any additional code or plotting adjustments)


% %--Beta (slope coefficients)--%
% 
% % Plotting the coefficients for gvt health expenditure for each i
% figure(1)
% ylim([-2,2]);
% B1=plot(1:23,b1,'o','color','red');
% title('GHE');
% xlabel('high income countries');
% ylabel('β_{GHEi}');
% ylim([-2,2]);
% hold on
% xlim([0,22]);
% plot(1:23,betaMG(1,1),'_','Color','blue');
% hold off
% 
% % Plotting the coefficients for GDP for each i
% figure(2)
% ylim([-2,2]);
% B2=plot(1:23,b2,'o','color','red');
% title('GDP');
% xlabel('high income countries');
% ylabel('β_{GDPi}');
% ylim([-2,2]);
% hold on
% xlim([0,22]);
% plot(1:23,betaMG(2,1),'_','Color','blue');
% hold off
% 
% % Plotting the coefficients for DR65
% figure(3)
% ylim([-2,2]);
% B2=plot(1:23,b3,'o','color','red');
% title('DR65');
% xlabel('high income countries');
% ylabel('β_{DR65i}');
% ylim([-2,2]);
% hold on
% xlim([0,22]);
% plot(1:23,betaMG(3,1),'_','Color','blue');
% hold off
% 
% % Plotting the coefficients for DR15
% figure(4)
% ylim([-2,2]);
% B2=plot(1:23,b4,'o','color','red');
% title('DR15');
% xlabel('high income countries');
% ylabel('β_{DR15i}');
% ylim([-2,2]);
% hold on
% xlim([0,22]);
% plot(1:23,betaMG(4,1),'_','Color','blue');
% hold off


%--Independent variables--%

% % Plotting GHE of each country against time 
% for i=1:28
% %     i
%     X2ti=X(1,:,i)';
%     V1=plot(linspace(1990,2019,20),X2ti);
%     title('GHE');
%     xlabel('year');
%   hold on
%  
% end

% % Plotting GDP of each country against time 
% for i=1:28
% %     i
%     X2ti=X(2,:,i)';
%     V1=plot(linspace(1990,2019,20),X2ti);
%     title('GDP');
%     xlabel('year');
%   hold on
%  
% end

% % Plotting DR65 of each country against time 
% for i=1:28
% %     i
%     X3ti=X(3,:,i)';
%     ylim([0,30]);
%     V1=plot(linspace(1990,2019,20),X3ti);
%     title('DR65');
%     xlabel('year');
%   hold on
%  
% end

% % Plotting DR15 of each country against time 
% for i=1:28
%    i
%     X4ti=X(4,:,i)'
%     V1=plot(linspace(1990,2019,20),X4ti);
%     title('DR15');
%     xlabel('year');
%   hold on
%  
% end

% % Plotting HE for each country against time
% for i=1:28
%     y1ti=Y(:,i)';
%     v1=plot(linspace(1990,2019,20),y1ti);
%     title('HE');
%     xlabel('year');
%     hold on 
% end