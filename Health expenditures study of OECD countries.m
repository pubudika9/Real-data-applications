%======Confirmation Empirical Application======%
%-----------------OECD HE Data-----------------%

clc;
clear all;
clf("reset");
format longG 

%----Import data----%

[lhe]=xlsread("HEdata/HE.xlsx",'B2:U35'); % per capita health expenditure (dependent variable)
[lghe]=xlsread("HEdata/GHE.xlsx",'B2:U35'); % per capita gvt health expenditure 
[lgdp]=xlsread("HEdata/GDP.xlsx",'B2:U35'); % per capita GDP
[dr65]=xlsread("HEdata/DR65.xlsx",'B2:U35'); % Dependency rate: proportion of population aged above 65
[dr15]=xlsread("HEdata/DR15.xlsx",'B2:U35'); % Dependency rate: proportion of population aged below 15

%gdp14 = xlsread("HEdata/GDP.xlsx",'U2:U35');
%he14 = xlsread("HE.xlsx",'U2:U35');



Z=zeros(5,20,34); % to store Zi
X=zeros(4,20,34); % to store Xi
Y=zeros(20,34); % to store Yi

%-------------------------------------------------------------------------------------------------------------------%

%----Estimation----%

for i=1:34

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

for j=1:34

    z=z+Z(:,:,j);
    

end

Zw=(1/34)*z;

% Creating the projection matrix Mw

Hw=Zw'; % This matrix stores all the observed data

Mw=eye(20)-Hw*inv(Hw'*Hw)*Hw';

%-------------------------------------------------------------------------------------------------------------------%


% Estimating the CCE estimator betai for each i=1:N

 betahat=zeros(4,1); % To store the four estimated coefficients for each i

for i=1:34
    xi=X(:,:,i)';
    yi=Y(:,i);
i;
    betaihat = inv(xi'*Mw*xi)*(xi'*Mw*yi);
    betahat =  [betahat betaihat];
    betaihat = zeros(4,1);
    
  
end

betahat';

% Estimating the CCEMG estimator for each four coefficients
b1 = betahat(1,:)';
b2 = betahat(2,:)';
b3 = betahat(3,:)';
b4 = betahat(4,:)';

% % Plotting the kernal function for GDP coefficients
% [f,xi,bw] = ksdensity(b2)
% figure
% plot(xi,f);
% title("CCE estimation");
% ylabel('Density');

%-------------------------------------------------------------------------------------------------------------------%

weight = (1/34)*ones(1,35);
betaMG =[weight*b1; weight*b2; weight*b3; weight*b4];

betamean1 = mean(betahat(:,2:13)');
betasd1 = std(betahat(:,2:13)');

betamean2 = mean(betahat(:,14:35)');
betasd2 = std(betahat(:,14:35)');

%-------------------------------------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------------------------------------%


% Estimating the covariance matrix of CCEMG (Omega_{v}) using the formula given in the paper for each iteration

    covMG=zeros(4,4);
  
    for i=1:34
     mat1i = betahat(:,i)-betaMG(:,1);
   
     covarianceMGi = covMG + (mat1i)*(mat1i)';
     covMG = covarianceMGi;
      
    end
    estcovMG = (1/33)*covMG;
    
% Estimating the covariance matrix of the CCEMG estimator under the new adjusted assumptions
   newcovMGniter = zeros(4,4);
   
   covMGnew = zeros(4,4);
   for i=1:34
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

   estcovMGnew =  (20/34)*estcovMG + (1/34)*covMGnew;
   
%-------------------------------------------------------------------------------------------------------------------%

   % Confidence interval

   % Critical value at 32 degrees of freedom for each beta

   tc = 2.036933;

%-------------------------------------------------------------------------------------------------------------------%

%----Estimating factors and the loadings----%

u = Y -lghe'* b1(2:35) -lgdp'* b2(2:35) - dr65'*b3(2:35) - dr15'*b4(2:35);
% results=NbFactors(u,5);

%-------------------------------------------------------------------------------------------------------------------%

%-----------------------------Results-------------------------------------%

display(estcovMG);

%-------------------------------------------------------------------------------------------------------------------%

%----Data Visualization----%


%
%--Beta (slope coefficients)--%

% % Plotting the coefficients for gvt health expenditure for each i
% figure(1)
% ylim([-2,2]);
% B1=plot(1:35,b1,'o','color','red');
% title('GHE');
% xlabel('countries');
% ylabel('β_{GHEi}');
% ylim([-2,2]);
% hold on
% plot(1:35,betaMG(1,1),'_','color','blue');
% hold off

%  Plotting the coefficients for GDP for each i

figure(2)
ylim([-2,2]);
B2=plot(1:35,b2,'o','color','red');
title('Estimated Income Elasticities foreach country');
xlabel('countries');
ylabel('β_{GDPi}');
 ylim([-2,2]);
hold on
plot(1:35,betaMG(2,1),'_','color','blue');
hold off;

% % Plotting the coefficients for DR65
% figure(3)
% ylim([-2,2]);
% B2=plot(1:35,b3,'o','color','red');
% title('DR65');
% xlabel('countries');
% ylabel('β_{DR65i}');
% ylim([-2,2]);
% hold on
% plot(1:35,betaMG(3,1),'_','color','blue');
% hold off
% 
% % Plotting the coefficients for DR15
% figure(4)
% ylim([-2,2]);
% B2=plot(1:35,b4,'o','color','red');
% title('DR15');
% xlabel('countries');
% ylabel('β_{DR15i}');
% ylim([-2,2]);
% hold on
% plot(1:35,betaMG(4,1),'_','color','blue');


%--Independent variables--%

%plot(1:35)
% % Plotting GHE of each country against time 
% figure(1);
% for i=1:34
% %     i
%     X2ti=X(1,:,i)'
%     
%     V1=plot(linspace(1990,2019,20),X2ti);
%     title('GHE');
%     xlabel('year');
%   hold on
%  
% end
% hold off
% 
% Plotting GDP of each country against time 
% figure(1)
% for i=1:34
% %     i
%     X2ti=X(2,:,i)';
%     V1=plot(linspace(1996,2014,20),X2ti);
%     title('GDP');
%     xlabel('year');
%   hold on
%  
% end
% hold off
% 
% % Plotting DR65 of each country against time 
% figure(3)
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
% hold off
% 
% % Plotting DR15 of each country against time
% figure(4)
% for i=1:28
%    i
%     X4ti=X(4,:,i)'
%     V1=plot(linspace(1990,2019,20),X4ti);
%     title('DR15');
%     xlabel('year');
%   hold on
%  
% end
% hold off 
% 
% % Plotting HE for each country against time
% figure(5)
% for i=1:28
%     y1ti=Y(:,i)';
%     v1=plot(linspace(1990,2019,20),y1ti);
%     title('HE');
%     xlabel('year');
%     hold on 
% end

% figure(1)
% meanplot=plot(1:34,sort(meanlgdp),'o','color','red');
% title('Average log-GDP');
% xlabel('countries');
% ylabel('mean_{GDPi}');
% % ylim([-2,2]);
% hold on
% 
% plot(1:35,10.30*ones(1,35),'_','color','blue');


%-----------------------------------------------------------------------------------------%














