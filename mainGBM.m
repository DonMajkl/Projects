clear; clc; close all;  % clear workspace, clear commands, close all figures

%% Define the Seed
% if you use random numbers at some stage, define a seed 
% (the same random numbers will be drawn after every 'clear' command)
%Picking "mt19937ar" and seed 123456
par.seed.s = RandStream('mt19937ar', 'Seed', 123456);
RandStream.setGlobalStream(par.seed.s);

tic;        % initialize time counter
%% Load data
Stock_Data = xlsread ( "data_pfInsurance.csv", 'B2:B7100');
Bonds_Data = xlsread ( "data_pfInsurance.csv", 'D2:D7100');

% calculate stock returns
stock_returns = (Stock_Data(2:end,1)) ./ (Stock_Data(1:end-1,1))-1;

%% Parameters of Monte Carlo using geometric brownian motion

S0 = Stock_Data(7099,1); % Initial price of DAX
mu = mean(stock_returns); % expected return on DAX
sig = std(stock_returns); % volatility of DAX
dt = 1; % size of time steps
steps = 252-1; % number of days -1 because the initial price S0 is first day
nsims = 1e5; % number of simulated paths 100 000

S = GBM(S0,sig, mu, dt, steps, nsims);

runtime = toc; % stop time counter, save as variable 'runtime'

%% Strategies
Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
Multiplicator_sig = Multiplicator - 1;

% 1= Buy & Hold
% 2= Constant mix
% 3= Stop & loss
% 4= CPPI
% 5= TIPP
% 6= Synthetic put  

Invested_Amount = 100; % units invested
rf = mean(Bonds_Data)/100; % risk-free rate
sig = std(Multiplicator_sig); % volatility of returns
th = 250; % trading days 251 (Elapsed_Time)
ts = 1; % for 250 days
strategy = 6; % type the number of the strategy you wanna use

%% Parameters for strategies
% Adjust parameters for the strategy you wanna use

if strategy==1 % Buy & Hold
    Prct = 30; % percentage of stocks in portfolio (the rest is bonds)
end

if strategy==2 % Constant mix
    Prct = 30; % percentage of stocks in portfolio (the rest is bonds)
end

if strategy==3 % Stop & loss
    Floor = 85; % 85%, 95%
    Prct = 100; % all units invested to stocks at t=0
end

if strategy==4 % CPPI
    Floor = 90; % 85%, 95% 
    Multiple_CPPI = 3; % 3 or 5
end

if strategy==5 % TIPP
    Floor = 85 ; % 85%, 95% 
    Multiple_CPPI = 5; % 3 or 5
end

if strategy==6 % Synthetic put
   Floor = 85; % 85%, 95% 
end

%% Initial values

Stock_Index = Invested_Amount;  
Start = Stock_Index;
Elapsed_Time = 0;
Portfolio_Value1 = Invested_Amount;

if strategy==1 % Buy & Hold
    Stocks_in_Portfolio_Prct = Prct/100;
end

if strategy==2 % Constant mix
    Stocks_in_Portfolio_Prct = Prct/100;
    Bonds_in_Portfolio_Prct = 1 - Stocks_in_Portfolio_Prct;
end

if strategy==3 % Stop & loss
    Stocks_in_Portfolio_Prct = Prct/100;
    Floor = (zeros(nsims,1) + Floor)';
    Portfolio_Value_new = zeros(nsims,1)';
end

if strategy==4 % CPPI
    Cushion = Portfolio_Value1 - Floor;
    Stocks_in_Portfolio = Multiple_CPPI .* Cushion;
    Bonds_in_Portfolio = Portfolio_Value1 - Stocks_in_Portfolio;
    Stocks_in_Portfolio_Prct = Stocks_in_Portfolio ./ Portfolio_Value1;
end

if strategy==5 % TIPP
    Floor = (zeros(nsims,1) + Floor)';
    Floor_prct = Floor / 100;
    Cushion = Portfolio_Value1 - Floor;
    Stocks_in_Portfolio = Multiple_CPPI .* Cushion;
    Bonds_in_Portfolio = Portfolio_Value1 - Stocks_in_Portfolio;
    Stocks_in_Portfolio_Prct = Stocks_in_Portfolio ./ Portfolio_Value1;
end

if strategy==6 % Synthetic put 
    Floor_prct = Floor / 100;
    Stocks_in_Portfolio_Prct = Floor_prct;
end

Stocks_in_Portfolio = Portfolio_Value1 .* Stocks_in_Portfolio_Prct;
Bonds_in_Portfolio = Portfolio_Value1 - Stocks_in_Portfolio;

%% Evolution of assets and strategy rebalancing

while Elapsed_Time <= th
    % Time elapses
    Elapsed_Time = Elapsed_Time + ts;  
    
    % Rebalancing

    if strategy==1 % Buy & Hold

        % New values of assets

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,:);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,:);

        if Elapsed_Time == 62 
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy

        Stocks_in_Portfolio_Prct = Stocks_in_Portfolio ./ Portfolio_Value1;
    end


    if strategy==2 % Constant mix
        
        % New values of assets

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,1:end);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,1:end);

        if Elapsed_Time == 62 
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy
        
        Bonds_in_Portfolio = Portfolio_Value1 .* Bonds_in_Portfolio_Prct;
        Stocks_in_Portfolio = Portfolio_Value1 .* Stocks_in_Portfolio_Prct;
        Stocks_in_Portfolio_Prct = Stocks_in_Portfolio ./ Portfolio_Value1;
    end

    if strategy==3 % Stop & loss

        % New values of assets

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,:);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,:);

        if Elapsed_Time == 62 
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy
          
          Iid_SL = Portfolio_Value1 (1, 1:end) <= Floor;
          Portfolio_Value_new (1,Iid_SL) = 0 + Portfolio_Value1 (1,Iid_SL);
    
            if Elapsed_Time == 62 
            Portfolio_Value_new = Portfolio_Value_new .* (ts+rf); % rf only every 62  days and 4 times during the simulation
            end
            if Elapsed_Time == 124
            Portfolio_Value_new = Portfolio_Value_new .* (ts+rf); % rf only every 62  days and 4 times during the simulation
            end
            if Elapsed_Time == 186
            Portfolio_Value_new = Portfolio_Value_new .* (ts+rf); % rf only every 62  days and 4 times during the simulation
            end
            if Elapsed_Time == 248
            Portfolio_Value_new = Portfolio_Value_new .* (ts+rf); % rf only every 62  days and 4 times during the simulation
            end
        
         Portfolio_Value1 (1:end,Iid_SL) = Portfolio_Value_new (1:end, Iid_SL);
    
         Stocks_in_Portfolio_Prct = Stocks_in_Portfolio ./ Portfolio_Value1;
        
    end
         
    if strategy==4 % CPPI

        % New values of assets

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,:);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,:);

        if Elapsed_Time == 62 
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy
        
        Cushion = Portfolio_Value1 - Floor;
        Stocks_in_Portfolio = Multiple_CPPI .* Cushion;
        
        Iid_CPPI = Stocks_in_Portfolio > Portfolio_Value1; % no short sale constraint
        Stocks_in_Portfolio(1,Iid_CPPI) = Portfolio_Value1(1,Iid_CPPI); % no short sale constraint

        Bonds_in_Portfolio = Portfolio_Value1 - Stocks_in_Portfolio;
        Stocks_in_Portfolio_Prct = Stocks_in_Portfolio./Portfolio_Value1;
    end

    if strategy==5 % TIPP

        % New values of assets

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,:);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,:);

        if Elapsed_Time == 62 
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Floor = Floor .* (ts+rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy
    
        Iid_TIPP = (Portfolio_Value1 (1, 1:end)) .* Floor_prct > Floor;
        Floor_new = ((Portfolio_Value1 (1, 1:end)) .* Floor_prct);
        Floor(1:end,Iid_TIPP) = Floor_new(1:end,Iid_TIPP);  

        Cushion = Portfolio_Value1 - Floor;
        Stocks_in_Portfolio = Multiple_CPPI .* Cushion;

        Iid_TIPP2 = Stocks_in_Portfolio > Portfolio_Value1; % no short sale constraint
        Stocks_in_Portfolio(1,Iid_TIPP2) = Portfolio_Value1(1,Iid_TIPP2); % no short sale constraint

        Bonds_in_Portfolio = Portfolio_Value1 - Stocks_in_Portfolio;
        Stocks_in_Portfolio_Prct = Stocks_in_Portfolio./Portfolio_Value1;
    end

    if strategy==6 % Synthetic put

        Multiplicator = (S(2:end,1:end)) ./ (S(1:end-1,1:end));
        Stock_Index = Stock_Index .* Multiplicator(0+Elapsed_Time,:);
        Stocks_in_Portfolio = Stocks_in_Portfolio .* Multiplicator(0+Elapsed_Time,:);

        if Elapsed_Time == 62 
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 124
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 186
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end
        if Elapsed_Time == 248
        Bonds_in_Portfolio = Bonds_in_Portfolio .* (ts + rf); % rf only every 62  days and 4 times during the simulation
        end

        Portfolio_Value1 = Stocks_in_Portfolio + Bonds_in_Portfolio;

        % rebalacing of strategy

        Time_to_Maturity = (th - Elapsed_Time + 1) / 365;

        d = (log(Stocks_in_Portfolio ./ Floor) + (rf + 0.5 .* sig.^2) .* Time_to_Maturity) ./ sig .* sqrt(Time_to_Maturity);

        Stocks_in_Portfolio_Prct = (Stocks_in_Portfolio .* normcdf(d)) ./ (Stocks_in_Portfolio .* normcdf(d) + Floor .* exp(-rf * Time_to_Maturity ) .* normcdf(sig .* sqrt(Time_to_Maturity) - d));
        Bonds_in_Portfolio_Prct = 1 - Stocks_in_Portfolio_Prct; 

        Bonds_in_Portfolio = Portfolio_Value1 .* Bonds_in_Portfolio_Prct;
        Stocks_in_Portfolio = Portfolio_Value1 .* Stocks_in_Portfolio_Prct;    
    end

end

%% Evaluation of the results

% Wealth

Portfolio_Value1 = Portfolio_Value1';
End_of_Period_Wealth = mean(Portfolio_Value1);
Average_Portfolio_Return_Prct = mean(Portfolio_Value1) -100;

% Standard deviation

Stdev_End_of_Period_Wealth = std(Portfolio_Value1);

% plot
PR_plot = mean(Multiplicator_sig);
figure;
subplot(1,2,1)
h=histogram(Multiplicator_sig,50);
h.FaceColor = [0 0 0];
    set(gca, 'FontSize', 14); 
    title('Distribution of the returns');
    xlabel('Returns GBM')
    ylabel('Volume')

plot_x = S(252, 1:end); 

figure;
scatter(plot_x, Portfolio_Value1, 10, 'red')
      set(gca, 'FontSize', 14);
      title('Distribution of the returns, CPPI');
      xlabel('DAX price at time T')
      ylabel('Portfolio value')


