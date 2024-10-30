function EarlyWorld_example_simulation(Simulation_Name, alpha0, AT_arr, Nstart, Nsfin,interval)

%% General inputs
A='aguc'; % Alphabet nucleotides
D=['au';'gc']; % Allowed nucleobase complementarities
X=2; % Variable used to calculate the initial number of nts of each type
L=200; % Number of simulated clay positions
tmax=20000; % Simulation time steps
alpha0=str2num(alpha0); % Value used to compute the array containing the value of alpha for each time step
beta0=0; % Beta is 0 in CompI simulations, because template-dependent replication is not allowed
Nstart=str2num(Nstart);
Nsfin=str2num(Nsfin);
A_a_T_a=str2num(AT_arr); % Array that stores the amplitude and period of the alpha-oscillations (supposing alpha oscillates with time)
A_a=A_a_T_a(1); % Extracting the amplitude
Ta=A_a_T_a(2); % Extracting the period
A_b=0; % Beta related oscillation: irrelevant for CompI simulations.
Tb=0; % Beta related oscillation: irrelevant for CompI simulations.
%fi0_a=0;
%fi0_b=0;
                                               

%% Random number generator

% Code generates a new seed each simulation 
    R=0;
% but a selected seed can be introduce (uncomment text below)
    % R=1; 
    % Seed=2;

%% Create the folder to save the simulated data 

path_data='./Data/';
mkdir(path_data)
cd(path_data)
mkdir(Simulation_Name)
cd .. 

%% Computes the simulation according to the parameters 

N=X*L*ones(1,size(A,2)); % Number of initial nts of each type

for n=Nstart:Nsfin 

        disp(['Simulation ',num2str(n),'alpha_0',num2str(alpha0)]);

        tic;

        % -- Set random number generator

        if R==0 % Condition to use a different seed for every realization
            seed=sum(1000*clock)+str2num(interval)*100;
        end
        
        rand('seed',seed);

        % -- Compute alpha and beta arrays that store their values at each time step
        
        beta_arr=[];
        beta_arr(1:tmax)=0; % Beta is 0 in CompI simulations

        omega_a=(2*pi)/Ta; 
        alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax)); % Sinusoidal alpha oscillation

        % -- Run the simulation
        [Lt,Lt_OCR,Lt_bl,tf,Ut,Ut_OCR,Ct,Ct_OCR,Ct_bl,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=main(A,D,N,L,tmax,beta_arr,alpha_arr);

        % -- Save the simulation
        simu_time=toc;
        File_Name = strcat('Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat');
        path=strcat(path_data,Simulation_Name,'/');
        save(strcat(path,File_Name),'Lt','Lt_OCR','Lt_bl','tf','Ut','Ut_OCR','Ct','Ct_OCR','Ct_bl','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')
                  
       
end
