function EarlyWorld_example_simulation(Simulation_Name, alpha0, AT_arr, Nstart, Nsfin,interval)
%% RNA world 2 is the same of RNA world but is valid for alphabet mixted as 4'
% RNA world is a model to study the replication of a RNA template chain attached to a clay surface.
% clear
tSart=tic;
%% General inputs
% Simulation_Name='Prueba_betaVariable';
A='aguc'; % Alphabet %A='au'; %A='agucbd';
D=['au';'gc']; % Alphabet combinations matrix %D=['au']; %D=['au';'gc';'bd'];
Ns=100; %1 Number of simulations
X=2; % Proportion of available nucleotides (X) acoording to the length l of the tempalte chain
% Vis='off'; % % Make figures vissible or not (just for cluster running)

%% Template Chain
% L can be a scalar (1 Length) or a vector (study of different lengths)
L=200; % Length of the total available space

%% Time (T_max) can be a scalar (time of simulation) or a vector (study of different times)

tmax=20000; % time of simulation 

%% Beta can be constant or a function of time
% If beta is constant uncomment the lines below and write the scalar or the values of beta that you want to study.
%      beta=0.10:0.15:1.50; cte=1;       
       %beta=1; cte=1;


% If beta is a function of time uncomment the lines below and set the function beta(t)
        cte=0;

        %--  Constant values
        alpha0=str2num(alpha0);
        beta0=0;
        Nstart=str2num(Nstart);
        Nsfin=str2num(Nsfin);

        %-- Introduce values from command line

        A_a_T_a=str2num(AT_arr);
        A_a=A_a_T_a(1); 
        Ta=A_a_T_a(2);
        A_b=0;
        Tb=0;
        fi0_a=0;
        fi0_b=0;
                                               

%% Random number

% Code generates a new seed each simulation 
    R=0;
% but a selected seed can be introduce (uncomment text below)
    % R=1; 
    % Seed=2;

%% Software Code

%% Create the folder to save the data simulation 

%if not(isfolder('Data'))
path_data='./Data/';
mkdir(path_data)
%end
cd(path_data)
%if not(isfolder(Simulation_Name))
mkdir(Simulation_Name)
%end
% cd(Simulation_Name);
% File_Name = strcat('_General_Inputs_a0_',num2str(alpha0),'_b0_',num2str(beta0),'_ATarr_',AT_arr,'.mat');
% save(File_Name,'A','D','beta0','alpha0','A_a','A_b','Ta','Tb','tmax','L','Ns','polymer_size');
%path_run='C:\Users\carla\Desktop\CAB-Garantia\CODE_Github_07_10_2024\CompI\';
%cd(path_run);
cd .. 

%% Computes the simulation according to the parameters 

N=X*L*ones(1,size(A,2));

for n=Nstart:Nsfin 
    %disp(strcat('Simulation_',num2str(n)));


        disp(['Simulation ',num2str(n),'alpha_0',num2str(alpha0)]);

        tic;

        %-- Set random number generator

        if R==0 % Condition of use a different seed for random numbers
            seed=sum(1000*clock)+str2num(interval)*100;
        end
        %reset(RandStream.getGlobalStream,seed) % Control of random number in each simulation by seed generated with clock
        rand('seed',seed);

        %-- Compute alpha and beta arrays 
%                     alpha_arr=[];
%                     alpha_arr(1:tmax)=alpha0;
        
        beta_arr=[];
        beta_arr(1:tmax)=0;

        omega_a=(2*pi)/Ta; 
        alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax));

        [Lt,Lt_OCR,Lt_bl,tf,Ut,Ut_OCR,Ct,Ct_OCR,Ct_bl,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=main(A,D,N,L,tmax,beta_arr,alpha_arr);

        % Save the simulation
        simu_time=toc;
        %File_Name = strcat('Sim_','beta_',num2str(beta0),'_alpha_',num2str(alpha0),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat'); 
        File_Name = strcat('Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_fi0a_',num2str(fi0_a),'_fi0b_',num2str(fi0_b),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat');
        path=strcat(path_data,Simulation_Name,'/');
        save(strcat(path,File_Name),'Lt','Lt_OCR','Lt_bl','tf','Ut','Ut_OCR','Ct','Ct_OCR','Ct_bl','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')
                  
       
end
