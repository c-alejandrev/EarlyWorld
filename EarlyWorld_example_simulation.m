function EarlyWorld_example_simulation(Simulation_Name, alpha0, beta0, A_arr, T_arr, polymer_size, Ns_start, Ns_final, interval, tmax)

%tSart=tic;
%% General inputs

A='aguc'; % Alphabet %A='au'; %A='agucbd';
D=['au';'gc']; % Alphabet combinations matrix %D=['au']; %D=['au';'gc';'bd'];
X=2; % Proportion of available nucleotides (X) acoording to the length l of the tempalte chain
% Vis='off'; % % Make figures vissible or not (just for cluster running)
L=200; % Length of the total available space
tmax=str2num(tmax); % time of simulation 

%-- Introduce values from command line
Ns_start=str2num(Ns_start);
Ns_final=str2num(Ns_final);
alpha0=str2num(alpha0);
beta0=str2num(beta0);
T_arr=str2num(T_arr);
A_arr=str2num(A_arr);
A_a=A_arr(1); 
A_b=A_arr(2);
Ta=T_arr(1);
Tb=T_arr(2);
        
                                               

%% Random number

% Code generates a new seed each simulation 
    R=0;
% but a selected seed can be introduce (uncomment text below)
    % R=1; 
    % Seed=2;

%% Software Code

%% Create the folder to save the data simulation 

%if not(isfolder('Data'))
%path_data='/lustre/home/cab/calejandre/Data/';
path_data='C:\Users\carla\Desktop\CAB-Garantia\CODE_Github_07_10_2024\Data\';
mkdir(path_data)
%end
cd(path_data)
%if not(isfolder(Simulation_Name))
mkdir(Simulation_Name)
%end
% cd(Simulation_Name);
% File_Name = strcat('_General_Inputs_a0_',num2str(alpha0),'_b0_',num2str(beta0),'_ATarr_',AT_arr,'.mat');
% save(File_Name,'A','D','beta0','alpha0','A_a','A_b','Ta','Tb','tmax','L','Ns','polymer_size');
%path_run='/lustre/home/cab/calejandre/';
path_run='C:\Users\carla\Desktop\CAB-Garantia\CODE_Github_07_10_2024\';
cd(path_run);
%% Computes the simulation according to the parameters 

N=X*L*ones(1,size(A,2));

for n=Ns_start:Ns_final
    %disp(strcat('Simulation_',num2str(n)));
    
    disp(['Simulation ',num2str(n), ', tmax ',num2str(tmax),', a0 ',num2str(alpha0),', b0 ',num2str(beta0),', Aa ',num2str(A_a),', Ab ',num2str(A_b),', Ta ',num2str(Ta),', Tb ',num2str(Tb)]);
    
    tic;

    %-- Set random number generator

    if R==0 % Condition of use a different seed for random numbers
        seed=sum(1000*clock)+str2num(interval)*10000;
    end
    %reset(RandStream.getGlobalStream,seed) % Control of random number in each simulation by seed generated with clock
    rand('seed',seed);

    %-- Compute alpha and beta arrays, can be sinusoidal or constant values
    % Notice that beta=0 correspond to Comp. I simulations
    omega_a=(2*pi)/Ta; 
    alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax));
    omega_b=(2*pi)/Tb;
    beta_arr=beta0+A_b*sin(omega_b*(1:tmax));
    

    [Lt,Lt_OCR,Lt_bl,Lt_E,tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=main(A,D,N,L,tmax,beta_arr,alpha_arr,str2num(polymer_size));

    % Save the simulation
    simu_time=toc;
    
    File_Name = strcat('Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat'); 
    path=strcat(path_data,Simulation_Name,'/');
    save(strcat(path,File_Name),'Lt','Lt_OCR','Lt_bl','Lt_E','tf','Ut','Ut_OCR','Ut_E','Ct','Ct_OCR','Ct_bl','Ct_E','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')
    
    

end
    
