function EarlyWorld_example_alphabet(Simulation_Name, alpha0, beta0, f_array, T_arr, polymer_size, Ns, Nstart, interval)
%% RNA world 2 is the same of RNA world but is valid for alphabet mixted as 4'
% RNA world is a model to study the replication of a RNA template chain attached to a clay surface.
% clear
%tSart=tic;
%% General inputs
% Simulation_Name='Prueba_betaVariable';
alphabet='2_3';
A='gc'; % Alphabet %A='au'; %A='agucbd';
D=['gc']; % Alphabet combinations matrix %D=['au']; %D=['au';'gc';'bd'];
Nstart=str2num(Nstart);
Ns=str2num(Ns); %1 Number of simulations
if strcmp(alphabet,'2_2') 
    X=4; % 800 of each nt
    label='A2au__';
elseif strcmp(alphabet,'2_3') 
    X=4; % 800 of each nt
    label='A2gc__';
elseif strcmp(alphabet,'4')
    X=2; % 400 of each nt
    label='A4__';
elseif strcmp(alphabet,'4*') 
    X=2; % 400 of each nt
    label='A4ast__';
elseif strcmp(alphabet,'6_2')
    X=1.335; % 267 of each nt
    label='A6_2__';
elseif strcmp(alphabet,'6_3')
    X=1.335; % 267 of each nt
    label='A6_3__';
else
end
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
        beta0=str2num(beta0);
        f_array=str2num(f_array);

        %-- Introduce values from command line

        T_arr=str2num(T_arr);
                                               

%% Random number

% Code generates a new seed each simulation 
    R=0;
% but a selected seed can be introduce (uncomment text below)
    % R=1; 
    % Seed=2;

%% Software Code

%% Create the folder to save the data simulation 

%if not(isfolder('Data'))
path_data='/lustre/home/cab/calejandre/Data/';
mkdir(path_data)
%end
cd(path_data)
%if not(isfolder(Simulation_Name))
mkdir(Simulation_Name)
%end
% cd(Simulation_Name);
% File_Name = strcat('_General_Inputs_b0_',num2str(beta0),'_ATarr_',T_arr,'.mat');
% save(File_Name,'A','D','alpha0','beta_array','A_a','A_b','Ta','Tb','tmax','L','Ns','polymer_size');
path_run='/lustre/home/cab/calejandre/';
cd(path_run);
%% Computes the simulation according to the parameters 

N=X*L*ones(1,size(A,2));

for n=Nstart:Ns 
    %disp(strcat('Simulation_',num2str(n)));
    for T=T_arr
        Ta=T;
        Tb=T;
            for f=f_array
                A_a=f*alpha0;
                A_b=f*beta0;                

                disp(['ALPHABET ',alphabet,' Simulation ',num2str(n),' T ',num2str(Ta),' f ',num2str(f)]);

                tic;

                %-- Set random number generator

                if R==0 % Condition of use a different seed for random numbers
                    seed=sum(1000*clock)+str2num(interval)*100;
                end
                %reset(RandStream.getGlobalStream,seed) % Control of random number in each simulation by seed generated with clock
                rand('seed',seed);

                %-- Compute alpha and beta arrays 
                omega_a=(2*pi)/Ta;
                alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax));
                omega_b=(2*pi)/Tb;
                beta_arr=beta0+A_b*sin(omega_b*(1:tmax));

                [tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=main_alphabets(A,D,N,L,tmax,beta_arr,alpha_arr,str2num(polymer_size),alphabet);

                % Save the simulation
                simu_time=toc;                    
                File_Name = strcat(label,'Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_f_',num2str(f),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat'); 
                path=strcat(path_data,Simulation_Name,'/');
                save(strcat(path,File_Name),'tf','Ut','Ut_OCR','Ut_E','Ct','Ct_OCR','Ct_bl','Ct_E','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')

            end
    end
end