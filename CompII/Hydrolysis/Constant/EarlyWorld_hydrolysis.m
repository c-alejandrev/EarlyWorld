function EarlyWorld_hydrolysis(Simulation_Name, alpha0, beta0, A_arr, T_arr, polymer_size, Ns_start, Ns_final, interval, tmax, k_hyd_arr)



%% General inputs
A='aguc'; % Alphabet nucleotides
D=['au';'gc']; % Allowed nucleobase complementarities
X=2; % Variable used to calculate the initial number of nts of each type
L=200; % Number of simulated clay positions
tmax=str2num(tmax); % Simulation time steps
Ns_start=str2num(Ns_start); % First realization number
Ns_final=str2num(Ns_final); % Last realization number
alpha0=str2num(alpha0); % Value used to compute the array containing the value of alpha for each time step
beta0=str2num(beta0); % Value used to compute the array containing the value of beta for each time step
T_arr=str2num(T_arr); % Array that stores the periods of the alpha-oscillaitons and the beta-oscillations
A_arr=str2num(A_arr); % Array that stores the amplitudes of the alpha-oscillaitons and the beta-oscillations
A_a=A_arr(1);  % Alpha's amplitude
A_b=A_arr(2);  % Beta's amplitude
Ta=T_arr(1);   % Alpha's period
Tb=T_arr(2);   % Beta's period  
k_hyd_arr=str2num(k_hyd_arr);

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
for k_hyd=k_hyd_arr
    Simulation_Name=num2str(k_hyd);
    mkdir(Simulation_Name);    
end
cd .. 

%% Computes the simulation according to the parameters 

N=X*L*ones(1,size(A,2)); % Number of initial nts of each type

for n=Ns_start:Ns_final
        
    for k_hyd=k_hyd_arr
        disp(['Simulation ',num2str(n), ', k_hyd ',num2str(k_hyd), ', tmax ',num2str(tmax),', a0 ',num2str(alpha0),', b0 ',num2str(beta0),', Aa ',num2str(A_a),', Ab ',num2str(A_b),', Ta ',num2str(Ta),', Tb ',num2str(Tb)]);
        
        tic;

        Simulation_Name=num2str(k_hyd);
    
        % -- Set random number generator
    
        if R==0 % Condition of use a different seed for random numbers
            seed=sum(1000*clock)+str2num(interval)*10000;
        end
        
        rand('seed',seed);
    
        % -- Compute alpha and beta arrays that store their values at each time step (in this case they are sinusoidal oscillations)
        omega_a=(2*pi)/Ta; 
        alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax));
        omega_b=(2*pi)/Tb;
        beta_arr=beta0+A_b*sin(omega_b*(1:tmax));
        
        % -- Run the simulation
        [Lt,Lt_OCR,Lt_bl,Lt_E,tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL,HYDROLYZED,N_pbonds,H_pbonds]=main_hydrolysis(A,D,N,L,tmax,beta_arr,alpha_arr,str2num(polymer_size),k_hyd);
    
        % -- Save the simulation
        simu_time=toc;    
        File_Name = strcat('Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat'); 
        path=strcat(path_data,Simulation_Name,'/');
        save(strcat(path,File_Name),'N_pbonds','H_pbonds','HYDROLYZED','Lt','Lt_OCR','Lt_bl','Lt_E','tf','Ut','Ut_OCR','Ut_E','Ct','Ct_OCR','Ct_bl','Ct_E','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')
    end
    

end
    
