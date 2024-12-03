function EarlyWorld_example_alphabet(Simulation_Name, alpha0, beta0, A_arr, T_arr, polymer_size, Ns_start, Ns_final, interval, tmax)

%% General inputs
alphabet='2_3'; % String code that represents a particular alphabet (in this case a 2-letter alphabet of G-C (3 hydrogen bonds)).
A='gc'; % Alphabet nucleotides
D=['gc']; % Allowed nucleobase complementarities
L=200; % Number of simulated clay positions
tmax=str2num(tmax); % Simulation time steps
Ns_start=str2num(Ns_start); % First realization number
Ns_final=str2num(Ns_final); % Last realization number
alpha0=str2num(alpha0); % Value used to compute the array containing the value of alpha for each time step
beta0=str2num(beta0); % Value used to compute the array containing the value of beta for each time step
A_arr=str2num(A_arr); % Array that stores the amplitudes of the alpha-oscillaitons and the beta-oscillations
T_arr=str2num(T_arr); % Array that stores the periods of the alpha-oscillaitons and the beta-oscillations
A_a=A_arr(1);  % Alpha's amplitude
A_b=A_arr(2);  % Beta's amplitude
Ta=T_arr(1);   % Alpha's period
Tb=T_arr(2);   % Beta's period  

% Initial variables that depend on the alphabet
if strcmp(alphabet,'2_2') 
    X=4; % Variable used in main_alphabets.m to calculate the initial number of nts of each type, so that the total number of initial nts is 1600.
    label='A2au__'; % Label created for file saving purposes
elseif strcmp(alphabet,'2_3') 
    X=4;
    label='A2gc__';
elseif strcmp(alphabet,'4')
    X=2; 
    label='A4__';
elseif strcmp(alphabet,'4*') 
    X=2; 
    label='A4ast__';
elseif strcmp(alphabet,'6_2')
    X=1.335; 
    label='A6_2__';
elseif strcmp(alphabet,'6_3')
    X=1.335; 
    label='A6_3__';
else
end

                                               

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

for n=Ns_start:Ns_final 
                 
    disp(['ALPHABET ',alphabet,', Simulation ',num2str(n), ', tmax ',num2str(tmax),', a0 ',num2str(alpha0),', b0 ',num2str(beta0),', Aa ',num2str(A_a),', Ab ',num2str(A_b),', Ta ',num2str(Ta),', Tb ',num2str(Tb)]);
    
    tic;

    %-- Set random number generator

    if R==0 % Condition of use a different seed for random numbers
        seed=sum(1000*clock)+str2num(interval)*100;
    end
    
    rand('seed',seed);

     % -- Compute alpha and beta arrays that store their values at each time step (in this case they are sinusoidal oscillations)
    omega_a=(2*pi)/Ta;
    alpha_arr=alpha0+A_a*sin(omega_a*(1:tmax));
    omega_b=(2*pi)/Tb;
    beta_arr=beta0+A_b*sin(omega_b*(1:tmax));

    % -- Run the simulation
    [Lt,Lt_OCR,Lt_bl,Lt_E,tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=main_alphabets(A,D,N,L,tmax,beta_arr,alpha_arr,str2num(polymer_size),alphabet);

    % -- Save the simulation
    simu_time=toc;                    
    File_Name = strcat(label,'Sim_','beta0_',num2str(beta0),'_alpha0_',num2str(alpha0),'_AaAbTaTb_',num2str(A_a),'_',num2str(A_b),'_',num2str(Ta),'_',num2str(Tb),'_tmax_',num2str(tmax),'_L_',num2str(L),'_Ns_',num2str(n),'.mat'); path=strcat(path_data,Simulation_Name,'/');
    save(strcat(path,File_Name),'Lt','Lt_OCR','Lt_bl','Lt_E','tf','Ut','Ut_OCR','Ut_E','Ct','Ct_OCR','Ct_bl','Ct_E','Ot_C','Ot_U','Activity','R_CELL','C_CELL','seed','simu_time')

end
