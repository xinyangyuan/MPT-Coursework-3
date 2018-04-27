 function [] = ageing_module(STAMPING_VELOCITY,POISSON_RATIO,temperature,time,LEAN_MODE)

%{
Ageing Model
File Format:
You should put all of your (1) major strain file and (2) minor strain file  
at all stages in one folder, you can choose any folder name you want. 

But the major strain file has to be names as 'major_strain_' plus stage 
number, e.g., for stage one it will be 'major_strain_1'. Similarly for
minor strain, it will be 'minor_strain_1'.

Arguments:
STAMPING_VELOCITY
            -- velocity at stamping stage, a scalar
POISSON_RATIO
            -- poisson ratio of material, should in range [-0.5,1]
               (default for alumninium is 0.33), a scalar
temperature
            -- oven temperatures in Kelvin of ageing stages, a string

               For example, if the formed part is first aged in 500K oven
               for 600 seconds, then transfered to an oven with 450K and
               aged for 36000 seconds. The input string should be:
               '500,450' or '500 450'

               At least one temperature has to be given, and the
               temperature should be in range of [440,520] for accurate
               results
time
            -- ageing in seconds of ageing stages, a string

               For example, if the formed part is first aged in 500K oven
               for 600 seconds, then transfered to an oven with 450K and
               aged for 36000 seconds. The input string should be:
               '600,36000' or '600 36000'
LEAN_MODE
            -- lean mode setting, a boolean input (i.e. 0 or 1)
               recommended setting is turned on lean monde (i.e.
               LEAN_MODE=1). Module runs faster in lean mode, because
               several viscoplastic material values are not stored under
               lean mode (~2 minute for 10,000 elements)
Returns:
Yield_Strength.txt
            -- text file contains element number, yield strength [HV] of
               each element after ageing
Final_Yield_Strength_Distribution.fig
            -- matlab plot for material post ageing yield strength
               distribution
Yield_Strength_Distribution_Evolution.fig
            -- material yield strength distribution at eight different
               times during ageing process
Yield_Strength_Evolution.fig
            -- yield strength evolution of all elements during aging

%**********************************************************************
Author: https://github.com/xinyangyuan
%**********************************************************************
%}


%--------------------------------------------------------------------------
% Input Arguments Handle
%--------------------------------------------------------------------------
try
    rmdir ('Viscoplastic','s');
end

ERROR_FLAG = 0;
temperature = str2num(temperature);
time = str2num(time);

if numel(temperature) == 0
    fprintf(['Argument Error | Error argument name(s): temperature',...
            '\n', ...
            'Error info: you have input at least one temperature value',...
            '\n']);
    ERROR_FLAG = 1;
elseif numel(temperature) == 1
    MULTI_STAGE = 0;
elseif numel(temperature) > 1
    MULTI_STAGE = 1;
end

if numel(time) == 0
    fprintf(['Argument Error | Error argument name(s): time',...
            '\n', ...
            'Error info: you have input at least one time value',...
            '\n']);
    ERROR_FLAG = 1;
end

if numel(time) ~= numel(temperature)
    fprintf(['Argument Error | Error argument name(s): time,temperature',...
            '\n', ...
            'Error info: the number of temperature has to be equal to number of time',...
            '\n']);
    ERROR_FLAG = 1;
end

if LEAN_MODE ~= 0 && LEAN_MODE ~= 1
    fprintf(['Argument Error | Error argument name(s): LEAN_MODE',...
            '\n', ...
            'Error info: the input argument LEAN_MODE has to be 1 or 0',...
            '\n']);
    ERROR_FLAG = 1;
end


if (POISSON_RATIO > 0.5) || (POISSON_RATIO < -1)
    fprintf(['Argument Error | Error argument name(s): POISSON_RATIO',...
            '\n', ...
            'Error info: POISSON_RATIO has to be in the range of [-1,0.5]',...
            '\n']);
    ERROR_FLAG = 1;
end

if (max(temperature) > 520) || (max(temperature) < 440)
    fprintf(['Warning | Warning argument name(s): temperature',...
            '\n', ...
            'Warning info:  temperature should be in the range of [440,520] for accurate results',...
            '\n']);
end

% Stop function due to error
if ERROR_FLAG == 1
    error('Function Error Occured')
end

%--------------------------------------------------------------------------
% Data Load
%--------------------------------------------------------------------------
[NUMBER_STAGES,element_no,major_strain,minor_strain,closing_distance] = ...
                                                                loadData();
[strain] = ...
           cal_equivilent_strain(major_strain,minor_strain,POISSON_RATIO);
strain = strain';

% problem dimension
LEAN_MODE        = boolean(LEAN_MODE);
MULTI_STAGE      = boolean(MULTI_STAGE);
STAMPING_SPEED   = STAMPING_VELOCITY;
NUMBER_STAGES    = size(strain,1);  % n_t
NUMBER_ELEMENTS  = size(strain,2);  % m

% forming information
time_between_stages = diff(closing_distance) ./ STAMPING_SPEED;
time_between_stages = ones((NUMBER_STAGES-1),NUMBER_ELEMENTS) .* ...
                      time_between_stages; % eliminate future broadcasting
strain_rate = diff(strain,1,1) ./ time_between_stages; % calculate strain rate

%--------------------------------------------------------------------------
% Viscoplastic
%--------------------------------------------------------------------------
% viscoplastic model material constants

% REMOVED DUE TO IP
% w_viscoplastic = [******];
%

% main viscoplastic
if LEAN_MODE == false
    % calculate step size
    step_size = ones(size(strain_rate)) .* 0.00001;
    step_size(abs(strain_rate) <= 5.0)  = 0.00005;
    step_size(abs(strain_rate) <= 2.0)  = 0.0001;
    step_size(abs(strain_rate) <= 1.0)  = 0.0002;
    step_size(abs(strain_rate) <= 0.5)  = 0.0004;
    step_size(abs(strain_rate) <= 0.5)  = 0.0005;
    step_size(abs(strain_rate) <= 0.1)  = 0.001;
    step_size(abs(strain_rate) <= 0.05) = 0.003;
    step_size(abs(strain_rate) <= 0.01) = 0.008;
    step_size(abs(strain_rate) <= 0.001)= 0.05;

    % preallocation
    rho_final     = zeros(1,NUMBER_ELEMENTS);
    stress_history = cell(1,NUMBER_ELEMENTS);
    strain_history = cell(1,NUMBER_ELEMENTS);
    ep_history     = cell(1,NUMBER_ELEMENTS);
    rho_history    = cell(1,NUMBER_ELEMENTS);
    R_history      = cell(1,NUMBER_ELEMENTS);

    [stress_history{:}] = deal([]);
    [strain_history{:}] = deal([]);
    [ep_history{:}] = deal([]);
    [rho_history{:}] = deal([]);
    [R_history{:}] = deal([]);

    % main loop
    parfor element_idx = 1:NUMBER_ELEMENTS
        initial_condition = [0,0,0,0];
        for time_idx = 1:(NUMBER_STAGES-1)
            [stress_stage, strain_stage, ep_stage, rho_stage, R_stage] = ...
            viscoplastic_model_predict(w_viscoplastic,...
                                       initial_condition,...
                                       strain_rate(time_idx,element_idx),...
                                       time_between_stages(time_idx,element_idx),...
                                       step_size(time_idx,element_idx)...
                                       );
            % update next stage initial condition to the end of previous
            % stage
            initial_condition = ...
                [rho_stage(end),R_stage(end),ep_stage(end),strain_stage(end)];

            % record final stamping normalized dislocation
            if time_idx == (NUMBER_STAGES-1)
                rho_final(1,element_idx) = rho_stage(end);
            end

            % recod history
            stress_history{element_idx} = [stress_history{element_idx};stress_stage];
            strain_history{element_idx} = [strain_history{element_idx};strain_stage];
            ep_history{element_idx}     = [ep_history{element_idx};ep_stage];
            rho_history{element_idx}    = [rho_history{element_idx};rho_stage];
            R_history{element_idx}      = [R_history{element_idx};R_stage];
        end
    end

else
    % calculate step size
    step_size = ones(size(strain_rate)) .* 0.00001;
    step_size(abs(strain_rate) <= 5.0)  = 0.00005;
    step_size(abs(strain_rate) <= 2.0)  = 0.0001;
    step_size(abs(strain_rate) <= 1.0)  = 0.0002;
    step_size(abs(strain_rate) <= 0.5)  = 0.0004;
    step_size(abs(strain_rate) <= 0.5)  = 0.0005;
    step_size(abs(strain_rate) <= 0.1)  = 0.0008;
    step_size(abs(strain_rate) <= 0.05) = 0.002;
    step_size(abs(strain_rate) <= 0.01) = 0.005;
    step_size(abs(strain_rate) <= 0.001)= 0.02;

    % preallocation
    rho_final  = zeros(1,NUMBER_ELEMENTS);

    % main loop
    parfor element_idx = 1:NUMBER_ELEMENTS
        initial_condition = [0,0,0,0];
        for time_idx = 1:(NUMBER_STAGES-1)

            [stress_stage_end,...
             strain_stage_end, ...
             ep_stage_end, ...
             rho_stage_end, ...
             R_stage_end] = viscoplastic_model_predict_lean ...
                           (w_viscoplastic,...
                            initial_condition,...
                            strain_rate(time_idx,element_idx),...
                            time_between_stages(time_idx,element_idx),...
                            step_size(time_idx,element_idx)...
                            );

            % set next stage initial condition to the end of previous
            % stage
            initial_condition = ...
                [rho_stage_end,R_stage_end,ep_stage_end,strain_stage_end];

            % record final stamping normalized dislocation
            if time_idx == (NUMBER_STAGES-1)
                rho_final(1,element_idx) = rho_stage_end;
            end
        end
    end
end
%save('rho_viscoplastic','rho_final')

%--------------------------------------------------------------------------
% Ageing
%--------------------------------------------------------------------------
% load data
rho_0 = rho_final;

% ageing model material constants

% REMOVED DUE TO IP
% w_ageing = [*********];
% 

% numerical time step
TIMESTEP = 5; % default to 5 seconds

if MULTI_STAGE == false
    [time_steps,sig_y] = ...
            ageing_model_predict(rho_0,...
                                 temperature,...
                                 time,...
                                 TIMESTEP,...
                                 w_ageing);

%--------------------------------------------------------------------------
% Multi-stage Ageing
%--------------------------------------------------------------------------
elseif MULTI_STAGE == true
    temperature_multi = temperature(1:(end-1));
    temperature_last  = temperature(end);
    time_multi = time(1:(end-1));
    time_last  = time(end);

    [time_steps,sig_y] = ...
            multistage_ageing_model_predict(rho_0,...
                                            temperature_multi,...
                                            time_multi,...
                                            temperature_last,...
                                            time_last,...
                                            TIMESTEP,...
                                            w_ageing);
end
%--------------------------------------------------------------------------
% Data Output
%--------------------------------------------------------------------------
% save final yield strength data
Yield_Strength = sig_y(end,:)';
Element_No = element_no;
output_data = table(Element_No,Yield_Strength);
% text format
writetable(output_data,'Yield_Strength.txt','Delimiter',' ')
% csv format
writetable(output_data,'Yield_Strength.csv','Delimiter',',','QuoteStrings',true)

% yield strength evolution plot
figure_1 = figure; % ('visible','off')
hold on
for i = 1:400
    plot3(time_steps,i+ones(numel(time_steps),1),sig_y(:,i));
end
xlabel('Time [sec]')
ylabel('Element No.')
zlabel('Yield Strength [HV]')
view(3)
saveas(figure_1,'Yield_Strength_Evolution','fig');

% yield strength evolution histgram
time_length = numel(time_steps);
t_1 = ceil(time_length/16);
t_2 = ceil(time_length/12);
t_3 = ceil(time_length/8);
t_4 = ceil(time_length/4);
t_5 = ceil(time_length*3/8);
t_6 = ceil(time_length/2);
t_7 = ceil(time_length*3/4);
t_8 = time_steps(end);

figure_2 = figure;
hold on
h1 = histogram(sig_y(t_1,:),30);
h2 = histogram(sig_y(t_2,:),30);
h3 = histogram(sig_y(t_3,:),30);
h4 = histogram(sig_y(t_4,:),30);
h5 = histogram(sig_y(t_5,:),30);
h6 = histogram(sig_y(t_6,:),30);
h7 = histogram(sig_y(t_7,:),30);
h8 = histogram(Yield_Strength,30);

h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h4.Normalization = 'probability';
h5.Normalization = 'probability';
h6.Normalization = 'probability';
h7.Normalization = 'probability';
h8.Normalization = 'probability';

h1.BinWidth = 0.25;
h2.BinWidth = 0.25;
h3.BinWidth = 0.25;
h4.BinWidth = 0.25;
h5.BinWidth = 0.25;
h6.BinWidth = 0.25;
h7.BinWidth = 0.25;
h8.BinWidth = 0.25;

xlabel('Yield Strength [HV]')
ylabel('Propability')
legend(strcat(num2str(t_1),' sec'),...
       strcat(num2str(t_2),' sec'),...
       strcat(num2str(t_3),' sec'),...
       strcat(num2str(t_4),' sec'),...
       strcat(num2str(t_5),' sec'),...
       strcat(num2str(t_6),' sec'),...
       strcat(num2str(t_7),' sec'),...
       strcat(num2str(t_8),' sec')...
       );

saveas(figure_2,'Yield_Strength_Distribution_Evolution','fig');

% final yield strength distribution histgram
figure_3 = figure;
histogram(Yield_Strength,30);
xlabel('Yield Strength [HV]')
ylabel('Number of Elements')
saveas(figure_3,'Final_Yield_Strength_Distribution','fig');

% save viscoplatic information
if LEAN_MODE ==false
    mkdir('Viscoplastic')
    save Viscoplastic/Stress.mat stress_history
    save Viscoplastic/Strain.mat strain_history
    save Viscoplastic/Plastic_Strain.mat ep_history
    save Viscoplastic/Normalized_Dislocation.mat rho_history
    save Viscoplastic/Dislocation_Strengthening.mat R_history
end

end




%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------
% load data
function [NUMBER_STAGES,element_no,major_strain,minor_strain,closing_distance] = loadData()

    % loacate pamstamp data folder directory
    Files = dir;
    dirFlags = [Files.isdir];
    Folders = Files(dirFlags);
    Folders(1:2)=[];% remove '.', '..' folder directories

    % find files in the pamstamp data folder
    files = dir(Folders.name);%files in the folder
    files(1:2)=[];%remove meaningless data
    NUMBER_STAGES = size(files,1)/2;%number of simulation steps. Assume only major and minor strain data of each step

    % extract element numbers (from the first ASCII file)
    dataname=strcat(Folders.name,'/major_strain_1.asc');
    data=importdata(dataname,' ',9);  %ignore 9 lines of data, load data
    value_of_data=data.data(:,1);     %First column of data, element number
    value_of_data(end,:)=[];          %remove last 0

    NUMBER_ELEMENTS=numel(value_of_data);%number of elements = number of array elements
    element_no=value_of_data;         %Element numbers
    clear value_of_data


    % Preallocation
    major_strain = zeros(NUMBER_ELEMENTS,NUMBER_STAGES);
    minor_strain = zeros(NUMBER_ELEMENTS,NUMBER_STAGES);
    closing_distance = zeros(NUMBER_STAGES,1);

    %extract the element number

    for n=1:NUMBER_STAGES
       % obtain major strain
       dataname=strcat(Folders.name,'/major_strain_',num2str(n),'.asc');%read major strain of step n
       data=importdata(dataname,' ',9);%ignore 9 lines of data, load data
       value_of_data_1=data.data(:,2);
       value_of_data_1(end,:)=[];     %remove last NaN
       major_strain(:,n) = value_of_data_1;

       % obtain minor strain
       dataname=strcat(Folders.name,'/minor_strain_',num2str(n),'.asc');
       data=importdata(dataname,' ',9);
       value_of_data_2=data.data(:,2);
       value_of_data_2(end,:)=[];
       minor_strain(:,n) = value_of_data_2;

       % obtain preogressive stamping distance
       dataname=strcat(Folders.name,'/major_strain_',num2str(n),'.asc');%read major strain of step n
       data=importdata(dataname,' ',7);%ignore 7 lines of data, load data
       value_of_data_3=data.data;
       closing_distance(n,1) = value_of_data_3;

       clear value_of_data_1 value_of_data_2 value_of_data_3
    end
     clear value_of_data_1 value_of_data_2 value_of_data_3 data dataname
end

% equivilent strain
function [equivilent_strain] = cal_equivilent_strain(e1,e2,POISSON_RATIO)
    e3 = 0 - e1 - e2;
    equivilent_strain = ...
       1/(1+POISSON_RATIO).*(0.5.*((e1-e2).^2+(e2-e3).^2+(e3-e1).^2)).^0.5;
end

%----------------------------Material Models-------------------------------
% viscoplastic model
function [stress, strain, ep, rho, R] = viscoplastic_model_predict...
                                           (w,...
                                            initial_condition,...
                                            STRAIN_RATE,...
                                            TOTAL_TIME,...
                                            TIMESTEP)
    %{
    Vectorized Implementation of the Viscoplastic Model

    Arguments:
    w      -- model parameters, an array (n_w,1)
    initial_condition
           -- initial values for dislocation, dis strengthening, plastic
              strain, total strain, [rho_0,R_0,ep_0,strain_0], an array
              (1, 4)
    STRAIN_RATE
           -- the strain rate between the two stages, a scalar
    TOTAL_TIME
           -- total time between the two stages, a scalar
    TIMESTEP
           -- timestep for numerical integration, a scalar

    Returns:
    stress -- stress (TOTAL_TIME/TIMESTEP,1)
    strain -- total strain (TOTAL_TIME/TIMESTEP,1)
    ep     -- plastic strain (TOTAL_TIME/TIMESTEP,1)
    rho    -- normalized dislocation (TOTAL_TIME/TIMESTEP,1)
    R      -- dislocation/plastic hardening (TOTAL_TIME/TIMESTEP,1)
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    % Retrive iniitial conditions
    initial_condition = num2cell(initial_condition);
    [rho_0,R_0,ep_0,strain_0] = deal(initial_condition{:});


    % Retrive model parameters
    n1 = w(1,:);
    n2 = w(2,:);
    A = w(3,:);
    B = w(4,:);
    C = w(5,:);
    E = w(6,:);
    k = w(7,:);
    K = w(8,:);


    % Pobelem dimension
    time = (0:TIMESTEP:ceil(TOTAL_TIME/TIMESTEP)*TIMESTEP)';
    strain = (time .* STRAIN_RATE) + strain_0;
    TEMPORAL_SIZE = length(time);  % n_t

    % Initilize Variables
    rho = zeros(TEMPORAL_SIZE, 1);
    R   = zeros(TEMPORAL_SIZE, 1); % R = B*rho^(0.5)
    ep  = zeros(TEMPORAL_SIZE, 1);

    if rho_0 ==0
        rho(1,:) = realmin;% realmin to approx 0, since 0 -> math error
    else
        rho(1,:) = rho_0;
    end
    R(1,:) = R_0;
    ep(1,:) = ep_0;

    % Main loop
    for i = 1:length(strain)-1
        % rate update
        dep = ((E.*(strain(i)-ep(i,:))-R(i,:)-k)./K).^(n1); % clipping dep
        dep = (STRAIN_RATE >= 0) * min(dep, STRAIN_RATE) + ...
              (STRAIN_RATE <  0) * max(dep, STRAIN_RATE);

        drho = A.*(1-rho(i,:)).*dep - (C.*rho(i,:).^(n2));

        % value update
        ep(i+1,:)  = eulerForward(ep(i,:),dep,TIMESTEP);
        rho(i+1,:) = max(...
                         eulerForward(rho(i,:),drho,TIMESTEP),...
                         realmin ...  % normalized rho cannot be negative
                         );
        rho(i+1,:) = min(rho(i+1,:) , 1);
        R(i+1,:)   = B .* rho(i+1,:).^(0.5);

    end

    % Output
    ratio = (TOTAL_TIME - time(end-1))/(time(end) - time(end-1));

    ep(end,:)  = ep(end-1,:) + ratio * (ep(end,:) - ep(end-1,:));
    rho(end,:) = rho(end-1,:)+ ratio * (rho(end,:) - rho(end-1,:));
    R(end,:)   = R(end-1,:)  + ratio * (R(end,:) - R(end-1,:));
    strain(end,:) = strain(end-1,:)+ ratio*(strain(end,:)-strain(end-1,:));
    stress = E .* (strain - ep);

end
% viscoplastic lean
function [stress_, strain_, ep_, rho_, R_] = viscoplastic_model_predict_lean...
                                            (w,...
                                             initial_condition,...
                                             STRAIN_RATE,...
                                             TOTAL_TIME,...
                                             TIMESTEP)
    %{
    Vectorized Implementation of the Viscoplastic Model, Lean Verzion

    Arguments:
    w      -- model parameters, an array (n_w,1)
    initial_condition
           -- initial values for dislocation, dis strengthening, plastic
              strain, total strain, [rho_0,R_0,ep_0,strain_0], an array
              (1, 4)
    STRAIN_RATE
           -- the strain rate between the two stages, a scalar
    TOTAL_TIME
           -- total time between the two stages, a scalar
    TIMESTEP
           -- timestep for numerical integration, a scalar

    Returns:
    stress_ -- stress, a scalar
    strain_ -- total strain, a scalar
    ep_     -- plastic strain, a scalar
    rho_    -- normalized dislocation, a scalar
    R_      -- dislocation/plastic hardening, a scalar
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    % Retrive iniitial conditions
    initial_condition = num2cell(initial_condition);
    [rho_0,R_0,ep_0,strain_0] = deal(initial_condition{:});


    % Retrive model parameters
    n1 = w(1,:);
    n2 = w(2,:);
    A = w(3,:);
    B = w(4,:);
    C = w(5,:);
    E = w(6,:);
    k = w(7,:);
    K = w(8,:);


    % Pobelem dimension
    time = (0:TIMESTEP:ceil(TOTAL_TIME/TIMESTEP)*TIMESTEP)';
    strain = (time .* STRAIN_RATE) + strain_0;
    TEMPORAL_SIZE = length(time);  % n_t


    % Initilize Variables
    R_ = R_0;
    ep_ = ep_0;
    rho_ = max(realmin,rho_0);
    ep_pre = zeros(1,1); % pre-defined variable, otherwiser parfor error
    rho_pre= zeros(1,1);
    R_pre  = zeros(1,1);



    % Main loop
    for i = 1:TEMPORAL_SIZE-1
        % rate update
        dep = ((E.*(strain(i)-ep_)-R_-k)./K).^(n1); % clipping dep value
        dep = (STRAIN_RATE >= 0) * min(dep, STRAIN_RATE) + ...
              (STRAIN_RATE <  0) * max(dep, STRAIN_RATE);

        drho = A.*(1-rho_).*dep - (C.*rho_.^(n2));

        % value update
        ep_  = eulerForward(ep_,dep,TIMESTEP);
        rho_ = max(...
                  eulerForward(rho_,drho,TIMESTEP),...
                  realmin...     % normalized rho cannot be negative
                  );
        rho_ = min(rho_ , 1);
        R_   = B .* rho_.^(0.5);

        % recod the second last step values
        if i == (TEMPORAL_SIZE-2)
            ep_pre = ep_;
            rho_pre = rho_;
            R_pre = R_;
        end
    end

    % Output (intropolation of the value at t = TOTAL_TIME)
    ratio = (TOTAL_TIME - time(end-1))/(time(end) - time(end-1));

    ep_ = ep_pre + (ep_-ep_pre) * ratio;
    rho_ = rho_pre + (rho_- rho_pre) * ratio;
    R_ = R_pre + (R_-R_pre) * ratio;
    strain_ = strain(end-1) + (strain(end)-strain(end-1)) * ratio;
    stress_ = E * (strain_ - ep_);

end
% single stage ageing
function [time,sig_y] = ageing_model_predict...
                      (rho_0,TEMPERATURE,TOTAL_TIME,TIMESTEP,w)
    %{
    Implementation of the Ageing Model
    Arguments:
    rho_0 -- initial normalized dislocation, an array (1,ELEMENT_SIZE)
    C_0   -- initial solute fraction, an array (1,ELEMENT_SIZE)
    f_0   -- initial volume fraction, an array (1,ELEMENT_SIZE)
    r_0   -- initial mean precipitate size, an array (1,ELEMENT_SIZE)
    T     -- temperature of the ageing, a scalar
    TOTAL_TIME -- total length of the ageing process, a scalar
    tp    -- peak time, a scalar
    TIMESTEP -- numerical integration time step, a scalar
    w      -- model parameters, an array (n_w,1)

    Returns:
    sig_y
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    % Problem Dimension
    ELEMENT_SIZE = size(rho_0,2); % m
    TEMPORAL_SIZE = ceil(TOTAL_TIME/TIMESTEP) + 1;  % n_t

    % Physical Constant
    R = 8.314;

    % Retrive Training Parameters
    A1 = w(1,:);
    A2 = w(2,:);
    A_d = w(3,:);
    B1 = w(4,:);
    B2 = w(5,:);
    C1 = w(6,:);
    C2 = w(7,:);
    C3 = w(8,:);
    C4 = w(9,:);
    C_s = w(10,:);
    C_i = w(11,:);
    C_age = w(12,:);
    f_max = w(13,:);
    K1 = w(14,:);
    n2 = w(15,:);
    Q_a = w(16,:);
    Q_s = w(17,:);
    T_s = w(18,:);
    sig_i = w(19,:);

    % Derived Physical Parameter
    C_e = C_s .* exp(-(Q_s./R) .* (1/TEMPERATURE - 1/T_s));
    f_e = f_max .* ((C_s - C_e) ./ C_s);
    tp = TEMPERATURE * exp(12653/TEMPERATURE - 24.475);
    tau = K1 .* tp;

    % Initial Conditions
    C_0 = C_i;
    f_0 = 0;
    r_0 = 0;

    %======================NUMERICAL INTEGRATION===========================
    % Initiallization
    time = (0:TIMESTEP:(TEMPORAL_SIZE-1)*TIMESTEP)' ;

    % rho calculation (RK4)
    rho = rungeKutta4(rho_0,C_age,n2,TIMESTEP,TEMPORAL_SIZE);

    % C calculation (since dC is calculated from rho, also RK4)
    dC = A1 .* ((C_e-C_i)./tau) .* exp(-time./tau) + B1.*rho; % rate step
    C = cumtrapz(dC,1) .* TIMESTEP + C_0;  % trapezoidal integration step


    % f calculation (since df is calculated from rho, also RK4)
    df =  - (f_e./(C_i - C_e)) .* dC;     % rate step
    f = cumtrapz(df,1) .* TIMESTEP + f_0; % trapezoidal integration step
    f(f<0) = 0;


    % r calculation
    dr = zeros(TEMPORAL_SIZE, ELEMENT_SIZE);
    dr(2:end,:) = A2.* 1/3 .* ...
                 (C1 .* (exp(-Q_a./(R*TEMPERATURE))./TEMPERATURE)).^(1/3)...
                 .*  time(2:end,:).^(-2/3) + B2 .* rho(2:end,:);
    r_corr = ...
           rCorrection(rho_0,TIMESTEP,TEMPERATURE,A2,B2,C1,Q_a,R,C_age,n2);
    r = cumtrapz(dr,1) .* TIMESTEP + r_0 + r_corr;


    % time at final step correction
    ratio = (TOTAL_TIME - time(end-1)) / (time(end) - time(end-1));
    time(end) = TOTAL_TIME;
    rho(end)  = rho(end-1) + ratio * (rho(end) - rho(end-1));
    f(end)    = f(end-1) + ratio * (f(end) - f(end-1));
    r(end)    = r(end-1) + ratio * (r(end) - r(end-1));


    % Calulate the strengthes/stresses, (TEMPORAL_SIZE,ELEMENT_SIZE)
    sig_sh  = C2 .* f.^(1/2) .* r.^(1/2);
    sig_by  = C3 .* f.^(1/2) .* r.^(-1);
    sig_ppt = (sig_by .* sig_sh) ./ (sig_by + sig_sh);
    sig_ss  = C4 .* C.^(2/3);
    sig_dis = A_d .* rho.^(0.5);

    sig_y = sig_dis + sig_ss + sig_i + sig_ppt;
    %======================================================================
end

% multi-stage ageing
function [time,sig_y] = ...
          multistage_ageing_model_predict(rho_0,...
                                          temperature_multi,...
                                          time_multi,...
                                          TEMPERATURE,...
                                          TOTAL_TIME,...
                                          TIMESTEP,...
                                          w)
    %{
    Implementation of the Multistage Ageing Model
    Arguments:
    rho_0 -- initial normalized dislocation, an array (1,ELEMENT_SIZE)
    C_0   -- initial solute fraction, an array (1,ELEMENT_SIZE)
    f_0   -- initial volume fraction, an array (1,ELEMENT_SIZE)
    r_0   -- initial mean precipitate size, an array (1,ELEMENT_SIZE)
    T     -- temperature of the ageing, a scalar
    TOTAL_TIME -- total length of the ageing process, a scalar
    tp    -- peak time, a scalar
    TIMESTEP -- numerical integration time step, a scalar
    w      -- model parameters, an array (n_w,1)

    Returns:
    sig_y
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    % Problem Dimension
    ELEMENT_SIZE = size(rho_0,2);         % m
    TEMPORAL_SIZE = ceil(TOTAL_TIME/TIMESTEP) + 1;  % n_t

    % Physical Constant
    R = 8.314;

    % Retrive Training Parameters
    A1 = w(1,:);
    A2 = w(2,:);
    A_d = w(3,:);
    B1 = w(4,:);
    B2 = w(5,:);
    C1 = w(6,:);
    C2 = w(7,:);
    C3 = w(8,:);
    C4 = w(9,:);
    C_s = w(10,:);
    C_i = w(11,:);
    C_age = w(12,:);
    f_max = w(13,:);
    K1 = w(14,:);
    n2 = w(15,:);
    Q_a = w(16,:);
    Q_s = w(17,:);
    T_s = w(18,:);
    sig_i = w(19,:);

    % Derived Physical Parameter
    C_e = C_s .* exp(-(Q_s./R) .* (1/TEMPERATURE - 1/T_s));
    f_e = f_max .* ((C_s - C_e) ./ C_s);
    tp = TEMPERATURE * exp(12653/TEMPERATURE - 24.475);
    tau = K1 .* tp;

    %======================MULTISTAGE CORRECTION===========================
    % Multistage constants
    tp_multi = temperature_multi .* exp(12653./temperature_multi - 24.475);
    tau_multi = K1 .* tp_multi;
    C_e_multi = C_s .* exp(-(Q_s./R) .* (1./temperature_multi - 1/T_s));

    % Timeshift
    t_eq_C = sum(...
                 -tau .* log(...
                            (C_e_multi - C_e)./(C_i - C_e) +...
                            (C_i - C_e)./(C_i - C_e) .* exp(-time_multi./tau_multi)...
                        )...
                 );
    t_eq_r = sum(TEMPERATURE./temperature_multi .* ...
                 exp(Q_a./(R.*TEMPERATURE) .* exp(-Q_a./(R.*temperature_multi))));
    
    
    % Shifted initial condition
    C_0 = C_e + (C_i - C_e) * exp(-t_eq_C/tau);
    r_0 = C1 .* (t_eq_r/TEMPERATURE) .* exp(-Q_a/(R*TEMPERATURE));
    f_0 = ((C_i - C_0)/(C_i - C_e)) * f_e;
    rho_0 = rho_0 + ...
            rungeKutta4_lean(rho_0,sum(time_multi),C_age,n2,TIMESTEP);

    %======================NUMERICAL INTEGRATION===========================
    % Initiallization
    time = (0:TIMESTEP:(TEMPORAL_SIZE-1)*TIMESTEP)' ;

    % rho calculation (RK4)
    rho = rungeKutta4(rho_0,C_age,n2,TIMESTEP,TEMPORAL_SIZE);

    % C calculation (since dC is calculated from rho, also RK4)
    dC = A1 .* ((C_e-C_i)./tau) .* exp(-(time+t_eq_C)./tau) + B1.*rho; % rate step
    C = cumtrapz(dC,1) .* TIMESTEP + C_0;  % trapezoidal integration step

    % f calculation (since df is calculated from rho, also RK4)
    df =  - (f_e./(C_i - C_e)) .* dC;     % rate step
    f = cumtrapz(df,1) .* TIMESTEP + f_0; % trapezoidal integration step
    f(f<0) = 0;

    % r calculation
    dr = A2.* 1/3 .* ...
         (C1 .* (exp(-Q_a./(R*TEMPERATURE))./TEMPERATURE)).^(1/3)...
         .*  (time+t_eq_r).^(-2/3) + B2 .* rho;
    r = cumtrapz(dr,1) .* TIMESTEP + r_0;


    % time at final step correction
    ratio = (TOTAL_TIME - time(end-1)) / (time(end) - time(end-1));
    time(end) = TOTAL_TIME;
    rho(end)  = rho(end-1) + ratio * (rho(end) - rho(end-1));
    f(end)    = f(end-1) + ratio * (f(end) - f(end-1));
    r(end)    = r(end-1) + ratio * (r(end) - r(end-1));


    % Calulate the strengthes/stresses, (TEMPORAL_SIZE,ELEMENT_SIZE)
    sig_sh  = C2 .* f.^(1/2) .* r.^(1/2);
    sig_by  = C3 .* f.^(1/2) .* r.^(-1);
    sig_ppt = (sig_by .* sig_sh) ./ (sig_by + sig_sh);
    sig_ss  = C4 .* C.^(2/3);
    sig_dis = A_d .* rho.^(0.5);

    sig_y = sig_dis + sig_ss + sig_i + sig_ppt;
    %======================================================================
end

%----------------------Integration Methods---------------------------------

function [x] = eulerForward(x,dx,step)
    x = dx*step + x;
end

function [rho] = rungeKutta4(rho_0,C_age,n2,TIMESTEP,TEMPORAL_SIZE)
    %{
    Implementation of RungeKutta-4 Numerical  Integration

    Arguments:
    x_0  -- initial value of x at time zero, a vector (1,ELEMENT_SIZE)
    t    -- total ageing time, a scalar
    C_age-- ageing coefficient of dislocation growth rate, a scalar
    n2   -- exponential coefficient of dislocation growth rate, a scalar
    time_step -- time step in numerical integration, a scalar


    Returns:
    x  -- value of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    dx -- derivative of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    rho  = zeros(TEMPORAL_SIZE,length(rho_0));
    rho(1,:) = rho_0;

    iteration = TEMPORAL_SIZE;
    h = TIMESTEP;

    function drho = f(rho)
        drho = - C_age .* rho .^ n2;
    end
    fun = @f;

    for i = 1:iteration-1
        k1 = fun(rho(i,:));
        x1 = rho(i,:) + k1 .* h/2;
        k2 = fun(x1);
        x2 = rho(i,:) + k2 .* h/2;
        k3 = fun(x2);
        x3 = rho(i,:) + k3 .* h;
        k4 = fun(x3);

        % update the derivative
        drho = (1/6).*h.*(k1+(2.*k2)+(2.*k3)+k4);

        % update the value
        rho(i+1,:) = rho(i,:) + drho;
    end
end

function [r_corr] = ...
              rCorrection(rho_0,TIMESTEP,TEMPERATURE,A2,B2,C1,Q_a,R,C_age,n2)
    % Gaussian quadrature points
    q1 = - sqrt(3/7 + 2/7*sqrt(6/5));
    q2 = - sqrt(3/7 - 2/7*sqrt(6/5));
    q3 =   sqrt(3/7 - 2/7*sqrt(6/5));
    q4 =   sqrt(3/7 + 2/7*sqrt(6/5));

    % Gaussian quadrature weights
    w1 = (18 - sqrt(30))/36;
    w2 = (18 + sqrt(30))/36;
    w3 = (18 + sqrt(30))/36;
    w4 = (18 - sqrt(30))/36;

    % Evaluation
    time = TIMESTEP/2 .* [-1,q1,q2,q3,q4] + TIMESTEP/2;
    timestep = diff(time);
    rho  = rungeKutta4_variableTimestep(rho_0,timestep,C_age,n2);

    % Calculate r_corr
    g = A2.* 1/3.*(C1 .* (exp(-Q_a./(R*TEMPERATURE))./TEMPERATURE)).^(1/3);
    r_corr = g/3 * TIMESTEP^(1/3) + B2 .* sum([0;w1;w2;w3;w4] .* rho, 1);
end

function [rho] = rungeKutta4_variableTimestep(rho_0,timestep,C_age,n2)
    %{
    Implementation of RungeKutta-4 Numerical  Integration

    Arguments:
    x_0  -- initial value of x at time zero, a vector (1,ELEMENT_SIZE)
    t    -- total ageing time, a scalar
    C_age-- ageing coefficient of dislocation growth rate, a scalar
    n2   -- exponential coefficient of dislocation growth rate, a scalar
    time_step -- time step in numerical integration, a scalar


    Returns:
    x  -- value of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    dx -- derivative of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    rho  = zeros(length(timestep)+1,length(rho_0));
    rho(1,:) = rho_0;

    iteration = length(timestep);
    h = timestep;

    function drho = f(rho)
        drho = - C_age .* rho .^ n2;
    end
    fun = @f;

    for i = 1:iteration
        k1 = fun(rho(i,:));
        x1 = rho(i,:) + k1 .* h(i)/2;
        k2 = fun(x1);
        x2 = rho(i,:) + k2 .* h(i)/2;
        k3 = fun(x2);
        x3 = rho(i,:) + k3 .* h(i);
        k4 = fun(x3);

        % update the derivative
        drho = (1/6).*h(i).*(k1+(2.*k2)+(2.*k3)+k4);

        % update the value
        rho(i+1,:) = rho(i,:) + drho;
    end
end

function [rho] = rungeKutta4_lean(rho_0,t,C_age,n2,TIMESTEP)
    %{
    Lean Implementation of RungeKutta-4 Numerical Integration

    Arguments:
    x_0  -- initial value of x at time zero, a vector (1,ELEMENT_SIZE)
    t    -- total ageing time, a scalar
    C_age-- ageing coefficient of dislocation growth rate, a scalar
    n2   -- exponential coefficient of dislocation growth rate, a scalar
    time_step -- time step in numerical integration, a scalar


    Returns:
    x  -- value of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    dx -- derivative of x, an array (TEMPORAL_SIZE, ELEMENT_SIZE)
    %**********************************************************************
    Author: https://github.com/xinyangyuan
    %**********************************************************************
    %}

    rho = rho_0;
    rho_pre = zeros(size(rho));

    iteration = ceil(t/TIMESTEP);
    h = TIMESTEP;

    function drho = f(rho)
        drho = - C_age .* rho .^ n2;
    end
    fun = @f;

    for i = 1:iteration
        k1 = fun(rho);
        x1 = rho + k1 .* h/2;
        k2 = fun(x1);
        x2 = rho + k2 .* h/2;
        k3 = fun(x2);
        x3 = rho + k3 .* h;
        k4 = fun(x3);

        % update the derivative
        drho = (1/6).*h.*(k1+(2.*k2)+(2.*k3)+k4);

        % update value
        rho = rho + drho;

        if i == (iteration-1)
            rho_pre = rho;
        end

    end

    % time at final step correction
    ratio = (t - (iteration-1)*TIMESTEP) ...
            /(iteration*TIMESTEP - (iteration-1)*TIMESTEP);
    rho  = rho_pre + ratio .* (rho - rho_pre);

end
