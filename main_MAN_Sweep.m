%% General QNM Computation


clear

addpath("SRC\")
addpath("SRC\maxdistcolor\")

solver = 'eig'; 
% Change name
name = 'QNMEig_const_Sphere.mph#'; 
% Change Modes
N_modes = 10; 

wl_min = 425e-9;
wl_max = 455e-9;

% Set Calculation Boundaries
W.lambda_list = linspace(wl_min,wl_max,201);
omegas = 2*pi./W.lambda_list.*retconstantes("c"); 




folder0='RESULTS';
    if ~exist(folder0, 'dir')
        mkdir(folder0);
    end
 


RESULTS.sweep = linspace(50,100,3);

RESULTS.Models = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.general_solutions = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.resonator_solutions = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.fields.YZ.E = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.fields.YZ.mesh = repmat(struct('Y',[],'Z',[]), 1, length(RESULTS.sweep));
RESULTS.fields.XZ.E = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.fields.XZ.mesh = repmat(struct('X',[], 'Z', []), 1, length(RESULTS.sweep));
RESULTS.fields.XY.E = repmat(struct('data',[]), 1, length(RESULTS.sweep));
RESULTS.fields.XY.mesh = repmat(struct('X',[], 'Y', []), 1, length(RESULTS.sweep));


for ir=1:length(RESULTS.sweep)

    folder=num2str(RESULTS.sweep(ir));
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    movefile(folder,['RESULTS/', folder]);
    
    Model = loadModels(name,'solver',solver, 'modes',N_modes); 
    Model.comsol_model.param.set('r0', strcat(num2str(RESULTS.sweep(ir)),'[nm]'), 'radius');
    

    % Particle name LIKE SELECTION in COMSOL
    scatter_label_pattern = "Sphere"; 
    scatter_selections = getSelections(Model, scatter_label_pattern, ... 
        'inlabel', true); 
    % For Dispersive Materials
    dispersive_materials = getDispersiveMaterial(Model); 
     
    materials = getMaterialType(Model, scatter_selections, ... 
        dispersive_materials); 
    % type in wl cycle: epsw_sca = getEpsilons(omega,dispersive_materials{1});
    
    % For constant-Dispersive Materials
    epsw_sca = 36;
    
    RESULTS.Models(ir).data = Model;

    % Start EIGENSOLVER
    QNM = getEigQNMs(Model);
    RESULTS.general_solutions(ir).data = QNM;
    
   
    % Get Fields
    QNM_resonator = getModelResonatorFields(Model, QNM, scatter_selections);
    RESULTS.resonator_solutions(ir).data = QNM_resonator;

    Q_factors = real(QNM.omega) ./ (2 * imag(QNM.omega));

    % Save Eigenfrequencies and Q-factors
    Matrix_EFQ = [...
        linspace(1,length(QNM_resonator{1}.omega),length(QNM_resonator{1}.omega))',...
        real(QNM_resonator{1}.omega), imag(QNM_resonator{1}.omega), Q_factors];
    writematrix(Matrix_EFQ, fullfile(folder0,folder,strcat("EigenMode_",num2str(RESULTS.sweep(ir)),"nm_Data.txt")),'Delimiter','space');
    

    
    % YZ plane in x=0 cut. Anoter ones in analogy.
    y = linspace(-350e-9,350e-9,101);
    z = linspace(-350e-9,350e-9,101);
    [Y, Z] = meshgrid(y,z);
    
    % Change cut-plane
    x = 0e-9; 
    
    % Making mesh
    coord = [x*ones(1,numel(Y));reshape(Y,1,[]);reshape(Z,1,[])];
    RESULTS.fields.YZ.mesh(ir).Y = Y;
    RESULTS.fields.YZ.mesh(ir).Z = Z;
    RESULTS.fields.YZ.E(ir).data = getModalFields(Model, QNM, 'coord',coord);
    

    % XZ plane in x=0 cut. Anoter ones in analogy.
    x = linspace(-350e-9,350e-9,101);
    z = linspace(-350e-9,350e-9,101);
    [X, Z] = meshgrid(x,z);
    
    % Change cut-plane
    y = 0e-9; 
    
    % Making mesh
    coord = [reshape(X,1,[]);y*ones(1,numel(X));reshape(Z,1,[])];
    RESULTS.fields.XZ.mesh(ir).X = X;
    RESULTS.fields.XZ.mesh(ir).Z = Z;
    RESULTS.fields.XZ.E(ir).data = getModalFields(Model, QNM, 'coord',coord);

    % XY plane in x=0 cut. Anoter ones in analogy.
    x = linspace(-350e-9,350e-9,101);
    y = linspace(-350e-9,350e-9,101);
    [X, Y] = meshgrid(x,y);
    
    % Change cut-plane
    z = 0e-9; 
    
    % Making mesh
    coord = [reshape(X,1,[]);reshape(Y,1,[]);z*ones(1,numel(X))];
    RESULTS.fields.XY.mesh(ir).X = X;
    RESULTS.fields.XY.mesh(ir).Y = Y;
    RESULTS.fields.XY.E(ir).data = getModalFields(Model, QNM, 'coord',coord);

end
%% Modal Fields

for ir=1:length(RESULTS.sweep)
    folder=num2str(RESULTS.sweep(ir));
    for modeId=1:N_modes
        
            QNM_fields = RESULTS.fields.YZ.E(ir).data;
            Y = RESULTS.fields.YZ.mesh(ir).Y;
            Z = RESULTS.fields.YZ.mesh(ir).Z;
    
            % Reshape the array (norm E -> other ones in analogy)
            Ex = reshape(QNM_fields.Ex(modeId,:),size(Y));
            Ey = reshape(QNM_fields.Ey(modeId,:),size(Y));
            Ez = reshape(QNM_fields.Ez(modeId,:),size(Y));
        
            E = abs(Ex.^2 + Ey.^2 + Ez.^2);
        
            % Plot the norm field
            figure('visible','off')
            surf(Y,Z,E,'LineStyle','none');
            colorbar;view(2);
            title("YZ norm E of mode#" + modeId);
            xlabel('y, m');
            ylabel('z, m');
            colormap(hot);
        
            folder1 = 'Fields YZ plane';
            if ~exist(folder1, 'dir')
                mkdir(folder1);
            end
            
        
            % For Save .jpg
            
            ax = gca;
            exportgraphics(ax,fullfile(folder1,['Mode',num2str(RESULTS.sweep(ir)),'_',num2str(modeId),'.jpg']),'Resolution',300);
            
    end 
    movefile(folder1,[strcat(folder0,'/',folder,'/'), folder1]);
end
%% Extinction  Calculation for every mode

%--- 

% Based on DOI: 10.1002/lpor.201700113
% Formulas 4.7, 5.5 and 9.3

%---

for ir=1:length(RESULTS.sweep)
    
    folder=num2str(RESULTS.sweep(ir));
    QNM_resonator = RESULTS.resonator_solutions(ir).data;

    normalization = (RESULTS.sweep(ir)*1e-9)^2 * pi;

    % Get Mesh
    X = QNM_resonator{1}.coord(1,:);
    Y = QNM_resonator{1}.coord(2,:);
    Z = QNM_resonator{1}.coord(3,:);
    
    % Poynting Normalization
    W.S0=0.5*sqrt(retconstantes("eps0")/retconstantes("mu0")); 
    
     
    ext_m = zeros(length(W.lambda_list), N_modes);
    abs_m = zeros(length(W.lambda_list), N_modes); 
    scs_m = zeros(length(W.lambda_list), N_modes); 
    
    % Start Calculations
    for ilambda = 1:length(W.lambda_list) 
        lambda = W.lambda_list(ilambda); 
        omega = omegas(ilambda); 
        
        % Set Incl. Wave Params
        k0 = 2*pi/lambda; 
        E0 = 1; 
         
        kx=k0;  
        ky=0; 
        kz=0; 
     
        phase = exp(1i*kx*X+1i*ky*Y + 1i*kz*Z); 
     
        Exincl=E0*phase*(0); 
        Eyincl=E0*phase*(1); 
        Ezincl=0.*phase; 
     
     
        fac0 = conj(epsw_sca-1)*retconstantes("eps0")*... 
                            omega/2/W.S0; 
        fac1 = (epsw_sca-1) * retconstantes("eps0"); 
    
        % Go through every mode
        for isol=1:N_modes
    
            % Get Normalized Mode
            Emx = QNM_resonator{1}.Ex(isol,:); 
            Emy = QNM_resonator{1}.Ey(isol,:); 
            Emz = QNM_resonator{1}.Ez(isol,:); 
            
    
            % Calculating extinced coefficient
            coef_w = omega/(QNM.omega(isol)-omega); 
            alpha_m = sum( ... 
                fac1*coef_w.*QNM_resonator{1}.mesh_vol.*( ... 
                Emx.*Exincl + Emy.*Eyincl + Emz.*Ezincl)); 
            % Calculating sigma ext
            ext_m(ilambda, isol)= imag( ... 
                sum((conj(Emx).*Exincl+conj(Emy).*Eyincl + ... 
                conj(Emz).*Ezincl).*QNM_resonator{1}.mesh_vol.*fac0)*... 
                conj(alpha_m)) + ext_m(ilambda, isol); 
            
            Es = transpose([(alpha_m .* Emx)',...
                (alpha_m .* Emy)', (alpha_m .* Emz)']);
            Eb = transpose([Exincl', Eyincl', Ezincl']);
            % For k-partice use carefully
            abs_m(ilambda, isol)= sum(...
                imag(epsw_sca) .* sum((Es+Eb) .* conj(Es+Eb)) ...
                .*QNM_resonator{1}.mesh_vol.*retconstantes("eps0")*... 
                omega/2/W.S0) + abs_m(ilambda, isol);
    
            scs_m(ilambda, isol)= imag(...
                sum((conj(Emx).*(alpha_m.*Emx+Exincl)+...
                conj(Emy).*(alpha_m.*Emy+Eyincl) + ... 
                conj(Emz).*(alpha_m.*Emz+Ezincl))...
                .*QNM_resonator{1}.mesh_vol.*(fac0))*... 
                conj(alpha_m) + scs_m(ilambda, isol));
            
    
        end 
    end 
    % Normalization in geometry shadow
    ext_m = ext_m/normalization;
    abs_m = abs_m/normalization;
    scs_m = scs_m/normalization;
    
    rgb = maxdistcolor(N_modes,@sRGB_to_OKLab);
    
    %Plotting SCS: Modes + Summary Spectrum
    figure();
    
    for isol=1:N_modes
        plot(W.lambda_list.*1e9, ext_m(:,isol), "DisplayName", ...
            [strcat("Mode #",num2str(isol))],'Color',rgb(isol,:));
        hold on;
        scatter(2*pi*retconstantes("c")./real(QNM_resonator{1}.omega(isol))*1e9,...
            0,"DisplayName", ...
            [strcat("Mode #",num2str(isol))], 'MarkerEdgeColor',rgb(isol,:));
        
    end
    
    scs_sum = sum(scs_m(:,:),2);
    ext_sum = sum(ext_m(:,:),2);
    abs_sum = sum(abs_m(:,:),2);
    plot(W.lambda_list.*1e9, ext_sum, "DisplayName", "SUM");
    
    xlabel('wavelength');
    ylabel("\sigma_{SCS}");
    legend();
    
    
    folder2='EXTIN';
    if ~exist(folder2, 'dir')
        mkdir(folder2);
    end
    
    
    folder3='SCA';
    if ~exist(folder3, 'dir')
        mkdir(folder3);
    end
    
    folder4='ABS';
    if ~exist(folder4, 'dir')
        mkdir(folder4);
    end
    
    
    Matrix_ext = [W.lambda_list', ext_m, ext_sum];
    writematrix(Matrix_ext, fullfile(folder2,strcat('Mode_excinction_Data.txt')),'Delimiter','space');
    
    
    Matrix_scs = [W.lambda_list', scs_m, scs_sum];
    writematrix(Matrix_scs,fullfile(folder3,strcat('Mode_scattering_Data.txt')),'Delimiter','space');
    
    
    Matrix_abs = [W.lambda_list', abs_m, abs_sum];
    writematrix(Matrix_abs,fullfile(folder4,strcat('Mode_absorption_Data.txt')),'Delimiter','space');
    
    movefile(folder2,[strcat(folder0,'/',folder,'/'), folder2]);
    movefile(folder3,[strcat(folder0,'/',folder,'/'), folder3]);
    movefile(folder4,[strcat(folder0,'/',folder,'/'), folder4]);
end
