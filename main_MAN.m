
% if you want to make parametric sweep use this:
% make a variable and for-cycle with repeat of next sections what you want
% ... success!

% Set parameters in Global Definition: name, value, comment
% Model.comsol_model.param.set('r0', '50[nm]', 'radius');


%% General QNM Computation

clear
solver = 'eig'; 
% Change name
name = 'QNMEig_const_Sphere.mph#'; 
% Change Modes
N_modes = 100; 

wl_min = 400e-9;
wl_max = 1000e-9;
 
addpath("SRC\") 
 
Model = loadModels(name,'solver',solver, 'modes',N_modes); 




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
epsw_sca = 16; 


% Set Calculation Boundaries
W.lambda_list = linspace(wl_min,wl_max,201);
omegas = 2*pi./W.lambda_list.*retconstantes("c"); 


% Geometry Shadow Cross-Section Normalization
normalization = 100e-9^2 * pi;


% Start EIGENSOLVER
QNM = getEigQNMs(Model); 

% Plot Modes in Re-Im plane with Q-factor and other parameters
plotComplexPlane(Model,QNM);
Q_factors = real(QNM.omega) ./ (2 * imag(QNM.omega));

% Get Fields
QNM_resonator = getModelResonatorFields(Model, QNM, scatter_selections);



%% Modal Fields

for modeId=1:N_modes

    % YZ plane in x=0 cut. Anoter ones in analogy.
    y = linspace(-350e-9,350e-9,101);
    z = linspace(-350e-9,350e-9,101);
    [Y, Z] = meshgrid(y,z);
    
    % Change cut-plane
    x = 0e-9; 

    % Making mesh
    coord = [x*ones(1,numel(Y));reshape(Y,1,[]);reshape(Z,1,[])];
    
    QNM_fields = getModalFields(Model, QNM, 'coord',coord);
    
    % Reshape the array (norm E -> other ones in analogy)
    Ex = reshape(QNM_fields.Ex(modeId,:),size(Y));
    Ey = reshape(QNM_fields.Ey(modeId,:),size(Y));
    Ez = reshape(QNM_fields.Ez(modeId,:),size(Y));

    E = abs(Ex.^2 + Ey.^2 + Ez.^2);

    % Plot the norm field
    figure();
    surf(Y,Z,E,'LineStyle','none');
    colorbar;view(2);title("norm E of mode#" + modeId);
    
    % For Save .jpg
    
    %ax = gca;
    %exportgraphics(ax,num2str(modeId),'Resolution',300)

end 

%% Extinction  Calculation for every mode

%--- 

% Based on DOI: 10.1002/lpor.201700113
% Formulas 4.7, 5.5 and 9.3

%---

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
    for isol=1:N_modes-1

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

%colors = ["red", "green", "blue", "black", "yellow", "magenta", "#112333", "#DDA123", "#FF1232", "#ABCF88"];

%Plotting SCS: Modes + Summary Spectrum
figure;
for isol=1:N_modes-1
    plot(W.lambda_list.*1e9, scs_m(:,isol), "DisplayName", ...
        [strcat("Mode #",num2str(isol))]);
    hold on;
    scatter(2*pi*retconstantes("c")./real(QNM_resonator{1}.omega(isol))*1e9,...
        0,"DisplayName", ...
        [strcat("Mode #",num2str(isol))]);
    
end

scs_sum = sum(scs_m(:,:),2);
ext_sum = sum(ext_m(:,:),2);
plot(W.lambda_list.*1e9, scs_sum, "DisplayName", "SUM");
plot(W.lambda_list.*1e9, ext_sum, "DisplayName", "EXT");
xlabel('wavelength');
ylabel("\sigma_{EXT}");
legend();




%% Multipole Decomposition of full field



W.lambda_list_MD = linspace(wl_min, wl_max, length(X));
omegas = 2*pi./W.lambda_list_MD.*retconstantes("c"); 
k_list = 2 * pi ./ W.lambda_list_MD;

MULTIPOLES = zeros(length(W.lambda_list_MD), 4);

% Calculating Scattering Field
for ilambda = 1:length(W.lambda_list_MD) 
    lambda = W.lambda_list_MD(ilambda); 
    omega = omegas(ilambda); 
    
    % Set Incl. Wave Params
    k0 = k_list(ilambda); 
    E0 = 1; 
     
    kx=k0;  
    ky=0; 
    kz=0; 
 
    phase = exp(1i*kx*X + 1i*ky*Y + 1i*kz*Z); 
 
    Exincl=E0*phase*(0); 
    Eyincl=E0*phase*(1); 
    Ezincl=0.*phase; 
    
    fac2 = (epsw_sca-1) * retconstantes("eps0"); 
    fac3 = -1i * omega;
    
    Eb = transpose([Exincl', Eyincl', Ezincl']);
    Es = zeros(3,length(Exincl));
    %zeros(3,length(Exincl));
    for isol=1:N_modes-1
        Emx = QNM_resonator{1}.Ex(isol,:); 
        Emy = QNM_resonator{1}.Ey(isol,:); 
        Emz = QNM_resonator{1}.Ez(isol,:); 

 
        coef_w = omega/(QNM.omega(isol)-omega); 
        alpha_m = sum( ... 
            fac2*coef_w.*QNM_resonator{1}.mesh_vol.*( ... 
            Emx.*Exincl + Emy.*Eyincl + Emz.*Ezincl));

        Es(1,:) = Es(1,:) + alpha_m .* Emx;
        Es(2,:) = Es(2,:) + alpha_m .* Emy;
        Es(3,:) = Es(3,:) + alpha_m .* Emz;
    end

    J = fac2 .* fac3 .* conj(Es);
    R = transpose([X', Y', Z']);

    % Additional variables:
    r2 = X.*X + Y.*Y + Z.*Z;
    kr  = k0 .* sqrt(r2);
        

    j0 = sqrt(pi./(2*kr)).*besselj(0+1/2,kr);  
    j1 = sqrt(pi./(2*kr)).*besselj(1+1/2,kr);  
    j2 = sqrt(pi./(2*kr)).*besselj(2+1/2,kr);  
    j3 = sqrt(pi./(2*kr)).*besselj(3+1/2,kr);  

    RJ = R(1,:).*J(1,:) + R(2,:).*J(2,:) +  R(3,:).*J(3,:);
    RxJ = transpose([(R(2,:).*J(3,:) - R(3,:).*J(2,:))',...
        (R(3,:).*J(1,:) - R(1,:).*J(3,:))',...
        (R(1,:).*J(2,:) - R(2,:).*J(1,:))']);

    const = k0^4./(6 * pi * retconstantes("eps0")^2 * abs(E0)^2);
    
    % ED
    p = zeros(1,3);
    for ic=1:3

        % Momentum
        % p(ic) = -1/(1i*omega) * (sum(J(ic,:).*j0 ...
        %     .*QNM_resonator{1}.mesh_vol) + ...
        %     (k0^2)/2 * sum((3.*RJ.*R(ic,:) - r2 .* J(ic,:)) .*...
        %      j2./ (kr.^2) .* QNM_resonator{1}.mesh_vol));
        % 
        % p(ic) = -1/(1i .* omega) .* ((k0^2)/2 .* sum(...
        %     (3.*RJ.*R(ic,:) - r2.*J(ic,:)) .* (j2 ./ (kr.^2)) ...
        %     .* QNM_resonator{1}.mesh_vol) + sum(...
        %     J(ic,:).*j0.*QNM_resonator{1}.mesh_vol));

        p(ic) = -1/(1i*omega) * sum(J(ic,:) ...
            .*QNM_resonator{1}.mesh_vol);


        % SCS
        MULTIPOLES(ilambda, 1) = MULTIPOLES(ilambda, 1) ...
            +  const .* (abs(p(ic)).^2);
    end
       

    % MD
    m = zeros(1,3);
    for ic=1:3

        % Momentum
        m(ic) = 3/2 * sum((RxJ(ic,:)).*...
            (j1./kr).*QNM_resonator{1}.mesh_vol);
            
        % SCS
        MULTIPOLES(ilambda, 2) = MULTIPOLES(ilambda, 2) ...
            + const * ...
            (retconstantes("mu0") * retconstantes("eps0")).*...
            (m(ic) .* conj(m(ic))) ;
    end


    % EQ
    EQ = zeros(3,3);

        
    for ic0=1:3
        for ic1=1:3
                
         % Kroneker  delta Definition
            if ic0==ic1
                kronek = 1;
            else
                kronek = 0;
            end

            % Momentum
            EQ(ic0, ic1) = -3/(1i*omega) * (...
                sum((3*(R(ic1,:).*J(ic0,:) + R(ic0,:).*J(ic1,:))...
                - 2 .* RJ .* kronek ) .* j1 ./ kr ...
                .* QNM_resonator{1}.mesh_vol) + ...
                2 .*(k0^2) .* sum((5.*R(ic0,:).*R(ic1,:).*RJ ...
                - (R(ic0,:).*J(ic1,:) + R(ic1,:).*J(ic0,:)) .* r2...
                - r2 .* (RJ) .* kronek) .*j3./(kr.^3) ...
                .* QNM_resonator{1}.mesh_vol));

           % SCS
           MULTIPOLES(ilambda, 3) = ...
               MULTIPOLES(ilambda, 3) + ...
               const ./ 120 .* (k0.^2) .* (EQ(ic0, ic1) .* ...
               conj(EQ(ic0, ic1)));
        end
    end
       

    % MQ
    MQ = zeros(3,3);

    for ic0=1:3
        for ic1=1:3
                
        % Momentum
            MQ(ic0, ic1) = 15 * (...
                sum((R(ic0,:).*RxJ(ic1,:) ...
                +R(ic1,:).*RxJ(ic0,:)) ...
                .* j2./(kr.^2) ...
                .* QNM_resonator{1}.mesh_vol));
               
           % SCS
           MULTIPOLES(ilambda, 4) = ...
               MULTIPOLES(ilambda, 4) + ...
               const *retconstantes("mu0") * retconstantes("eps0") ./ (120) .* ...
               (k0.^2) .* (MQ(ic0, ic1) .* conj(MQ(ic0, ic1)));
        end
    end

end
multipole_names = ["ED" "MD" "EQ" "MQ"];
figure;
for i=1:4
    plot(W.lambda_list_MD.*1e9, transpose(MULTIPOLES(:,i))...
        ./normalization, 'DisplayName',multipole_names(i));
    hold on;
end
% ./(W.lambda_list.^2/(2*pi))
TOTAL_MD = sum(MULTIPOLES(:,:),2);
 
plot(W.lambda_list_MD.*1e9, transpose(TOTAL_MD)...
    ./normalization,'DisplayName', 'SUM MD');

plot(W.lambda_list.*1e9, scs_sum, "DisplayName", "scs");
xlabel('wavelength');
ylabel("\sigma_{SCS}");
legend();


%% MD TEST

c = retconstantes("c");
eps0 = retconstantes("eps0");
W.lambda_list = linspace(wl_min, wl_max, length(X));
omegas = 2*pi./W.lambda_list.*retconstantes("c"); 
k_list = 2 * pi ./ W.lambda_list;

MULTIPOLES = zeros(length(W.lambda_list), 4);

% Calculating Scattering Field
for ilambda = 1:length(W.lambda_list) 
    lambda = W.lambda_list(ilambda); 
    omega = omegas(ilambda); 
    
    % Set Incl. Wave Params
    k0 = k_list(ilambda); 
    E0 = 1; 
     
    kx=k0;  
    ky=0; 
    kz=0; 
 
    phase = exp(1i*kx*X + 1i*ky*Y + 1i*kz*Z); 
 
    Exincl=E0*phase*(0); 
    Eyincl=E0*phase*(1); 
    Ezincl=0.*phase; 
    
    fac2 = (epsw_sca-1) * retconstantes("eps0"); 
    fac3 = 1i * omega;
    
    Eb = transpose([Exincl', Eyincl', Ezincl']);
    Es = zeros(3,length(Exincl));
    
    for isol=1:N_modes
        Emx = QNM_resonator{1}.Ex(isol,:); 
        Emy = QNM_resonator{1}.Ey(isol,:); 
        Emz = QNM_resonator{1}.Ez(isol,:); 

 
        coef_w = omega/(QNM.omega(isol)-omega); 
        alpha_m = sum( ... 
            fac2*coef_w.*QNM_resonator{1}.mesh_vol.*( ... 
            Emx.*Exincl + Emy.*Eyincl + Emz.*Ezincl));

        Es(1,:) = Es(1,:) + alpha_m .* Emx;
        Es(2,:) = Es(2,:) + alpha_m .* Emy;
        Es(3,:) = Es(3,:) + alpha_m .* Emz;
    end

    J = fac2 .* fac3 .* (Es);
    Jx = J(1,:);
    Jy = J(2,:);
    Jz = J(3,:);

    x4d = X;
    y4d = Y;
    z4d = Z;
    k4d = k_list;
    k = k0;
    d3r = QNM_resonator{1}.mesh_vol;

    % calculate often used values
    % constant for scattering cross section
    const = k.^4/(6*pi*eps0^2*1);  % E0 = 1
    
    % scalar product
    rJ = x4d.*Jx + y4d.*Jy + z4d.*Jz;  % product r,J(r)
    rr = x4d.*x4d + y4d.*y4d + z4d.*z4d;  % product r,r
    r = sqrt(rr);  % norm(r)
    
    % cross product r x J = (ry*Jz-rz*Jy, rz*Jx-rx*Jz, rx*Jy-ry*Jx)
    rxJx = (y4d.*Jz - z4d.*Jy);
    rxJy = (z4d.*Jx - x4d.*Jz);
    rxJz = (x4d.*Jy - y4d.*Jx);
    
    % spherical bessel functions
    sbj0 = sqrt(pi./(2*k4d.*r)).*besselj(0+1/2,k4d.*r);  %j0
    sbj1 = sqrt(pi./(2*k4d.*r)).*besselj(1+1/2,k4d.*r);  %j1
    sbj2 = sqrt(pi./(2*k4d.*r)).*besselj(2+1/2,k4d.*r);  %j2
    sbj3 = sqrt(pi./(2*k4d.*r)).*besselj(3+1/2,k4d.*r);  %j3
    
    % calculate multipole moments and cross sections
    % calculate electric dipole p,
    dpx = (3*rJ.*x4d-rr.*Jx).*sbj2./(k4d.*r).^2;
    dpy = (3*rJ.*y4d-rr.*Jy).*sbj2./(k4d.*r).^2;
    dpz = (3*rJ.*z4d-rr.*Jz).*sbj2./(k4d.*r).^2;
    px = -1./(1i*omega).*(sum(Jx.*sbj0.*d3r)+ ...
    k.^2/2.*sum(dpx.*d3r));
    py = -1./(1i*omega).*(sum(Jy.*sbj0.*d3r)+ ...
    k.^2/2.*sum(dpy.*d3r));
    pz = -1./(1i*omega).*(sum(Jz.*sbj0.*d3r)+ ...
    k.^2/2.*sum(dpz.*d3r));
    norm2_p = px.*conj(px)+py.*conj(py)+pz.*conj(pz);
    Cp = const.*norm2_p;
    MULTIPOLES(ilambda, 1) = Cp;
    
    % calculate magnetic dipole m,
    dmx = rxJx.*sbj1./(k4d.*r);
    dmy = rxJy.*sbj1./(k4d.*r);
    dmz = rxJz.*sbj1./(k4d.*r);
    mx = 3/2*sum(dmx.*d3r);
    my = 3/2*sum(dmy.*d3r);
    mz = 3/2*sum(dmz.*d3r);
    norm2_m = mx.*conj(mx)+my.*conj(my)+mz.*conj(mz);
    Cm = const.*norm2_m/c^2;
    MULTIPOLES(ilambda, 2) = Cm;
    
    % calculate electric quadrupole Qe
    dQe1xx = (3*2*x4d.*Jx - 2*rJ).*sbj1./(k4d.*r);
    dQe1xy = (3*(y4d.*Jx+x4d.*Jy)).*sbj1./(k4d.*r);
    dQe1xz = (3*(z4d.*Jx+x4d.*Jz)).*sbj1./(k4d.*r);
    dQe1yy = (3*2*y4d.*Jy - 2.*rJ).*sbj1./(k4d.*r);
    dQe1yx = (3*(x4d.*Jy+y4d.*Jx)).*sbj1./(k4d.*r);
    dQe1yz = (3*(z4d.*Jy+y4d.*Jz)).*sbj1./(k4d.*r);
    dQe1zz = (3*2*z4d.*Jz - 2*rJ).*sbj1./(k4d.*r);
    dQe1zx = (3*(x4d.*Jz+z4d.*Jx)).*sbj1./(k4d.*r);
    dQe1zy = (3*(y4d.*Jz+z4d.*Jy)).*sbj1./(k4d.*r);
    dQe2xx = (5*x4d.*x4d.*rJ - rr*2.*x4d.*Jx - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2xy = (5*x4d.*y4d.*rJ - rr.*(x4d.*Jy+y4d.*Jx)).*sbj3./(k4d.*r).^3;
    dQe2xz = (5*x4d.*z4d.*rJ - rr.*(x4d.*Jz+z4d.*Jx)).*sbj3./(k4d.*r).^3;
    dQe2yy = (5*y4d.*y4d.*rJ - rr*2.*y4d.*Jy - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2yx = (5*y4d.*x4d.*rJ - rr.*(y4d.*Jx+x4d.*Jy)).*sbj3./(k4d.*r).^3;
    dQe2yz = (5*y4d.*z4d.*rJ - rr.*(y4d.*Jz+z4d.*Jy)).*sbj3./(k4d.*r).^3;
    dQe2zz = (5*z4d.*z4d.*rJ - rr*2.*z4d.*Jz - rr.*rJ).*sbj3./(k4d.*r).^3;
    dQe2zx = (5*z4d.*x4d.*rJ - rr.*(z4d.*Jx+x4d.*Jz)).*sbj3./(k4d.*r).^3;
    dQe2zy = (5*z4d.*y4d.*rJ - rr.*(z4d.*Jy+y4d.*Jz)).*sbj3./(k4d.*r).^3;
    Qexx = -3./(1i*omega).*(sum(dQe1xx.*d3r)+ ...
    2*k.^2.*sum(dQe2xx.*d3r));
    Qexy = -3./(1i*omega).*(sum(dQe1xy.*d3r)+ ...
    2*k.^2.*sum(dQe2xy.*d3r));
    Qexz = -3./(1i*omega).*(sum(dQe1xz.*d3r)+ ...
    2*k.^2.*sum(dQe2xz.*d3r));
    Qeyy = -3./(1i*omega).*(sum(dQe1yy.*d3r)+ ...
    2*k.^2.*sum(dQe2yy.*d3r));
    Qeyx = -3./(1i*omega).*(sum(dQe1yx.*d3r)+ ...
    2*k.^2.*sum(dQe2yx.*d3r));
    Qeyz = -3./(1i*omega).*(sum(dQe1yz.*d3r)+ ...
    2*k.^2.*sum(dQe2yz.*d3r));
    Qezz = -3./(1i*omega).*(sum(dQe1zz.*d3r)+ ...
    2*k.^2.*sum(dQe2zz.*d3r));
    Qezx = -3./(1i*omega).*(sum(dQe1zx.*d3r)+ ...
    2*k.^2.*sum(dQe2zx.*d3r));
    Qezy = -3./(1i*omega).*(sum(dQe1zy.*d3r)+ ...
    2*k.^2.*sum(dQe2zy.*d3r));
    norm2_Qe = Qexx.*conj(Qexx)+Qexy.*conj(Qexy)+Qexz.*conj(Qexz)+ ...
       Qeyy.*conj(Qeyy)+Qeyx.*conj(Qeyx)+Qeyz.*conj(Qeyz)+ ...
       Qezz.*conj(Qezz)+Qezx.*conj(Qezx)+Qezy.*conj(Qezy);
    CQe = const/120.*k.^2.*norm2_Qe;
    MULTIPOLES(ilambda, 3) = CQe;
    
    % calculate magnetic quadrupole Qm
    dQmxx = (2*x4d.*rxJx).*sbj2./(k4d.*r).^2;
    dQmxy = (x4d.*rxJy+y4d.*rxJx).*sbj2./(k4d.*r).^2;
    dQmxz = (x4d.*rxJz+x4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmyy = (2*y4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmyx = (y4d.*rxJx+x4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmyz = (y4d.*rxJz+z4d.*rxJy).*sbj2./(k4d.*r).^2;
    dQmzz = (2*z4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmzx = (z4d.*rxJx+x4d.*rxJz).*sbj2./(k4d.*r).^2;
    dQmzy = (z4d.*rxJy+y4d.*rxJz).*sbj2./(k4d.*r).^2;
    Qmxx = 15*sum(dQmxx.*d3r);
    Qmxy = 15*sum(dQmxy.*d3r);
    Qmxz = 15*sum(dQmxz.*d3r);
    Qmyy = 15*sum(dQmyy.*d3r);
    Qmyx = 15*sum(dQmyx.*d3r);
    Qmyz = 15*sum(dQmyz.*d3r);
    Qmzz = 15*sum(dQmzz.*d3r);
    Qmzx = 15*sum(dQmzx.*d3r);
    Qmzy = 15*sum(dQmzy.*d3r);
    norm2_Qm = Qmxx.*conj(Qmxx)+Qmxy.*conj(Qmxy)+Qmxz.*conj(Qmxz)+ ...
       Qmyy.*conj(Qmyy)+Qmyx.*conj(Qmyx)+Qmyz.*conj(Qmyz)+ ...
       Qmzz.*conj(Qmzz)+Qmzx.*conj(Qmzx)+Qmzy.*conj(Qmzy);
    CQm = const./120.*(k/c).^2.*norm2_Qm;
    MULTIPOLES(ilambda, 4) = CQm;
   
end
multipole_names = ["ED" "MD" "EQ" "MQ"];
figure;
for i=1:4
    plot(W.lambda_list.*1e9, transpose(MULTIPOLES(:,i))...
        ./normalization, 'DisplayName', multipole_names(i));
    hold on;
end
% ./(W.lambda_list.^2/(2*pi))
TOTAL_MD = sum(MULTIPOLES(:,:),2);
 
plot(W.lambda_list.*1e9, transpose(TOTAL_MD)...
    ./normalization,'DisplayName', 'SUM MD');

xlabel('wavelength');
ylabel("\sigma_{SCS}");
legend();

plot(linspace(425e-9,455e-9,201).*1e9, scs_sum, "DisplayName", "scs");
xlabel('wavelength');
ylabel("\sigma_{SCS}");
legend();


%% Multipole Decomposition of Every Mode

%---

% Based on DOI: 10.1016/j.optcom.2017.08.064
% Formulas T2-1 -- T2-4

%---


MULTIPOLES = zeros(length(W.lambda_list), N_modes, 4);

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
 

    Eb = transpose([Exincl', Eyincl', Ezincl']);

    fac2 = (epsw_sca-1) * retconstantes("eps0"); 
    fac3 = 1i*omega;


    % Go through every mode
    for isol=1:N_modes-90

        % Get Normalized Mode
        Emx = QNM_resonator{1}.Ex(isol,:); 
        Emy = QNM_resonator{1}.Ey(isol,:); 
        Emz = QNM_resonator{1}.Ez(isol,:); 

        % Calculating extinced coefficient
        coef_w = omega/(QNM.omega(isol)-omega); 
        alpha_m = sum( ... 
            fac2*coef_w.*QNM_resonator{1}.mesh_vol.*( ... 
            Emx.*Exincl + Emy.*Eyincl + Emz.*Ezincl)); 
        
        % Calculating Mode Currents
        Jmx = fac2 .* fac3 .* alpha_m .* Emx;
        Jmy = fac2 .* fac3 .* alpha_m .* Emy;
        Jmz = fac2 .* fac3 .* alpha_m .* Emz;

        % Additional variables:
        r2 = X.*X + Y.*Y + Z.*Z;
        kr  = k0 .* sqrt(r2);
        

        j0 = sqrt(pi./(2*kr)).*besselj(0+1/2,kr);  
        j1 = sqrt(pi./(2*kr)).*besselj(1+1/2,kr);  
        j2 = sqrt(pi./(2*kr)).*besselj(2+1/2,kr);  
        j3 = sqrt(pi./(2*kr)).*besselj(3+1/2,kr);  
        
        J = transpose([Jmx', Jmy', Jmz']);
        R = transpose([X', Y', Z']);

        RJ = R(1,:) .* J(1,:) + R(2,:) .*J (2,:) +  R(3,:).*J(3,:);
        RxJ = transpose([(R(2,:).*J(3,:) - R(3,:).*J(2,:))',...
            (R(3,:).*J(1,:) - R(1,:).*J(3,:))',...
            (R(1,:).*J(2,:) - R(2,:).*J(1,:))']);

        const = k0^4/(6 * pi * retconstantes("eps0")^2 * abs(E0)^2);
        
        % ED
        p = zeros(3);
        for ic=1:3

            % Momentum
            p(ic) = -1/(1i*omega) * (sum((J(ic,:).*j0 ...
             + (k0^2)/2 * (3.*RJ.*R(ic,:) - r2 .* J(ic,:)) .*...
             j2./ (kr.^2) ).* QNM_resonator{1}.mesh_vol));

            % SCS
            MULTIPOLES(ilambda, 1) = MULTIPOLES(ilambda, 1) ...
                +  const .* (p(ic) .* conj(p(ic)));
        end
       

        % MD
        m = zeros(3);
        for ic=1:3

            % Momentum
            m(ic) = 3/2 * sum((RxJ(ic,:)).*...
                j1./kr.*QNM_resonator{1}.mesh_vol);
            
            % SCS
            MULTIPOLES(ilambda, isol, 2) = MULTIPOLES(ilambda, isol, 2) ...
                + const *(retconstantes("mu0") * retconstantes("eps0")).*abs( m(ic) .* conj(m(ic)));
        end


        % EQ
        EQ = zeros(3,3);

        
        for ic0=1:3
            for ic1=1:3
                
                % Kroneker  delta Definition
                if ic0==ic1
                    kronek = 1;
                else
                    kronek = 0;
                end

                % Momentum
                EQ(ic0, ic1) = -3/(1i*omega) * (...
                    sum((3*(R(ic1,:).*J(ic0,:) + R(ic0,:).*J(ic1,:))...
                    - 2 .* RJ .* kronek ) .* j1 ./ kr ...
                    .* QNM_resonator{1}.mesh_vol) + ...
                    2 .*(k0^2) .* sum((5.*R(ic0,:).*R(ic1,:).*RJ ...
                    - (R(ic0,:).*J(ic1,:) + R(ic1,:).*J(ic0,:)) .* r2...
                    - r2 .* (RJ) .* kronek) .*j3./(kr.^3) ...
                    .* QNM_resonator{1}.mesh_vol));

               % SCS
               MULTIPOLES(ilambda, isol, 3) = ...
                   MULTIPOLES(ilambda, isol, 3) + ...
                   const ./ 120 .* (k0.^2) .* abs( EQ(ic0, ic1) .* ...
                   conj(EQ(ic0, ic1)));
            end
        end
       

        % MQ
        MQ = zeros(3,3);
        
        for ic0=1:3
            for ic1=1:3
                
                % Momentum
                MQ(ic0, ic1) = 15 * (...
                    sum((R(ic0,:).*RxJ(ic1,:) ...
                    +R(ic1,:).*RxJ(ic0,:)) ...
                    .* j2./(kr.^2) ...
                    .* QNM_resonator{1}.mesh_vol));
               
               % SCS
               MULTIPOLES(ilambda, isol, 4) = ...
                   MULTIPOLES(ilambda, isol, 4) + ...
                   const *retconstantes("mu0") * retconstantes("eps0") ./ (30) .* ...
               (k0.^2) .* (MQ(ic0, ic1) .* conj(MQ(ic0, ic1)));
            end
        end


    end 
end 

% Plotting: For every mode Multipoles + Sum of M + Mode Extintion
multipole_names = ["ED" "MD" "EQ" "MQ"];
for isol=1:N_modes
    figure;
    plot(W.lambda_list.*1e9, scs_m(:,isol), "DisplayName", ...
        [strcat("Ext Mode #",num2str(isol))]);
    hold on;
    for imult=1:4
        plot(W.lambda_list.*1e9, transpose(MULTIPOLES(:,isol,imult))...
        ./(W.lambda_list.^2/(2*pi)),'DisplayName',multipole_names(imult));
    end
    plot(W.lambda_list.*1e9, transpose(sum(MULTIPOLES(:,isol,:),3))...
        ./(W.lambda_list.^2/(2*pi)),'DisplayName',"Sum of Mode");
    xlabel('wavelength');
    ylabel("\sigma_{EXT}, \sigma_{SCS}");
    legend();
end

% Plotting: Summmary Modes Multipoles + Total SCS + Extinction from prev. 
SUM_MD = sum(MULTIPOLES(:,:,:), 2);


figure;
for i=1:4
    plot(W.lambda_list.*1e9, transpose(SUM_MD(:,1,i))...
        ./(W.lambda_list.^2/(2*pi)),'DisplayName',multipole_names(i));
    hold on;
end
% 
TOTAL_MD = sum(SUM_MD(:,1,:),3);
 
plot(W.lambda_list.*1e9, transpose(TOTAL_MD)...
    ./(W.lambda_list.^2/(2*pi)),'DisplayName', 'SUM MD');

plot(W.lambda_list.*1e9, scs_sum, "DisplayName", "scs");
xlabel('wavelength');
ylabel("\sigma_{SCS}");
legend();

%% MD every mode test

c = retconstantes("c");
eps0 = retconstantes("eps0");
W.lambda_list = linspace(wl_min, wl_max, length(X));
omegas = 2*pi./W.lambda_list.*retconstantes("c"); 
k_list = 2 * pi ./ W.lambda_list;


MULTIPOLES = zeros(length(W.lambda_list), N_modes, 4);

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
 

    Eb = transpose([Exincl', Eyincl', Ezincl']);

    fac2 = (epsw_sca-1) * retconstantes("eps0"); 
    fac3 = 1i*omega;


    % Go through every mode
    for isol=1:N_modes

        % Get Normalized Mode
        Emx = QNM_resonator{1}.Ex(isol,:); 
        Emy = QNM_resonator{1}.Ey(isol,:); 
        Emz = QNM_resonator{1}.Ez(isol,:); 

        % Calculating extinced coefficient
        coef_w = omega/(QNM.omega(isol)-omega); 
        alpha_m = sum( ... 
            fac2*coef_w.*QNM_resonator{1}.mesh_vol.*( ... 
            Emx.*Exincl + Emy.*Eyincl + Emz.*Ezincl)); 
        
        % Calculating Mode Currents
        Jmx = fac2 .* fac3 .* alpha_m .* (Emx+Exincl);
        Jmy = fac2 .* fac3 .* alpha_m .* (Emy+Eyincl);
        Jmz = fac2 .* fac3 .* alpha_m .* (Emz+Ezincl);

        Jx = Jmx;
        Jy = Jmy;
        Jz = Jmz;
    
        x4d = X;
        y4d = Y;
        z4d = Z;
        k4d = k_list;
        k = k0;
        d3r = QNM_resonator{1}.mesh_vol;
    
        % calculate often used values
        % constant for scattering cross section
        const = k.^4/(6*pi*eps0^2*1);  % E0 = 1
        
        % scalar product
        rJ = x4d.*Jx + y4d.*Jy + z4d.*Jz;  % product r,J(r)
        rr = x4d.*x4d + y4d.*y4d + z4d.*z4d;  % product r,r
        r = sqrt(rr);  % norm(r)
        
        % cross product r x J = (ry*Jz-rz*Jy, rz*Jx-rx*Jz, rx*Jy-ry*Jx)
        rxJx = (y4d.*Jz - z4d.*Jy);
        rxJy = (z4d.*Jx - x4d.*Jz);
        rxJz = (x4d.*Jy - y4d.*Jx);
        
        % spherical bessel functions
        sbj0 = sqrt(pi./(2*k4d.*r)).*besselj(0+1/2,k4d.*r);  %j0
        sbj1 = sqrt(pi./(2*k4d.*r)).*besselj(1+1/2,k4d.*r);  %j1
        sbj2 = sqrt(pi./(2*k4d.*r)).*besselj(2+1/2,k4d.*r);  %j2
        sbj3 = sqrt(pi./(2*k4d.*r)).*besselj(3+1/2,k4d.*r);  %j3
        
        % calculate multipole moments and cross sections
        % calculate electric dipole p,
        dpx = (3*rJ.*x4d-rr.*Jx).*sbj2./(k4d.*r).^2;
        dpy = (3*rJ.*y4d-rr.*Jy).*sbj2./(k4d.*r).^2;
        dpz = (3*rJ.*z4d-rr.*Jz).*sbj2./(k4d.*r).^2;
        px = -1./(1i*omega).*(sum(Jx.*sbj0.*d3r)+ ...
        k.^2/2.*sum(dpx.*d3r));
        py = -1./(1i*omega).*(sum(Jy.*sbj0.*d3r)+ ...
        k.^2/2.*sum(dpy.*d3r));
        pz = -1./(1i*omega).*(sum(Jz.*sbj0.*d3r)+ ...
        k.^2/2.*sum(dpz.*d3r));
        norm2_p = px.*conj(px)+py.*conj(py)+pz.*conj(pz);
        Cp = const.*norm2_p;
        MULTIPOLES(ilambda, isol, 1) = Cp;
        
        % calculate magnetic dipole m,
        dmx = rxJx.*sbj1./(k4d.*r);
        dmy = rxJy.*sbj1./(k4d.*r);
        dmz = rxJz.*sbj1./(k4d.*r);
        mx = 3/2*sum(dmx.*d3r);
        my = 3/2*sum(dmy.*d3r);
        mz = 3/2*sum(dmz.*d3r);
        norm2_m = mx.*conj(mx)+my.*conj(my)+mz.*conj(mz);
        Cm = const.*norm2_m/c^2;
        MULTIPOLES(ilambda, isol, 2) = Cm;
        
        % calculate electric quadrupole Qe
        dQe1xx = (3*2*x4d.*Jx - 2*rJ).*sbj1./(k4d.*r);
        dQe1xy = (3*(y4d.*Jx+x4d.*Jy)).*sbj1./(k4d.*r);
        dQe1xz = (3*(z4d.*Jx+x4d.*Jz)).*sbj1./(k4d.*r);
        dQe1yy = (3*2*y4d.*Jy - 2.*rJ).*sbj1./(k4d.*r);
        dQe1yx = (3*(x4d.*Jy+y4d.*Jx)).*sbj1./(k4d.*r);
        dQe1yz = (3*(z4d.*Jy+y4d.*Jz)).*sbj1./(k4d.*r);
        dQe1zz = (3*2*z4d.*Jz - 2*rJ).*sbj1./(k4d.*r);
        dQe1zx = (3*(x4d.*Jz+z4d.*Jx)).*sbj1./(k4d.*r);
        dQe1zy = (3*(y4d.*Jz+z4d.*Jy)).*sbj1./(k4d.*r);
        dQe2xx = (5*x4d.*x4d.*rJ - rr*2.*x4d.*Jx - rr.*rJ).*sbj3./(k4d.*r).^3;
        dQe2xy = (5*x4d.*y4d.*rJ - rr.*(x4d.*Jy+y4d.*Jx)).*sbj3./(k4d.*r).^3;
        dQe2xz = (5*x4d.*z4d.*rJ - rr.*(x4d.*Jz+z4d.*Jx)).*sbj3./(k4d.*r).^3;
        dQe2yy = (5*y4d.*y4d.*rJ - rr*2.*y4d.*Jy - rr.*rJ).*sbj3./(k4d.*r).^3;
        dQe2yx = (5*y4d.*x4d.*rJ - rr.*(y4d.*Jx+x4d.*Jy)).*sbj3./(k4d.*r).^3;
        dQe2yz = (5*y4d.*z4d.*rJ - rr.*(y4d.*Jz+z4d.*Jy)).*sbj3./(k4d.*r).^3;
        dQe2zz = (5*z4d.*z4d.*rJ - rr*2.*z4d.*Jz - rr.*rJ).*sbj3./(k4d.*r).^3;
        dQe2zx = (5*z4d.*x4d.*rJ - rr.*(z4d.*Jx+x4d.*Jz)).*sbj3./(k4d.*r).^3;
        dQe2zy = (5*z4d.*y4d.*rJ - rr.*(z4d.*Jy+y4d.*Jz)).*sbj3./(k4d.*r).^3;
        Qexx = -3./(1i*omega).*(sum(dQe1xx.*d3r)+ ...
        2*k.^2.*sum(dQe2xx.*d3r));
        Qexy = -3./(1i*omega).*(sum(dQe1xy.*d3r)+ ...
        2*k.^2.*sum(dQe2xy.*d3r));
        Qexz = -3./(1i*omega).*(sum(dQe1xz.*d3r)+ ...
        2*k.^2.*sum(dQe2xz.*d3r));
        Qeyy = -3./(1i*omega).*(sum(dQe1yy.*d3r)+ ...
        2*k.^2.*sum(dQe2yy.*d3r));
        Qeyx = -3./(1i*omega).*(sum(dQe1yx.*d3r)+ ...
        2*k.^2.*sum(dQe2yx.*d3r));
        Qeyz = -3./(1i*omega).*(sum(dQe1yz.*d3r)+ ...
        2*k.^2.*sum(dQe2yz.*d3r));
        Qezz = -3./(1i*omega).*(sum(dQe1zz.*d3r)+ ...
        2*k.^2.*sum(dQe2zz.*d3r));
        Qezx = -3./(1i*omega).*(sum(dQe1zx.*d3r)+ ...
        2*k.^2.*sum(dQe2zx.*d3r));
        Qezy = -3./(1i*omega).*(sum(dQe1zy.*d3r)+ ...
        2*k.^2.*sum(dQe2zy.*d3r));
        norm2_Qe = Qexx.*conj(Qexx)+Qexy.*conj(Qexy)+Qexz.*conj(Qexz)+ ...
           Qeyy.*conj(Qeyy)+Qeyx.*conj(Qeyx)+Qeyz.*conj(Qeyz)+ ...
           Qezz.*conj(Qezz)+Qezx.*conj(Qezx)+Qezy.*conj(Qezy);
        CQe = const/120.*k.^2.*norm2_Qe;
        MULTIPOLES(ilambda, isol, 3) = CQe;
        
        % calculate magnetic quadrupole Qm
        dQmxx = (2*x4d.*rxJx).*sbj2./(k4d.*r).^2;
        dQmxy = (x4d.*rxJy+y4d.*rxJx).*sbj2./(k4d.*r).^2;
        dQmxz = (x4d.*rxJz+x4d.*rxJz).*sbj2./(k4d.*r).^2;
        dQmyy = (2*y4d.*rxJy).*sbj2./(k4d.*r).^2;
        dQmyx = (y4d.*rxJx+x4d.*rxJy).*sbj2./(k4d.*r).^2;
        dQmyz = (y4d.*rxJz+z4d.*rxJy).*sbj2./(k4d.*r).^2;
        dQmzz = (2*z4d.*rxJz).*sbj2./(k4d.*r).^2;
        dQmzx = (z4d.*rxJx+x4d.*rxJz).*sbj2./(k4d.*r).^2;
        dQmzy = (z4d.*rxJy+y4d.*rxJz).*sbj2./(k4d.*r).^2;
        Qmxx = 15*sum(dQmxx.*d3r);
        Qmxy = 15*sum(dQmxy.*d3r);
        Qmxz = 15*sum(dQmxz.*d3r);
        Qmyy = 15*sum(dQmyy.*d3r);
        Qmyx = 15*sum(dQmyx.*d3r);
        Qmyz = 15*sum(dQmyz.*d3r);
        Qmzz = 15*sum(dQmzz.*d3r);
        Qmzx = 15*sum(dQmzx.*d3r);
        Qmzy = 15*sum(dQmzy.*d3r);
        norm2_Qm = Qmxx.*conj(Qmxx)+Qmxy.*conj(Qmxy)+Qmxz.*conj(Qmxz)+ ...
           Qmyy.*conj(Qmyy)+Qmyx.*conj(Qmyx)+Qmyz.*conj(Qmyz)+ ...
           Qmzz.*conj(Qmzz)+Qmzx.*conj(Qmzx)+Qmzy.*conj(Qmzy);
        CQm = const./120.*(k/c).^2.*norm2_Qm;
        MULTIPOLES(ilambda, isol, 4) = CQm;
   

    end 
end 

% Plotting: For every mode Multipoles + Sum of M + Mode Extintion
multipole_names = ["ED" "MD" "EQ" "MQ"];
for isol=1:N_modes
    figure;
    plot(linspace(425e-9,455e-9,201).*1e9, scs_m(:,isol), "DisplayName", ...
        [strcat("Ext Mode #",num2str(isol))]);
    hold on;
    for imult=1:4
        plot(W.lambda_list.*1e9, transpose(MULTIPOLES(:,isol,imult))...
        ./normalization,'DisplayName',multipole_names(imult));
    end
    plot(W.lambda_list.*1e9, transpose(sum(MULTIPOLES(:,isol,:),3))...
        ./normalization,'DisplayName',"Sum of Mode");
    xlabel('wavelength');
    ylabel("\sigma_{EXT}, \sigma_{SCS}");
    legend();
end

% Plotting: Summmary Modes Multipoles + Total SCS + Extinction from prev. 
SUM_MD = sum(MULTIPOLES(:,:,:), 2);


figure;
for i=1:4
    plot(W.lambda_list.*1e9, transpose(SUM_MD(:,1,i))...
        ./normalization,'DisplayName',multipole_names(i));
    hold on;
end
% 
TOTAL_MD = sum(SUM_MD(:,1,:),3);
 
plot(W.lambda_list.*1e9, transpose(TOTAL_MD)...
    ./normalization,'DisplayName', 'SUM MD');

plot(linspace(425e-9,455e-9,201).*1e9, scs_sum, "DisplayName", "scs");
xlabel('wavelength');
ylabel("\sigma_{SCS}");
legend();

