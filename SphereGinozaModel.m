clc
clf
clear
Ang = char(197);

% 2022
% Code constructed by Ashley P. Williams, SMaCLab, Monash University, Australia

%--------------------------------------------------------------------------
%              Spherical Form Factor + Ginoza MSA Structure factor. 
%--------------------------------------------------------------------------
%
% This model is the spherical model with a Hayter-Penfold-Hansen-Ginoza
% structure factor to model charged spherical particles.
%
% If you use this model, please cite the paper below regarding the Ginoza
% structure factor:
% "A simple and accurate method for calculation of the structure factor of 
%  interacting charged spheres", J. Colloid Interface Sci. 2014.
%
% Load your file below in the appropriate section, where your data file is
% of the form of 3 column vectors: q, I, err. Also ensure that you cut any
% points with NaN values, to minimise any issues.
%
% To save the fit, open the vector "I_spherical_ginoza" from the workspace and
% copy the contents into a text file.
%
% To save the parameters, open the vector "parameters" and save as desired.
% The parameters are saved as a 3 column vector: parameters, values
% unnormalised and values normalised.
%
% The command window below will indicate whether you need to renormalise.
% To find an optimal renormalisation value, your MSA contact value should 
% be positive and within 3 decimals of zero.
%               e.g. 0 <= MSA contact value < 0.001.
% The plot of s (renormalisation factor) and MSA contact value will help
% you find this value and indicate where you are on the curve. You should
% choose a value of s closest to 1, where the curve crosses the x axis.
% 
%   ________  ___      ___       __       ______   ___            __       _______   
%  /"       )|"  \    /"  |     /""\     /" _  "\ |"  |          /""\     |   _  "\  
% (:   \___/  \   \  //   |    /    \   (: ( \___)||  |         /    \    (. |_)  :) 
%  \___  \    /\\  \/.    |   /' /\  \   \/ \     |:  |        /' /\  \   |:     \/  
%   __/  \\  |: \.        |  //  __'  \  //  \ _   \  |___    //  __'  \  (|  _  \\  
%  /" \   :) |.  \    /:  | /   /  \\  \(:   _) \ ( \_|:  \  /   /  \\  \ |: |_)  :) 
% (_______/  |___|\__/|___|(___/    \___)\_______) \_______)(___/    \___)(_______/ 
% 
%
% SMaCLab website can be found here:
% https://sites.google.com/view/smaclab
%
% SMaCLab github can be found here:
% https://github.com/SMaC-3
%            
%--------------------------------------------------------------------------
%                               Load File
%--------------------------------------------------------------------------

path = ""; % "../data/" for example, to locate your data.
filename = path + "data_filename";

load(filename)

data = data_filename;

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 5e-1;
xlower = 4e-3;
yupper = 1e1;
ylower = 3e-2;

ShowStructureFactorPlot = 0; % 1 = show SF plot, 0 = Don't show SF plot.

%--------------------------------------------------------------------------
%                          Sphere parameters
%--------------------------------------------------------------------------

volumeFraction = 0.034;
radius = 23.8;
SLD = -0.392;
SLD_solv = 6.393;
background = 0.033;
PD = 0.16;              %Schulz polydispersity              

%--------------------------------------------------------------------------
%                    Ginoza structure factor parameters
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%               Which parameter set would you like to use?
%         Dimensionless parameters = 0 or Physical parameters = 1
%--------------------------------------------------------------------------

whichParams = 1;

%--------------------------------------------------------------------------
%                 User Input for Dimensionless Parameters
if whichParams == 0   
%--------------------------------------------------------------------------

eta = volumeFraction; % Volume Fraction
k = 1.730086388;      % Debye Length * effective particle diameter.
gamma_0 = 117.56;     % Dimensionless energy parameter
Rad = 35;             % Particle radius
s = 1;                % scale for renormalisation

%--------------------------------------------------------------------------
%                   User input for Physical Parameters
elseif whichParams == 1
%--------------------------------------------------------------------------

chargePerParticle = 24;    % Charge per particle
dielectricConstant = 78;   % Dielectric constant
Rad = 23.8;                  % Particle radius
ionicStrength = 0.029;     % Ionic strength of solution
temperature = 298.15;         % Temperature (K)
counterIonCharge = 1;      % Counterion charge
s = 0.6725904;               % Scale for renormalisation

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%     Calculation of physical parameters into dimensionless parameters
%--------------------------------------------------------------------------
avogadrosNum = 6.0221415e23;
electricCharge = 1.60217657e-19;
permFreeSpace = 8.85418782e-12;
boltzmannConstant = 1.3806488e-23;

diameter = 2*Rad*1e-10;
volume_effective = pi*(diameter)^3/6;
n_0 = ((volumeFraction/volume_effective)/1000)/avogadrosNum; 
n_c = n_0*abs(chargePerParticle);
IonStr_n_c = 0.5*n_c*counterIonCharge^2;
IonStr_n_0 = 0.5*n_0*chargePerParticle^2;
IonStrIons = IonStr_n_c+ionicStrength;
IonStrAll = IonStrIons+IonStr_n_0;

debyeLength = sqrt(dielectricConstant*permFreeSpace*boltzmannConstant*temperature/(2000*avogadrosNum*electricCharge^2*IonStrIons));
kappa = 1/debyeLength;
potential = chargePerParticle*electricCharge/(pi*permFreeSpace*dielectricConstant*diameter*(2+(1/debyeLength)*diameter));

eta = volumeFraction;
k = diameter/debyeLength;
gamma_0 = pi*dielectricConstant*permFreeSpace*diameter*potential^2*exp(diameter/debyeLength)/(boltzmannConstant*temperature);
end

%--------------------------------------------------------------------------
%                          Spherical calculation
%--------------------------------------------------------------------------

volume = @ (r) (4*pi/3)*r^3;
contrast = (SLD-SLD_solv);

I_spherical = zeros(size(data,1),2);
Form = @(r,q) (1e-4)*(3*contrast*volume(r)*((sin(q*r)-q*r*cos(q*r))/(q*r)^3))^2;

%--------------------------------------------------------------------------
%                       Schulz-Zimm Polydispersity
%--------------------------------------------------------------------------

Z = (1-PD^2)/PD^2;
sz = @(y) (((Z+1)/radius)^(Z+1))*((y)^(Z))*exp(-((Z+1)/radius)*y)...
    /(gamma(Z+1));
N_steps = 80;
dist_r = [];
rvals = ((radius-3*PD*radius):(6*PD*radius)/N_steps:(radius+3*PD*radius)); %Number of points for PD calculation.
for j = 1:size(rvals,2)
    dist_r(j) = sz(rvals(j)); 
end
dist_r = dist_r/max(dist_r); %Normalise the distribution weights.

if PD ~= 0
    form = [];
    vol = [];
    for h = 1:(size(data(:,1),1))
        for i = 1:size(rvals,2)
           form(i,h)= (dist_r(i)*Form(rvals(i),data(h,1)));
        end
    end
    
    for i = 1:size(rvals,2)
       vol(i) = dist_r(i)*volume(rvals(i));
    end
  
    form = sum(form)/sum(vol);
 
    else
        form = [];
        for h = 1:(size(data(:,1),1))
            form(1,h) = Form(radius,data(h,1));
        end
end

I_spherical(:,2) = form';
I_spherical(:,1) = data(:,1); %[q,I]
%--------------------------------------------------------------------------
%                    Structure factor calculation
%--------------------------------------------------------------------------
eta_input = eta;
k_input = k;
gamma_0_input = gamma_0;

eta = eta/s^3;
k = k/s;
gamma_0 = gamma_0*s;

diameter = 2*Rad;
x = 1;
K_w = -gamma_0*exp(-k);
fy = 3*eta/(1-eta);
PSI1 = (1/k^3)*(x-(k/2)-(x+k/2)*exp(-k/x));
phi0 = (1/k)*(1-exp(-k/x));
PHI1 = x*phi0-4*fy*PSI1;
PHI0 = x+(fy*phi0)-(4*fy*PSI1)*(1+(k/2)+fy);

%--------------------------Quartic-Terms-----------------------------------

a = PHI0/(k*PHI1);
b = -6*eta*K_w/(k^4*PHI1^2);

p_qt = -0.25*(a^2-a+0.75);
x_plus = 0.5*(-(a/2)-0.75+sqrt(3/8-p_qt));

a4 = 1;
a3 = 2*a+1;
a2 = a^2+2*a;
a1 = a^2;
a0 = -b;

d3 = 4*a4;
d2 = 3*a3;
d1 = 2*a2;
d0 = a1;

x_sol(1) = x_plus+0.5;

% Newton-Raphson for finding x_sol,   x_(n+1) = x_n + f(x_n)/f'(x_n)
for i = 1:20
    f = a4*x_sol(i)^4+a3*x_sol(i)^3+a2*x_sol(i)^2+a1*x_sol(i)+a0;
    fdiv = d3*x_sol(i)^3+d2*x_sol(i)^2+d1*x_sol(i)+d0;
    x_sol(i+1) = x_sol(i)-(f/fdiv);
end

x_soln = x_sol(end);
Gamma = x_soln*k;

alpha_1 = 12*(eta/(1-eta))*(1/(k^2))*(1+(k/2));
alpha_0 = alpha_1*(1+(k/2)+fy)-fy;
v = 2*K_w*(alpha_0+(alpha_1-1)*Gamma)/(PHI0+PHI1*Gamma);

%--------------------------------------------------------------------------
%                        HP parameter calculation
%
% alpha will be reused as a variable but the above alpha aren't needed
% anymore, just need to find v. Added an underscore to above vars.
%--------------------------------------------------------------------------

alpha0 = 1+ (eta/2);
alpha1 = 1+2*eta;
alpha2 = 6*eta;
alpha4 = 24*eta;
alpha5 = 1-2*(eta^2)-8*eta;

beta0 = 1;
beta1 = 2;
beta2 = 2;
beta4 = 0;
beta5 = 6*eta;

gamma0 = K_w+(v/k)*(1-exp(-k))+v^2*(exp(-k)/(2*k^2*K_w))*(cosh(k)-1);
gamma1 = -K_w*k+v*exp(-k)+v^2*(exp(-k)/(2*k*K_w))*sinh(k);
gamma2 = K_w*k^2-v*k*exp(-k)+v^2*(exp(-k)/(2*K_w))*cosh(k);
gamma4 = -2*v*k^3+v^2*k^2*exp(-k)/K_w;
gamma5 = 12*(eta/k)*(1-(2/k^2)+2*(k+1)*exp(-k)/k^2)*v + (6*eta*exp(-k)/(K_w*k^2))*((2/k)*sinh(k)-(2/k^2)*cosh(k)+(2-k^2)/k^2)*v^2-24*eta*K_w*(k+1)/k^2+1;

m5 = beta5/alpha5;
n5 = gamma5/alpha5;
m0 = alpha0*m5+beta0;
n0 = alpha0*n5+gamma0;
m1 = alpha1*m5+beta1;
n1 = alpha1*n5+gamma1;
m2 = alpha2*m5+beta2;
n2 = alpha2*n5+gamma2;
B_linear = 2*(m1*n1-m0*n2-m2*n0)-m5;
C_linear = n1^2-2*n0*n2-n5-(gamma4/alpha4);

b_term = -C_linear/B_linear;
a_term = m5*b_term+n5;

A = -a_term;
B = -b_term;
C = -v/k;
F = v/k - (v^2*exp(-k))/(2*K_w*k^2);

MSA_contact_val = alpha0*a_term+beta0*b_term+gamma0; %This should be > 0. if not, renormalise parameters.

q = data(:,1);
k_vec = q*diameter/s;
I_spherical_ginoza = zeros(size(data,1),2);
I_spherical_ginoza(:,1) = data(:,1); 

for i = 1:size(k_vec,1)
    term1 = A*(sin(k_vec(i))-k_vec(i)*cos(k_vec(i)))/k_vec(i)^3;
    term2 = B*((2/k_vec(i)^2 - 1)*k_vec(i)*cos(k_vec(i))+2*sin(k_vec(i))-2/k_vec(i))/k_vec(i)^3;
    term3 = eta*A*(24/k_vec(i)^3+4*(1-6/k_vec(i)^2)*sin(k_vec(i))-(1-12/k_vec(i)^2+24/k_vec(i)^4)*k_vec(i)*cos(k_vec(i)))/(2*k_vec(i)^3);
    term4 = C*(k*cosh(k)*sin(k_vec(i))-k_vec(i)*sinh(k)*cos(k_vec(i)))/(k_vec(i)*(k_vec(i)^2+k^2));
    term5 = F*(k*sinh(k)*sin(k_vec(i))-k_vec(i)*(cosh(k)*cos(k_vec(i))-1))/(k_vec(i)*(k_vec(i)^2+k^2));
    term6 = F*(cos(k_vec(i))-1)/k_vec(i)^2;
    term7 = -gamma_0*exp(-k)*(k*sin(k_vec(i))+k_vec(i)*cos(k_vec(i)))/(k_vec(i)*(k_vec(i)^2+k^2));

    sum(i) = term1+term2+term3+term4+term5+term6+term7;
    S_K(i) = 1/(1-24*eta*sum(i));
    I_spherical_ginoza(i,2) = I_spherical(i,2).*S_K(i);
end
if PD ~= 0
    I_spherical_ginoza(:,2) = (volumeFraction).*I_spherical_ginoza(:,2)+background;
else
    I_spherical_ginoza(:,2) = (volumeFraction/volume(radius)).*I_spherical_ginoza(:,2)+background;
end
    Ginoza_structure_factor(:,1) = data(:,1);
    Ginoza_structure_factor(:,2) = S_K;

figure(2)
hold on
box on
xlabel("{\it q} ("+Ang+"^{-1})")
ylabel("Intensity (cm^{-1})")
set(gca,'Xscale', 'log','Yscale', 'log','fontweight','bold','Linewidth',2,'fontSize',14)
xlim([xlower, xupper])
ylim([ylower yupper])
errorbar(data(:,1),data(:,2),data(:,3),data(:,3),'o','MarkerFaceColor','auto','MarkerSize', 6,'Color','Black')
plot(I_spherical_ginoza(:,1),I_spherical_ginoza(:,2),'Color',[0 0 1],'Linewidth',3)
legend("Data","Sphere-Ginoza")

if ShowStructureFactorPlot == 1
    figure(3)
    hold on
    box on
    xlabel("{\it q} ("+Ang+"^{-1})")
    ylabel("Intensity (cm^{-1})")
    set(gca,'Xscale', 'log','Yscale', 'linear','fontweight','bold','Linewidth',2,'fontSize',14)
    xlim([xlower, xupper])
    plot(Ginoza_structure_factor(:,1),Ginoza_structure_factor(:,2),'Color',[0 0 0],'Linewidth',3)
    legend("GinozaSF","location","southeast")
end

%--------------------------------------------------------------------------
%                            Renormalisation
%--------------------------------------------------------------------------

if MSA_contact_val ~= 0
    s_vals = renorm(eta_input,k_input,gamma_0_input,radius);
    figure(1)
    hold on
    box on 
    set(gca,'Xscale', 'linear','fontweight','bold','Linewidth',2.6,'fontsize',14)
    set(gca,'Yscale', 'linear');
    xlim([0.05 1])
    ylim([-5 5])
    xlabel("s")
    ylabel("MSA contact value")
    f = @(x) 0;
    plot(s_vals(:,1),s_vals(:,2),'Color',[0 0 0],'linewidth',2)
    fplot(f,"--",'color','black')
    plot(s,MSA_contact_val,'color','red','marker','x','linewidth',2)
    
    if MSA_contact_val < 0
        error("Renormalise parameters, MSA contact value = "+MSA_contact_val+". MSA contact value should be positive and as close to 0 as possible (three decimal places).")
    end
end

if MSA_contact_val < 0.001
    "MSA contact value is "+MSA_contact_val+". This is within 3 decimal places of zero :)"
    "Your rescaled parameters are: eta = "+eta+", kappa = "+k+", gamma_0 = "+gamma_0+"."
else
    "MSA contact value is "+MSA_contact_val+". Please consider rescaling so that the MSA contact value is within 3 decimal places of zero."
end
    
parameters = ["RenormFactor" 1 s;"MSA Contact Factor" s_vals(1,2) MSA_contact_val;"volumeFraction" eta_input eta; "Kappa" k_input k; "gamma_0" gamma_0_input gamma_0; "Radius used in S(q)" Rad Rad; "Radius" radius radius; ; "SLD" SLD SLD; "SLDsolvent" SLD_solv SLD_solv; "Background" background background];
function s_vals = renorm(eta, k, gamma_0,effectiveRadius)

    s_vals = zeros(40,2);
    s_vals(:,1) = [1:-0.025:0.025];
    eta_orig = eta;
    k_orig = k;
    gamma_0_orig = gamma_0;
    for j = 1:40
        eta = eta_orig/s_vals(j,1)^3;
        k = k_orig/s_vals(j,1);
        gamma_0 = gamma_0_orig*s_vals(j,1);
    diameter = 2*effectiveRadius;
    x = 1;
    K_w = -gamma_0*exp(-k);
    fy = 3*eta/(1-eta);
    PSI1 = (1/k^3)*(x-(k/2)-(x+k/2)*exp(-k/x));
    phi0 = (1/k)*(1-exp(-k/x));
    PHI1 = x*phi0-4*fy*PSI1;
    PHI0 = x+(fy*phi0)-(4*fy*PSI1)*(1+(k/2)+fy);

    %--------------------------Quartic-Terms-----------------------------------

    a = PHI0/(k*PHI1);
    b = -6*eta*K_w/(k^4*PHI1^2);

    p_qt = -0.25*(a^2-a+0.75);
    x_plus = 0.5*(-(a/2)-0.75+sqrt(3/8-p_qt));

    a4 = 1;
    a3 = 2*a+1;
    a2 = a^2+2*a;
    a1 = a^2;
    a0 = -b;

    d3 = 4*a4;
    d2 = 3*a3;
    d1 = 2*a2;
    d0 = a1;

    x_sol(1) = x_plus+0.5;

    %Newton-Raphson for finding x_sol,   x_(n+1) = x_n + f(x_n)/f'(x_n)
    for i = 1:20
        f = a4*x_sol(i)^4+a3*x_sol(i)^3+a2*x_sol(i)^2+a1*x_sol(i)+a0;
        fdiv = d3*x_sol(i)^3+d2*x_sol(i)^2+d1*x_sol(i)+d0;
        x_sol(i+1) = x_sol(i)-(f/fdiv);
    end

    x_soln = x_sol(end);
    Gamma = x_soln*k;

    alpha_1 = 12*(eta/(1-eta))*(1/(k^2))*(1+(k/2));
    alpha_0 = alpha_1*(1+(k/2)+fy)-fy;
    v = 2*K_w*(alpha_0+(alpha_1-1)*Gamma)/(PHI0+PHI1*Gamma);

    %--------------------------HP---Parameters---------------------------------
    %alpha will be reused as a variable but the above alpha aren't needed
    %anymore, just need to find v. Added an underscore to above vars anyway.

    alpha0 = 1+ (eta/2);
    alpha1 = 1+2*eta;
    alpha2 = 6*eta;
    alpha4 = 24*eta;
    alpha5 = 1-2*(eta^2)-8*eta;

    beta0 = 1;
    beta1 = 2;
    beta2 = 2;
    beta4 = 0;
    beta5 = 6*eta;

    gamma0 = K_w+(v/k)*(1-exp(-k))+v^2*(exp(-k)/(2*k^2*K_w))*(cosh(k)-1);
    gamma1 = -K_w*k+v*exp(-k)+v^2*(exp(-k)/(2*k*K_w))*sinh(k);
    gamma2 = K_w*k^2-v*k*exp(-k)+v^2*(exp(-k)/(2*K_w))*cosh(k);
    gamma4 = -2*v*k^3+v^2*k^2*exp(-k)/K_w;
    gamma5 = 12*(eta/k)*(1-(2/k^2)+2*(k+1)*exp(-k)/k^2)*v + (6*eta*exp(-k)/(K_w*k^2))*((2/k)*sinh(k)-(2/k^2)*cosh(k)+(2-k^2)/k^2)*v^2-24*eta*K_w*(k+1)/k^2+1;

    m5 = beta5/alpha5;
    n5 = gamma5/alpha5;
    m0 = alpha0*m5+beta0;
    n0 = alpha0*n5+gamma0;
    m1 = alpha1*m5+beta1;
    n1 = alpha1*n5+gamma1;
    m2 = alpha2*m5+beta2;
    n2 = alpha2*n5+gamma2;
    B_linear = 2*(m1*n1-m0*n2-m2*n0)-m5;
    C_linear = n1^2-2*n0*n2-n5-(gamma4/alpha4);

    b_term = -C_linear/B_linear;
    a_term = m5*b_term+n5;

    A = -a_term;
    B = -b_term;
    C = -v/k;
    F = v/k - (v^2*exp(-k))/(2*K_w*k^2);

    s_vals(j,2) = alpha0*a_term+beta0*b_term+gamma0;
    end
end



