clc
clf
clear
Ang = char(197);

% 2021
% Code constructed by Ashley P. Williams, SMaCLab, Monash University, Australia

%--------------------------------------------------------------------------
%              Ellipsoid Form Factor + Ginoza MSA Structure factor. 
%--------------------------------------------------------------------------
%
% This model is the ellipsoid model with a Hayter-Penfold-Hansen-Ginoza
% structure factor to model charged ellipsoidal particles.
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
% To save the fit, open the vector "I_ellipsoid_ginoza" from the workspace and
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

filename = "YourFileNameHere";
load(filename);
data = YourFileNameHere;

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 4e-1;
xlower = 1.5e-3;
yupper = 1e1;
ylower = 1e-2;

ShowStructureFactorPlot = 0; % 1 = show SF plot, 0 = Don't show SF plot.

%--------------------------------------------------------------------------
%                          Ellipsoid parameters
%--------------------------------------------------------------------------

volumeFraction = 0.0225;
radius_polar = 16.5;
radius_equatorial = 27;
SLD = 0.381;
SLD_solv = 6.3;
background = 0.02;

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
radius = 35;          % particle radius
s = 1;                % scale for renormalisation

%--------------------------------------------------------------------------
%                   User input for Physical Parameters
elseif whichParams == 1
%--------------------------------------------------------------------------

chargePerParticle = 45;    % Charge per particle
dielectricConstant = 78;   % Dielectric constant
radius = 24.1;               % particle radius
ionicStrength = 0.007;     % Ionic strength of solution
temperature = 298;         % Temperature (K)
counterIonCharge = 1;      % Counterion charge
s = 0.48687;                % scale for renormalisation

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%     Calculation of physical parameters into dimensionless parameters
%--------------------------------------------------------------------------
avogadrosNum = 6.0221415e23;
electricCharge = 1.60217657e-19;
permFreeSpace = 8.85418782e-12;
boltzmannConstant = 1.3806488e-23;



diameter = 2*radius*1e-10;
volume = pi*(diameter)^3/6;
n_0 = ((volumeFraction/volume)/1000)/avogadrosNum; 
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
%                        Form factor calculation
%--------------------------------------------------------------------------

vol = (4*pi/3)*radius^3;
v_square_minus_one = (radius_polar/radius_equatorial)^2-1;
zm = 0.5;
zb = 0.5;

gauss_N = 76; %According to /lib/Gauss76.c
gauss_W = [.00126779163408536,.00294910295364247,.00462793522803742,.00629918049732845,.00795984747723973,.00960710541471375,.0112381685696677,.0128502838475101,.0144407317482767,.0160068299122486,.0175459372914742,.0190554584671906,.020532847967908,.0219756145344162,.0233813253070112,.0247476099206597,.026072164497986,.0273527555318275,.028587223650054,.029773487255905,.0309095460374916,.0319934843404216,.0330234743977917,.0339977794120564,.0349147564835508,.0357728593807139,.0365706411473296,.0373067565423816,.0379799643084053,.0385891292645067,.0391332242205184,.0396113317090621,.0400226455325968,.040366472122844,.0406422317102947,.0408494593018285,.040987805464794,.0410570369162294,.0410570369162294,.040987805464794,.0408494593018285,.0406422317102947,.040366472122844,.0400226455325968,.0396113317090621,.0391332242205184,.0385891292645067,.0379799643084053,.0373067565423816,.0365706411473296,.0357728593807139,.0349147564835508,.0339977794120564,.0330234743977917,.0319934843404216,.0309095460374916,.029773487255905,.028587223650054,.0273527555318275,.026072164497986,.0247476099206597,.0233813253070112,.0219756145344162,.020532847967908,.0190554584671906,.0175459372914742,.0160068299122486,.0144407317482767,.0128502838475101,.0112381685696677,.00960710541471375,.00795984747723973,.00629918049732845,.00462793522803742,.00294910295364247,.00126779163408536];
gauss_Z = [-.999505948362153,-.997397786355355,-.993608772723527,-.988144453359837,-.981013938975656,-.972229228520377,-.961805126758768,-.949759207710896,-.936111781934811,-.92088586125215,-.904107119545567,-.885803849292083,-.866006913771982,-.844749694983342,-.822068037328975,-.7980001871612,-.77258672828181,-.74587051350361,-.717896592387704,-.688712135277641,-.658366353758143,-.626910417672267,-.594397368836793,-.560882031601237,-.526420920401243,-.491072144462194,-.454895309813726,-.417951418780327,-.380302767117504,-.342012838966962,-.303146199807908,-.263768387584994,-.223945802196474,-.183745593528914,-.143235548227268,-.102483975391227,-.0615595913906112,-.0205314039939986,.0205314039939986,.0615595913906112,.102483975391227,.143235548227268,.183745593528914,.223945802196474,.263768387584994,.303146199807908,.342012838966962,.380302767117504,.417951418780327,.454895309813726,.491072144462194,.526420920401243,.560882031601237,.594397368836793,.626910417672267,.658366353758143,.688712135277641,.717896592387704,.74587051350361,.77258672828181,.7980001871612,.822068037328975,.844749694983342,.866006913771982,.885803849292083,.904107119545567,.92088586125215,.936111781934811,.949759207710896,.961805126758768,.972229228520377,.981013938975656,.988144453359837,.993608772723527,.997397786355355,.999505948362153];

I_ellipsoid = zeros(size(data,1),1);
for j = 1:size(data,1)
    for i = 1:gauss_N
       u = gauss_Z(i)*zm +zb;
       r = radius_equatorial*sqrt(1+u^2*v_square_minus_one);
       f = 3*(sin(data(j,1)*r)/(data(j,1)*r)-cos(data(j,1)*r))/(data(j,1)*r)^2;
       I_ellipsoid(j,1) = I_ellipsoid(j,1)+gauss_W(i)*f^2; 
    end
    I_ellipsoid(j,1) = zm.*I_ellipsoid(j,1);
end

cont_vol = (SLD-SLD_solv)*vol;
I_ellipsoid(:,1) = (cont_vol^2*1e-4).*I_ellipsoid(:,1); %F^2


%--------------------------------------------------------------------------
%                    Structure factor calculation
%--------------------------------------------------------------------------
eta_input = eta;
k_input = k;
gamma_0_input = gamma_0;

eta = eta/s^3;
k = k/s;
gamma_0 = gamma_0*s;

diameter = 2*radius;
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
I_ellipsoid_ginoza = zeros(size(data,1),2);
I_ellipsoid_ginoza(:,1) = data(:,1); 

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
    I_ellipsoid_ginoza(i,2) = I_ellipsoid(i,1).*S_K(i);
end
    I_ellipsoid_ginoza(:,2) = (volumeFraction/vol).*I_ellipsoid_ginoza(:,2)+background;
    Ginoza_structure_factor(:,1) = data(:,1);
    Ginoza_structure_factor(:,2) = S_K;

fig = figure(1);
set(fig, 'units', 'centimeters', 'position', [10 10 35 15])
subplot(1,2,1)
hold on
box on
set(gca,'Xscale', 'log','Yscale', 'log','fontweight','bold','Linewidth',2,'fontSize',14)
xlabel("{\it q} ("+Ang+"^{-1})")
ylabel("Intensity (cm^{-1})")
xlim([xlower, xupper])
ylim([ylower yupper])
errorbar(data(:,1),data(:,2),data(:,3),data(:,3),'o','MarkerFaceColor','auto','MarkerSize', 6,'Color','Black')
plot(I_ellipsoid_ginoza(:,1),I_ellipsoid_ginoza(:,2),'Color',[1 0.5 0],'Linewidth',3)
legend("data","Ellipsoid-Ginoza","Location","southwest")

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
    subplot(1,2,2)
    hold on
    box on
    set(gca,'Xscale', 'linear','fontweight','bold','Linewidth',2.6,'fontsize',14)
    set(gca,'Yscale', 'linear');
    xlim([0.05 1])
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
    
parameters = ["RenormFactor" 1 s;"MSA Contact Factor" s_vals(1,2) MSA_contact_val;"volumeFraction" eta_input eta; "Kappa" k_input k; "gamma_0" gamma_0_input gamma_0; "Radius used in S(q)" radius radius; "polarRadius" radius_polar radius_polar; "equatorialRadius" radius_equatorial radius_equatorial; "SLD" SLD SLD; "SLDsolvent" SLD_solv SLD_solv; "Background" background background];
function s_vals = renorm(eta, k, gamma_0,radius)

    s_vals = zeros(40,2);
    s_vals(:,1) = [1:-0.025:0.025];
    eta_orig = eta;
    k_orig = k;
    gamma_0_orig = gamma_0;
    for j = 1:40
        eta = eta_orig/s_vals(j,1)^3;
        k = k_orig/s_vals(j,1);
        gamma_0 = gamma_0_orig*s_vals(j,1);
    diameter = 2*radius;
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



