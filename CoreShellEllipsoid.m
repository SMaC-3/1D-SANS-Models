clc
clear
clf
Ang = char(197);

% 2022
% Code constructed by Ashley P. Williams, SMaCLab, Monash University, Australia
%--------------------------------------------------------------------------
%                       Core-Shell Ellipsoid Form Factor 
%--------------------------------------------------------------------------
%
% This model is the core-shell ellipsoid model to fit small–angle neutron 
% scattering data of core-shell ellipsoidal particles.
%
% Load your file below in the appropriate section, where your data file is
% of the form of 3 column vectors: q, I, err. Also ensure that you cut any
% points with NaN values, to minimise any issues.
%
% To save the fit, open the vector "I_ellipsoid" from the workspace and
% copy the contents into a text file.
%
% To save the parameters, open the vector "parameters" and save as desired.
% The parameters are saved as a 2 column vector: parameters, values
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
%                          Ellipsoid parameters
%--------------------------------------------------------------------------

scale = 1;
radius_polar = 20;
radius_equatorial = 30;
thickness = 10;            %Shell thickness, assumed to be equal everywhere.
SLD_core = 1;
SLD_shell = 3;
SLD_solv = 6.3;
background = 0.001;
PD = 0.16;                  % Dispersity in the polar radius.

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 4e-1;
xlower = 1.5e-3;
yupper = 2e1;
ylower = 1e-2;

%--------------------------------------------------------------------------
%                          Ellipsoid calculation
%--------------------------------------------------------------------------

contrast_core = (SLD_core-SLD_shell);
contrast_whole = (SLD_shell-SLD_solv);
volume_core = @ (r) (4*pi/3)*r*radius_equatorial^2; % r = polar radius
volume_whole = @(r) (4*pi/3)*(r+thickness)*(radius_equatorial+thickness)^2;% r = polar radius
v_square_minus_one_core = @(r) (r/radius_equatorial)^2-1; % r = polar radius
v_square_minus_one_shell = @(r) ((r+thickness)/(radius_equatorial+thickness))^2-1; % r = polar radius
zm = 0.5;
zb = 0.5;

gauss_N = 76; %According to /lib/Gauss76.c
gauss_W = [.00126779163408536,.00294910295364247,.00462793522803742,.00629918049732845,.00795984747723973,.00960710541471375,.0112381685696677,.0128502838475101,.0144407317482767,.0160068299122486,.0175459372914742,.0190554584671906,.020532847967908,.0219756145344162,.0233813253070112,.0247476099206597,.026072164497986,.0273527555318275,.028587223650054,.029773487255905,.0309095460374916,.0319934843404216,.0330234743977917,.0339977794120564,.0349147564835508,.0357728593807139,.0365706411473296,.0373067565423816,.0379799643084053,.0385891292645067,.0391332242205184,.0396113317090621,.0400226455325968,.040366472122844,.0406422317102947,.0408494593018285,.040987805464794,.0410570369162294,.0410570369162294,.040987805464794,.0408494593018285,.0406422317102947,.040366472122844,.0400226455325968,.0396113317090621,.0391332242205184,.0385891292645067,.0379799643084053,.0373067565423816,.0365706411473296,.0357728593807139,.0349147564835508,.0339977794120564,.0330234743977917,.0319934843404216,.0309095460374916,.029773487255905,.028587223650054,.0273527555318275,.026072164497986,.0247476099206597,.0233813253070112,.0219756145344162,.020532847967908,.0190554584671906,.0175459372914742,.0160068299122486,.0144407317482767,.0128502838475101,.0112381685696677,.00960710541471375,.00795984747723973,.00629918049732845,.00462793522803742,.00294910295364247,.00126779163408536];
gauss_Z = [-.999505948362153,-.997397786355355,-.993608772723527,-.988144453359837,-.981013938975656,-.972229228520377,-.961805126758768,-.949759207710896,-.936111781934811,-.92088586125215,-.904107119545567,-.885803849292083,-.866006913771982,-.844749694983342,-.822068037328975,-.7980001871612,-.77258672828181,-.74587051350361,-.717896592387704,-.688712135277641,-.658366353758143,-.626910417672267,-.594397368836793,-.560882031601237,-.526420920401243,-.491072144462194,-.454895309813726,-.417951418780327,-.380302767117504,-.342012838966962,-.303146199807908,-.263768387584994,-.223945802196474,-.183745593528914,-.143235548227268,-.102483975391227,-.0615595913906112,-.0205314039939986,.0205314039939986,.0615595913906112,.102483975391227,.143235548227268,.183745593528914,.223945802196474,.263768387584994,.303146199807908,.342012838966962,.380302767117504,.417951418780327,.454895309813726,.491072144462194,.526420920401243,.560882031601237,.594397368836793,.626910417672267,.658366353758143,.688712135277641,.717896592387704,.74587051350361,.77258672828181,.7980001871612,.822068037328975,.844749694983342,.866006913771982,.885803849292083,.904107119545567,.92088586125215,.936111781934811,.949759207710896,.961805126758768,.972229228520377,.981013938975656,.988144453359837,.993608772723527,.997397786355355,.999505948362153];

if PD ~= 0
    Z = (1-PD^2)/PD^2;
    sz = @(x) (((Z+1)/radius_polar)^(Z+1))*((x)^(Z))*exp(-((Z+1)/radius_polar)*x)...
        /(gamma(Z+1)); %Schulz distribution function
    N_steps = 80;
    dist_r = [];
    rvals = ((radius_polar-3*PD*radius_polar):(6*PD*radius_polar)/(N_steps):(radius_polar+3*PD*radius_polar)); %Number of points for PD calculation.
    for j = 1:size(rvals,2)
        dist_r(j) = sz(rvals(j)); 
    end
    dist_r = dist_r/sum(dist_r); %Normalise the distribution weights such that sum = 1.

    form = zeros(size(rvals,2),size(data,1));
    for k = 1:size(rvals,2)
        for j = 1:size(data,1)
            for i = 1:gauss_N
               u = gauss_Z(i)*zm +zb;
               r_core = (radius_equatorial)*sqrt(1+u^2*v_square_minus_one_core(rvals(k)));
               r_shell = (radius_equatorial+thickness)*sqrt(1+u^2*v_square_minus_one_shell(rvals(k)));
               f_core = contrast_core*volume_core(rvals(k))*3*(sin(data(j,1)*r_core)-(data(j,1)*r_core*cos(data(j,1)*r_core)))/(data(j,1)*r_core)^3;
               f_shell = contrast_whole*volume_whole(rvals(k))*3*(sin(data(j,1)*r_shell)-(data(j,1)*r_shell*cos(data(j,1)*r_shell)))/(data(j,1)*r_shell)^3;
               f = f_core+f_shell;
               form(k,j) = form(k,j)+gauss_W(i)*f^2; 
            end
            form(k,j) = zm.*form(k,j); %/sum(gaussian weights)
        end 
        form(k,:) = dist_r(k).*form(k,:); %weighting per radius
    end

     for n = 1:size(rvals,2)
         vol(n) = dist_r(n)*volume_whole(rvals(n));
     end

     form = sum(form)/sum(vol); 
else
    form = zeros(size(data,1),1);
    for j = 1:size(data,1)
        for i = 1:gauss_N
           u = gauss_Z(i)*zm +zb;
           r_core = radius_equatorial*sqrt(1+u^2*v_square_minus_one_core(radius_polar));
           r_shell = (radius_equatorial+thickness)*sqrt(1+u^2*v_square_minus_one_shell(radius_polar));
           f_core = contrast_core*volume_core(radius_polar)*(sin(data(j,1)*r_core)-(data(j,1)*r_core*cos(data(j,1)*r_core)))/(data(j,1)*r_core)^3;
           f_shell = contrast_whole*volume_whole(radius_polar)*(sin(data(j,1)*r_shell)-(data(j,1)*r_shell*cos(data(j,1)*r_shell)))/(data(j,1)*r_shell)^3;
           f = 3*(f_core+f_shell);
           form(j,1) = form(j,1)+gauss_W(i)*f^2; 
        end
        form(j,1) = zm.*form(j,1);
    end
end

if PD ~= 0
    I_ellipsoid(:,2) = (scale)*(1e-4).*form+background; %F^2
else
    I_ellipsoid(:,2) = (scale)/volume_whole(radius_polar)*(1e-4).*form+background; %F^2
end
I_ellipsoid(:,1) = data(:,1); %[q,I]

parameters = ["Scale" scale; "polarRadius" radius_polar; "equatorialRadius" radius_equatorial; "Thickness" thickness; "SLDcore" SLD_core; "SLDshell" SLD_shell; "SLDsolvent" SLD_solv; "Background" background];

figure(1)
hold on
box on
xlabel("{\it q} ("+Ang+"^{-1})")
ylabel("Intensity (cm^{-1})")
set(gca,'Xscale', 'log','Yscale', 'log','fontweight','bold','Linewidth',2,'fontSize',14)
xlim([xlower, xupper])
ylim([ylower yupper])
errorbar(data(:,1),data(:,2),data(:,3),data(:,3),'o','MarkerFaceColor','auto','MarkerSize', 6,'Color','Black')
plot(I_ellipsoid(:,1),I_ellipsoid(:,2),'Color',[0 0 1],'Linewidth',3)
legend("Data","Core-Shell Ellipsoid")

