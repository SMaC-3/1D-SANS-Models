clc
clear
clf
Ang = char(197);

% 2022
% Code constructed by Ashley P. Williams, SMaCLab, Monash University, Australia
%--------------------------------------------------------------------------
%                      Spherical Core-Shell Form Factor 
%--------------------------------------------------------------------------
%
% This model is the spherical core-shell model to fit small–angle neutron 
% scattering data of spherical core-shell particles.
%
% Load your file below in the appropriate section, where your data file is
% of the form of 3 column vectors: q, I, err. Also ensure that you cut any
% points with NaN values, to minimise any issues.
%
% To save the fit, open the vector "I_spherical" from the workspace and
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
%                     Core-Shell Sphere parameters
%--------------------------------------------------------------------------

scale = 0.0145;
radius = 23.8;
thickness = 10;
SLD_core = 6.393;
SLD_shell = 1;
SLD_solv = 6.393;
background = 0.01;
PD = 0.16;              %Schulz polydispersity              

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 7e-1; 
xlower = 1e-3;
yupper = 5e2;
ylower = 1e-1;

%--------------------------------------------------------------------------
%                          Spherical calculation
%--------------------------------------------------------------------------

volume_core = @ (r) (4*pi/3)*r^3;
volume_whole = @ (r) (4*pi/3)*(r+thickness)^3;
cont_core = (SLD_core-SLD_shell);
cont_whole = (SLD_shell-SLD_solv);

I_spherical = zeros(size(data,1),2);
Form = @(r,q) ((sin(q*r)-q*r*cos(q*r))/(q*r)^3);

%--------------------------------------------------------------------------
%                       Schulz-Zimm Polydispersity (core radius)
%--------------------------------------------------------------------------

Z = (1-PD^2)/PD^2;
sz = @(x) ((Z+1)^(Z+1))*((x/radius)^(Z))*exp(-(Z+1)*(x/radius))...
    /(radius*gamma(Z+1)); %Schulz distribution function
N_steps = 80;
dist_r = [];
rvals = ((radius-3*PD*radius):(6*PD*radius)/N_steps:(radius+3*PD*radius)); %Number of points for PD calculation.
for j = 1:size(rvals,2)
    dist_r(1,j) = sz(rvals(1,j)); 
end
dist_r = dist_r/max(dist_r); %Normalise the distribution weights.

if PD ~= 0
    form = [];
    vol = [];
    for h = 1:(size(data(:,1),1))
        for i = 1:size(rvals,2)
           form(i,h)= dist_r(i)*(3*(volume_core(rvals(i))*cont_core*Form(rvals(i),data(h,1))+volume_whole(rvals(i))*cont_whole*Form(rvals(1,i)+thickness,data(h,1))))^2;
        end
    end
    
    for i = 1:size(rvals,2)
       vol(i) = dist_r(i)*(volume_whole(rvals(i)));
    end
  
    form = sum(form)/sum(vol);   
  
else
    form = [];
    for h = 1:(size(data(:,1),1))
        form(1,h) = (3*(volume_core(radius)*cont_core*Form(radius,data(h,1))+volume_whole(radius)*cont_whole*Form(radius+thickness,data(h,1))))^2;
    end
end

%--------------------------------------------------------------------------
%                          Final file formatting
%--------------------------------------------------------------------------
if PD == 0
    I_spherical(:,2) = (scale/volume_whole(radius))*(1e-4)*form+background;
else
    I_spherical(:,2) = scale*(1e-4)*form+background;
end

I_spherical(:,1) = data(:,1); %[q,I]
parameters = ["Scale" scale; "Radius" radius; "Thickness" thickness; "SLDcore" SLD_core; "SLDshell" SLD_shell; "SLDsolvent" SLD_solv; "Background" background; "Polydispersity" PD];

%--------------------------------------------------------------------------
%                               Plotting
%--------------------------------------------------------------------------

figure(1)
hold on
box on
xlabel("{\it q} ("+Ang+"^{-1})")
ylabel("Intensity (cm^{-1})")
set(gca,'Xscale', 'log','Yscale', 'log','fontweight','bold','Linewidth',2,'fontSize',14)
xlim([xlower, xupper])
ylim([ylower yupper])
errorbar(data(:,1),data(:,2),data(:,3),data(:,3),'o','MarkerFaceColor','auto','MarkerSize', 6,'Color','Black')
plot(I_spherical(:,1),I_spherical(:,2),'Color',[0 0 1],'Linewidth',3)
legend("Data","Core-Shell Sphere Model")

