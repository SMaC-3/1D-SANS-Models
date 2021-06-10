clc
clear
clf
Ang = char(197);

% 2021
% Code constructed by Ashley P. Williams, SMaCLab, Monash University, Australia
%--------------------------------------------------------------------------
%                            Spherical Form Factor 
%--------------------------------------------------------------------------
%
% This model is the sphere model to fit small–angle neutron scattering 
% data of spherical particles.
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

filename = "YourFileNameHere";
load(filename);
data = YourFileNameHere;

%--------------------------------------------------------------------------
%                          Sphere parameters
%--------------------------------------------------------------------------

scale = 1;
radius = 20;
SLD = 1;
SLD_solv = 6.3;
background = 0.001;
PD = 0;              %Schulz polydispersity              

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 7e-1;
xlower = 1e-3;
yupper = 1e3;
ylower = 1e-3;

%--------------------------------------------------------------------------
%                          Spherical calculation
%--------------------------------------------------------------------------

volume = (4*pi/3)*radius^3;
s = (SLD-SLD_solv)*volume;

I_spherical = zeros(size(data,1),2);
Form = @(r,q) (3*(sin(q*r)/(q*r)-cos(q*r))/(q*r)^2)^2;

%--------------------------------------------------------------------------
%                       Schulz-Zimm Polydispersity
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
    for h = 1:(size(data(:,1),1))
        for i = 1:size(rvals,2)
        form(i,h)=dist_r(1,i)*Form(rvals(1,i),data(h,1));
        end
    end
    form = mean(form/sum(dist_r/(size(dist_r,2)))); 
else
    form = [];
    for h = 1:(size(data(:,1),1))
        form(1,h) = Form(radius,data(h,1));
    end
end

%--------------------------------------------------------------------------
%                          Final file formatting
%--------------------------------------------------------------------------

I_spherical(:,2) = (scale/volume)*(1e-4)*s^2.*form+background;
I_spherical(:,1) = data(:,1); %[q,I]
parameters = ["Scale" scale; "Radius" radius; "SLD" SLD; "SLDsolvent" SLD_solv; "Background" background; "Polydispersity" PD];

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
legend("Data","Sphere Model")

