clc
clear
clf
Ang = char(197);

% 2022
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
% Polydispersity here is given by a Schulz-Zimm distribution function of
% the form:
%
%       Schulz(x) = ((Z+1)/radius)^(Z+1) * x^Z 
%                   * exp[((Z+1)/radius)*x]/Gamma(Z+1)
%
% Where Z = (1-polydispersity^2)/polydispersity^2
%
% The final output Intensity is given as:
% 
% I(q) = scale*contrast^2*Sum[(w_i*F(q,r_i)^2)/Sum(w_i)] 
%                         /Sum[(w_i*volume(r_i)^2)/Sum(w_i)]+Background
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
%                          Sphere parameters
%--------------------------------------------------------------------------

scale = 0.0145;
radius = 23.5;
SLD = -0.392;
SLD_solv = 6.393;
background = 0.01;
PD = 0.16;              %Schulz polydispersity              

%--------------------------------------------------------------------------
%                          Plot preferences
%--------------------------------------------------------------------------

xupper = 5e-1;
xlower = 4e-3;
yupper = 1e1;
ylower = 3e-3;

%--------------------------------------------------------------------------
%                          Spherical calculation
%--------------------------------------------------------------------------

volume = @(r) (4*pi/3)*r^3;
s = (SLD-SLD_solv);

I_spherical = zeros(size(data,1),2);
Form = @(r,q) (1e-4)*(3*s*volume(r)*((sin(q*r)-q*r*cos(q*r))/(q*r)^3))^2;

%--------------------------------------------------------------------------
%                       Schulz-Zimm Polydispersity
%--------------------------------------------------------------------------

Z = (1-PD^2)/PD^2;
sz = @(x) (((Z+1)/radius)^(Z+1))*((x)^(Z))*exp(-((Z+1)/radius)*x)...
    /(gamma(Z+1)); %Schulz distribution function
N_steps = 80;
dist_r = [];
rvals = ((radius-3*PD*radius):(6*PD*radius)/(N_steps):(radius+3*PD*radius)); %Number of points for PD calculation.
for j = 1:size(rvals,2)
    dist_r(j) = sz(rvals(j)); 
end
dist_r = dist_r/sum(dist_r); %Normalise the distribution weights such that sum = 1.

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

%--------------------------------------------------------------------------
%                          Final file formatting
%--------------------------------------------------------------------------
if PD ~= 0
    I_spherical(:,2) = (scale).*form+background;
else
    I_spherical(:,2) = (scale/volume(radius)).*form+background;
end
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
legend("Data","Sphere")

