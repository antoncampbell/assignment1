%% Assignment 1
% By: Anton Campbell, 100975168
%% Question 1 - Electron Modeling
% A simple model was created that ran for 10,000 electrons. The electrons
% were all set to have the same speed (the thermal velocity) but in
% different directions.
%
% The thermal velocity was calculated using the equation:
% 
% $$v_{th}=\sqrt{\frac{2kT}{m_{n}}}$$
%
% The mean free path was calculated using the equation:
% 
% $$MFP=v_{th}\tau _{mn}$$
%
% The simulation was run:
electron_box_3modes(1)

%%
% a)
% The theoretical thermal velocity (as calculated above) is $187019.126$ m/s.
%
% b)
% The theoretical mean free path (as calculated above) is $$3.740\times
% 10^{-8}$ m.
%
% c)
%
% i)
% As seen in the Figure 1, the particles move in straight lines. They bounce off the top and
% bottom. The x-direction has a periodic boundary as electrons pass from
% one side to the other.
%
% ii)
% Figure 2 shows the temperature of the system over time. The temperature remains constant over all time.
% This makes sense because the velocities of the electrons are fixed.
%


%% Question 2 - Collisions with Mean Free Path
% Some changes where made to the initial simulation. Each particle was assigned a random velocity based on a normal
% distribution and the thermal velocity. The particles scatter randomly and
% are reassigned new velocities based on a normal
% distribution and the thermal velocity.
%
% The experimental mean time between was calculated was calculated using the equation:
% 
% $$\tau _{mn}=\frac{\sum_{i}^{\#collisions}(\#steps\ since\ last\ collsion)(time\ step)}{\#collisions}$$
%
% The thermal speed was calculated using the equation:
% 
% $$MFP=\frac{\sum_{i}^{\#collisions}(speed)(\#steps\ since\ last\ collsion)(time\ step)}{\#collisions}$$
% 
%
% The simulation was run:
electron_box_3modes(2)

%%
% a)
% Figure 3 shows a histogram of the different speeds. This histogram is a Maxwell-Boltzman distribution.
%
% b)
% Figure 4 shows some particle trajectories with the scattering functionality. The particles scatter after different distances traveled.
%
% c)
% Figure 5 shows the temperature of the system over time when the scattering functionality is present. The temperature fluctates slightly over time. It tends to 
% be close to the set temperature of 300K. The flucuation occurs due to the
% random scattering which does not account for conservation of energy.
% 
% d) The measured value for the time between collisions (above) is quite closed to
% the set value of $2\times 10^{-13}$ s.
% The measured MFP (above) is somewhat off from the actual value of $3.740\times 10^{-8}$ m
% (see question 1). This discrepancy may be due to approximations and
% random changes of the simulation.

%% Question 3 - Enhancements
% Futher changes where made to the simulation. A 'bottle neck' was added. The bottom boundary is diffusive and the top one is specular. 
%
% The simulation was run:
electron_box_3modes(3)

%%
% a)
% Figure 6 shows some particle trajectories with the added 'bottle neck'.
% The particles bounce with the same angle for the top boundary (specular).
% The particles bounce with a different angle off the bottom boundary
% (diffusive).
%
% c)
% Figure 7 shows a electron density map. The electrons appear to spread out
% in the area outside the boundaries.
%
% d)
% Figure 8 shows a temperature map. Due to the size of the bins, the
% local temperatures vary greatly. Some electrons happen to be travelling
% quite fast resulting in a much higher local temperature.
% For bins without any electrons, the temperature was set to be zero
% although there is not really any local temperature.
% 

