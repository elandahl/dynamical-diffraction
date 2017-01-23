% Program to save material properties into a matlab structure file called sample.dat
% Written by Eric Landahl, 12.22.2016
% Last revised by EL 12.28.2016
% Most properties from reference 'ioffe' : http://www.ioffe.ru/SVA/NSM/Semicond/
%
% Some hints:
%
% After running this program sample.dat is created and the variable sample is loaded
% After clearing memory, reload the data using:
% load sample.dat
%
% To find the number of a sample:
% index = find(strcmp({sample.name}, 'GaAs')==1)
%
% To assign the .var variable name to the value .val :
% v = genvarname(sample(1).bandGap.var); eval([v '=sample(1).bandGap.val']);
%
%
%% Program to load a particular sample's variables into memory
%% This program is also distributed as test_loader.m
%%
%clear all;
%load sample.dat;ID = find(strcmp({sample.name}, 'Si')==1); % ID is code for sample
%num_fields = length(fieldnames(sample(ID))); % Determine the number of fields
%for i = 2:num_fields % Skip first field (it is always the name)
%  vname=genvarname(getfield(getfield(sample(ID),fieldnames(sample(ID)){i}),'var'));
%  vval = getfield(getfield(sample(ID),fieldnames(sample(ID)){i}),'val');
%  eval([vname ' = vval;']);
%end


clear sample;


% Silicon
sample(1).name = 'Si';

  sample(1).bandGap.val = 1.5;
  sample(1).bandGap.unit = 'eV';
  sample(1).bandGap.var = 'Eg';
  sample(1).bandGap.ref = 'ioffe';

  sample(1).molarDensity.val = 5E22;
  sample(1).molarDensity.unit = 'atoms/cm^3';
  sample(1).molarDensity.var = 'rho_n';
  sample(1).molarDensity.ref = 'ioffe';
  
  sample(1).AugerN.val = 1.1E-30;
  sample(1).AugerN.unit = 'cm^6/s';
  sample(1).AugerN.var = 'Cn';
  sample(1).AugerN.ref = 'ioffe';
  
  sample(1).AugerP.val = 3E-31;
  sample(1).AugerP.unit = 'cm^6/s';
  sample(1).AugerP.var = 'Cp';
  sample(1).AugerP.ref = 'ioffe';
  
  sample(1).DebyeTemp.val = 640;
  sample(1).DebyeTemp.unit = 'K';
  sample(1).DebyeTemp.var = 'Theta_D';
  sample(1).DebyeTemp.ref = 'ioffe';

  sample(1).massDensity.val = 2.329;
  sample(1).massDensity.unit = 'g/cm^3';
  sample(1).massDensity.var = 'rho_m';
  sample(1).massDensity.ref = 'ioffe';
  
  sample(1).opticalPhonon.val = 0.063;
  sample(1).opticalPhonon.unit = 'eV';
  sample(1).opticalPhonon.var = 'E_O';
  sample(1).opticalPhonon.ref = 'ioffe'; 

  sample(1).diffusionElectrons.val = 36;
  sample(1).diffusionElectrons.unit = 'cm^2/s';
  sample(1).diffusionElectrons.var = 'D_e';
  sample(1).diffusionElectrons.ref = 'ioffe'; 
  
  sample(1).diffusionHoles.val = 12;
  sample(1).diffusionHoles.unit = 'cm^2/s';
  sample(1).diffusionHoles.var = 'D_h';
  sample(1).diffusionHoles.ref = 'ioffe';
  
  sample(1).radiativeRecombination.val = 1.1E-14;
  sample(1).radiativeRecombination.unit = 'cm^3/s';
  sample(1).radiativeRecombination.var = 'B';
  sample(1).radiativeRecombination.ref = 'ioffe';

  sample(1).bulkModulus.val = 9.8E11;
  sample(1).bulkModulus.unit = 'dyn/cm^2';
  sample(1).bulkModulus.var = 'K';
  sample(1).bulkModulus.ref = 'ioffe';
  
  sample(1).specificHeat.val = 0.7;
  sample(1).specificHeat.unit = 'J/(g K)';
  sample(1).specificHeat.var = 'C';
  sample(1).specificHeat.ref = 'ioffe';
  
  sample(1).thermalConductivity.val = 1.3;
  sample(1).thermalConductivity.unit = 'W/(cm K)';
  sample(1).thermalConductivity.var = 'kappa_t';
  sample(1).thermalConductivity.ref = 'ioffe';  

  sample(1).thermalDiffusion.val = 0.8;
  sample(1).thermalDiffusion.unit = 'cm^2/s';
  sample(1).thermalDiffusion.var = 'D_t';
  sample(1).thermalDiffusion.ref = 'ioffe';
  
  sample(1).thermalExpansion.val = 2.6E-6;
  sample(1).thermalExpansion.unit = '1/K';
  sample(1).thermalExpansion.var = 'alpha_t';
  sample(1).thermalExpansion.ref = 'ioffe';
  
  sample(1).PoissonRatio.val = 0.28;
  sample(1).PoissonRatio.unit = '';
  sample(1).PoissonRatio.var = 'nu';
  sample(1).PoissonRatio.ref = 'ioffe';
  
  sample(1).soundLongitudinal.val = 8.43;
  sample(1).soundLongitudinal.unit = 'um/ns';
  sample(1).soundLongitudinal.var = 'v_L';
  sample(1).soundLongitudinal.ref = 'ioffe';
  
  sample(1).soundTransverse.val = 5.84;
  sample(1).soundTransverse.unit = 'um/ns';
  sample(1).soundTransverse.var = 'v_T';
  sample(1).soundTransverse.ref = 'ioffe';

% Gallium Arsenide
sample(2).name = 'GaAs';

  sample(2).bandGap.val = 1.42;
  sample(2).bandGap.unit = 'eV';
  sample(2).bandGap.var = 'Eg';
  sample(2).bandGap.ref = 'ioffe';

  sample(2).molarDensity.val = 0;
  sample(2).molarDensity.unit = 'atoms/cm^3';
  sample(2).molarDensity.var = 'rho_n';
  sample(2).molarDensity.ref = 'ioffe';
  
  sample(2).AugerN.val = 0;
  sample(2).AugerN.unit = 'cm^6/s';
  sample(2).AugerN.var = 'Cn';
  sample(2).AugerN.ref = 'ioffe';
  
  sample(2).AugerP.val = 0;
  sample(2).AugerP.unit = 'cm^6/s';
  sample(2).AugerP.var = 'Cp';
  sample(2).AugerP.ref = 'ioffe';
  
  sample(2).DebyeTemp.val = 0;
  sample(2).DebyeTemp.unit = 'K';
  sample(2).DebyeTemp.var = 'Theta_D';
  sample(2).DebyeTemp.ref = 'ioffe';

  sample(2).massDensity.val = 5.32;
  sample(2).massDensity.unit = 'g/cm^3';
  sample(2).massDensity.var = 'rho_m';
  sample(2).massDensity.ref = 'ioffe';
  
  sample(2).opticalPhonon.val = 0;
  sample(2).opticalPhonon.unit = 'eV';
  sample(2).opticalPhonon.var = 'E_O';
  sample(2).opticalPhonon.ref = 'ioffe'; 

  sample(2).diffusionElectrons.val = 0;
  sample(2).diffusionElectrons.unit = 'cm^2/s';
  sample(2).diffusionElectrons.var = 'D_e';
  sample(2).diffusionElectrons.ref = 'ioffe'; 
  
  sample(2).diffusionHoles.val = 0;
  sample(2).diffusionHoles.unit = 'cm^2/s';
  sample(2).diffusionHoles.var = 'D_h';
  sample(2).diffusionHoles.ref = 'ioffe';
  
  sample(2).radiativeRecombination.val = 0;
  sample(2).radiativeRecombination.unit = 'cm^3/s';
  sample(2).radiativeRecombination.var = 'B';
  sample(2).radiativeRecombination.ref = 'ioffe';

  sample(2).bulkModulus.val = 0;
  sample(2).bulkModulus.unit = 'dyn/cm^2';
  sample(2).bulkModulus.var = 'K';
  sample(2).bulkModulus.ref = 'ioffe';
  
  sample(2).specificHeat.val = 3.30E-01;
  sample(2).specificHeat.unit = 'J/(g K)';
  sample(2).specificHeat.var = 'C';
  sample(2).specificHeat.ref = 'ioffe';
  
  sample(2).thermalConductivity.val = 0.55;
  sample(2).thermalConductivity.unit = 'W/(cm K)';
  sample(2).thermalConductivity.var = 'kappa_t';
  sample(2).thermalConductivity.ref = 'ioffe';  

  sample(2).thermalDiffusion.val = 0.31;
  sample(2).thermalDiffusion.unit = 'cm^2/s';
  sample(2).thermalDiffusion.var = 'D_t';
  sample(2).thermalDiffusion.ref = 'ioffe';
  
  sample(2).thermalExpansion.val = 5.73E-06;
  sample(2).thermalExpansion.unit = '1/K';
  sample(2).thermalExpansion.var = 'alpha_t';
  sample(2).thermalExpansion.ref = 'ioffe';
  
  sample(2).PoissonRatio.val = 0.31;
  sample(2).PoissonRatio.unit = '';
  sample(2).PoissonRatio.var = 'nu';
  sample(2).PoissonRatio.ref = 'ioffe';
  
  sample(2).soundLongitudinal.val = 4.73;
  sample(2).soundLongitudinal.unit = 'um/ns';
  sample(2).soundLongitudinal.var = 'v_L';
  sample(2).soundLongitudinal.ref = 'ioffe';
  
  sample(2).soundTransverse.val = 3.35;
  sample(2).soundTransverse.unit = 'um/ns';
  sample(2).soundTransverse.var = 'v_T';
  sample(2).soundTransverse.ref = 'ioffe';

% Germanium
sample(3).name = 'Ge';

  sample(3).bandGap.val = 0.661;
  sample(3).bandGap.unit = 'eV';
  sample(3).bandGap.var = 'Eg';
  sample(3).bandGap.ref = 'ioffe';

  sample(3).molarDensity.val = 0;
  sample(3).molarDensity.unit = 'atoms/cm^3';
  sample(3).molarDensity.var = 'rho_n';
  sample(3).molarDensity.ref = 'ioffe';
  
  sample(3).AugerN.val = 0;
  sample(3).AugerN.unit = 'cm^6/s';
  sample(3).AugerN.var = 'Cn';
  sample(3).AugerN.ref = 'ioffe';
  
  sample(3).AugerP.val = 0;
  sample(3).AugerP.unit = 'cm^6/s';
  sample(3).AugerP.var = 'Cp';
  sample(3).AugerP.ref = 'ioffe';
  
  sample(3).DebyeTemp.val = 0;
  sample(3).DebyeTemp.unit = 'K';
  sample(3).DebyeTemp.var = 'Theta_D';
  sample(3).DebyeTemp.ref = 'ioffe';

  sample(3).massDensity.val = 5.3234;
  sample(3).massDensity.unit = 'g/cm^3';
  sample(3).massDensity.var = 'rho_m';
  sample(3).massDensity.ref = 'ioffe';
  
  sample(3).opticalPhonon.val = 0;
  sample(3).opticalPhonon.unit = 'eV';
  sample(3).opticalPhonon.var = 'E_O';
  sample(3).opticalPhonon.ref = 'ioffe'; 

  sample(3).diffusionElectrons.val = 0;
  sample(3).diffusionElectrons.unit = 'cm^2/s';
  sample(3).diffusionElectrons.var = 'D_e';
  sample(3).diffusionElectrons.ref = 'ioffe'; 
  
  sample(3).diffusionHoles.val = 0;
  sample(3).diffusionHoles.unit = 'cm^2/s';
  sample(3).diffusionHoles.var = 'D_h';
  sample(3).diffusionHoles.ref = 'ioffe';
  
  sample(3).radiativeRecombination.val = 0;
  sample(3).radiativeRecombination.unit = 'cm^3/s';
  sample(3).radiativeRecombination.var = 'B';
  sample(3).radiativeRecombination.ref = 'ioffe';

  sample(3).bulkModulus.val = 0;
  sample(3).bulkModulus.unit = 'dyn/cm^2';
  sample(3).bulkModulus.var = 'K';
  sample(3).bulkModulus.ref = 'ioffe';
  
  sample(3).specificHeat.val = 0.31;
  sample(3).specificHeat.unit = 'J/(g K)';
  sample(3).specificHeat.var = 'C';
  sample(3).specificHeat.ref = 'ioffe';
  
  sample(3).thermalConductivity.val = 0.58;
  sample(3).thermalConductivity.unit = 'W/(cm K)';
  sample(3).thermalConductivity.var = 'kappa_t';
  sample(3).thermalConductivity.ref = 'ioffe';  

  sample(3).thermalDiffusion.val = 0.36;
  sample(3).thermalDiffusion.unit = 'cm^2/s';
  sample(3).thermalDiffusion.var = 'D_t';
  sample(3).thermalDiffusion.ref = 'ioffe';
  
  sample(3).thermalExpansion.val = 5.90E-06;
  sample(3).thermalExpansion.unit = '1/K';
  sample(3).thermalExpansion.var = 'alpha_t';
  sample(3).thermalExpansion.ref = 'ioffe';
  
  sample(3).PoissonRatio.val = 0.26;
  sample(3).PoissonRatio.unit = '';
  sample(3).PoissonRatio.var = 'nu';
  sample(3).PoissonRatio.ref = 'ioffe';
  
  sample(3).soundLongitudinal.val = 4.87;
  sample(3).soundLongitudinal.unit = 'um/ns';
  sample(3).soundLongitudinal.var = 'v_L';
  sample(3).soundLongitudinal.ref = 'ioffe';
  
  sample(3).soundTransverse.val = 3.57;
  sample(3).soundTransverse.unit = 'um/ns';
  sample(3).soundTransverse.var = 'v_T';
  sample(3).soundTransverse.ref = 'ioffe';

% Indium Antimonide
sample(4).name = 'InSb';

  sample(4).bandGap.val = 0.17;
  sample(4).bandGap.unit = 'eV';
  sample(4).bandGap.var = 'Eg';
  sample(4).bandGap.ref = 'ioffe';

  sample(4).molarDensity.val = 0;
  sample(4).molarDensity.unit = 'atoms/cm^3';
  sample(4).molarDensity.var = 'rho_n';
  sample(4).molarDensity.ref = 'ioffe';
  
  sample(4).AugerN.val = 0;
  sample(4).AugerN.unit = 'cm^6/s';
  sample(4).AugerN.var = 'Cn';
  sample(4).AugerN.ref = 'ioffe';
  
  sample(4).AugerP.val = 0;
  sample(4).AugerP.unit = 'cm^6/s';
  sample(4).AugerP.var = 'Cp';
  sample(4).AugerP.ref = 'ioffe';
  
  sample(4).DebyeTemp.val = 0;
  sample(4).DebyeTemp.unit = 'K';
  sample(4).DebyeTemp.var = 'Theta_D';
  sample(4).DebyeTemp.ref = 'ioffe';

  sample(4).massDensity.val = 5.77;
  sample(4).massDensity.unit = 'g/cm^3';
  sample(4).massDensity.var = 'rho_m';
  sample(4).massDensity.ref = 'ioffe';
  
  sample(4).opticalPhonon.val = 0;
  sample(4).opticalPhonon.unit = 'eV';
  sample(4).opticalPhonon.var = 'E_O';
  sample(4).opticalPhonon.ref = 'ioffe'; 

  sample(4).diffusionElectrons.val = 0;
  sample(4).diffusionElectrons.unit = 'cm^2/s';
  sample(4).diffusionElectrons.var = 'D_e';
  sample(4).diffusionElectrons.ref = 'ioffe'; 
  
  sample(4).diffusionHoles.val = 0;
  sample(4).diffusionHoles.unit = 'cm^2/s';
  sample(4).diffusionHoles.var = 'D_h';
  sample(4).diffusionHoles.ref = 'ioffe';
  
  sample(4).radiativeRecombination.val = 0;
  sample(4).radiativeRecombination.unit = 'cm^3/s';
  sample(4).radiativeRecombination.var = 'B';
  sample(4).radiativeRecombination.ref = 'ioffe';

  sample(4).bulkModulus.val = 0;
  sample(4).bulkModulus.unit = 'dyn/cm^2';
  sample(4).bulkModulus.var = 'K';
  sample(4).bulkModulus.ref = 'ioffe';
  
  sample(4).specificHeat.val = 0.2;
  sample(4).specificHeat.unit = 'J/(g K)';
  sample(4).specificHeat.var = 'C';
  sample(4).specificHeat.ref = 'ioffe';
  
  sample(4).thermalConductivity.val = 0.18;
  sample(4).thermalConductivity.unit = 'W/(cm K)';
  sample(4).thermalConductivity.var = 'kappa_t';
  sample(4).thermalConductivity.ref = 'ioffe';  

  sample(4).thermalDiffusion.val = 0.16;
  sample(4).thermalDiffusion.unit = 'cm^2/s';
  sample(4).thermalDiffusion.var = 'D_t';
  sample(4).thermalDiffusion.ref = 'ioffe';
  
  sample(4).thermalExpansion.val = 5.37E-06;
  sample(4).thermalExpansion.unit = '1/K';
  sample(4).thermalExpansion.var = 'alpha_t';
  sample(4).thermalExpansion.ref = 'ioffe';
  
  sample(4).PoissonRatio.val = 0.35;
  sample(4).PoissonRatio.unit = '';
  sample(4).PoissonRatio.var = 'nu';
  sample(4).PoissonRatio.ref = 'ioffe';
  
  sample(4).soundLongitudinal.val = 3.4;
  sample(4).soundLongitudinal.unit = 'um/ns';
  sample(4).soundLongitudinal.var = 'v_L';
  sample(4).soundLongitudinal.ref = 'ioffe';
  
  sample(4).soundTransverse.val = 2.29;
  sample(4).soundTransverse.unit = 'um/ns';
  sample(4).soundTransverse.var = 'v_T';
  sample(4).soundTransverse.ref = 'ioffe';

%save sample.dat sample




