%% saveCaravandata - save various Caravan datasets as struct files
%
%   This script loads various Caravan datasets and saves them as MATLAB
%   struct files for easy use in MATLAB. The timeseries are loaded based on
%   the provided period. Note that this script can be slow and requires
%   sufficient RAM. 
%   
%   References
%   Kratzert, F., Nearing, G., Addor, N. et al. Caravan - A global community
%   dataset for large-sample hydrology. Sci Data 10, 61 (2023).
%   
%   Copyright (C) 2023

close all
% clear all
clc

%% Data Locations and directories
% After downloading the Caravan datasets and upzip the file, you can
% download the code and set in the first layer of Caravan folder.
% The structure of folder can be seen as follows:
% Caravan/Caravan/attributes, code, licenses, shapefiles, timeseries/...
% The first Caravan folder will store all function scripts and this script.

% If we want to use a different folder structure, we have to adjust the
% paths. 
% We need to add the Caravan repository to the MATLAB path, so that we can
% work with relative paths. 
mydir = 'D:/2.dataset/14.Hydrogauge/Caravan/'; % Change to your path here
addpath(genpath(mydir));

% The resulting struct files will be stored in a folder named "Data". If
% this folder does not exist yet, we have to create it. 
if ~(exist(strcat(mydir, '/Data')) == 7)
    mkdir(strcat(mydir, '/Data'))
end

%% Caravan_camels
%   The Caravan_camels dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camels_0' for gauges with 9 digits in attributes_caravan_camels.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together. 
%save_struct = true;
%saveCaravanstruct_camels(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_camelsaus
%   The Caravan_camels dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelsaus_0' for gauges with 9 digits in attributes_caravan_camelsaus.csv.
%   gauge_id in camelsaus contain some characters. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together. 
%   BE CAREFUL: here we use Python (jupyter notebook) to firstly rename the
%   camelsaus timeseries csv files. We delete the Capital letters in the
%   file name. And you may want to rename 'camelsaus_006005.csv' file
%   separately since we numeric the gauge_id. 
%save_struct = true;
%saveCaravanstruct_camelsaus(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_camelsbr
%   The Caravan_camelsbr dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelsbr_' for gauge_id in attributes_caravan_camelsbr.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_camelsbr(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_camelscl
%   The Caravan_camelscl dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelscl_' for gauge_id in attributes_caravan_camelscl.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_camelscl(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_camelsgb
%   The Caravan_camelsgb dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelsgb_' for gauge_id in attributes_caravan_camelsgb.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_camelsgb(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_hysets
%   The Caravan_hysets is the largest subdatasets in the Caravan dataset,
%   you can dowload from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we generate the 'gauge_index' which is
%   numeric to load in MATLAB conveniently. 
%   You may want to see Python (jupyter notebook) code firstly, and then
%   run this code.
save_struct = true;
saveCaravanstruct_hysets(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_lamah
%   The Caravan_lamah dataset is included in the Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7540792
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'lamah_' for gauge_id in attributes_caravan_lamah.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_lamah(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_camelsdk
%   The Caravan_camels dataset is an extension of Caravan dataset, you can
%   download from: https://zenodo.org/record/7396466
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   You may concerned about the streamflow missing data in
%   Caravan_camelsdk, at least be careful when your analysis require avoid
%   null value.
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelsdk_' for gauge_id in attributes_caravan_camelsdk.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_camelsdk(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------

%% Caravan_il
%   The Caravan_il dataset is an extension of Caravan dataset, you can
%   download from: https://doi.org/10.5281/zenodo.7758516
%   For timeseries data, there are two types, csv and netcdf, we only use the
%   csv files here. 
%   You may concerned about the streamflow missing data in
%   Caravan_camelsdk, at least be careful when your analysis require avoid
%   null value.
%   Caravan_il has longer streamflow timeseries, from 1970-2021. Therefore,
%   it has different length with other Caravan subdatasets. And it also
%   contains many zero-value and null value here.
%   Due to the format of gauge_id, we conduct pre-processing to delete
%   'camelsdk_' for gauge_id in attributes_caravan_camelsdk.csv. 
%   We recommend you to run them separately, because it may require
%   sufficient RAM for your computer if you generate them together.
%save_struct = true;
%saveCaravanstruct_il(save_struct);
%   ----------------------------------------------------------------------
%   ----------------------------------------------------------------------