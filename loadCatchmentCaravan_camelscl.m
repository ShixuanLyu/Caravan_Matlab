function [Q, P, PET, Tmin, Tmax, Tmean, dewpoint_Tmin, dewpoint_Tmax, dewpoint_Tmean, surface_net_solar_radiation_min, surface_net_solar_radiation_max, surface_net_solar_radiation_mean, surface_net_thermal_radiation_min, surface_net_thermal_radiation_max, surface_net_thermal_radiation_mean, surface_pressure_min, surface_pressure_max, surface_pressure_mean, u_component_of_wind_10m_min, u_component_of_wind_10m_max, u_component_of_wind_10m_mean, v_component_of_wind_10m_min, v_component_of_wind_10m_max, v_component_of_wind_10m_mean, snow_depth_water_equivalent_min, snow_depth_water_equivalent_max, snow_depth_water_equivalent_mean, volumetric_soil_water_layer_1_min, volumetric_soil_water_layer_1_max, volumetric_soil_water_layer_1_mean, volumetric_soil_water_layer_2_min, volumetric_soil_water_layer_2_max, volumetric_soil_water_layer_2_mean, volumetric_soil_water_layer_3_min, volumetric_soil_water_layer_3_max, volumetric_soil_water_layer_3_mean, volumetric_soil_water_layer_4_min, volumetric_soil_water_layer_4_max, volumetric_soil_water_layer_4_mean] = loadCatchmentCaravan_camelscl(ID, path)
%loadCatchmentCaravan_camelscl loads timeseries (Q, P, PET, Tmin, Tmax, Tmean, 
%    dewpoint_Tmin, dewpoint_Tmax, dewpoint_Tmean, surface_net_solar_radiation_min, 
%    surface_net_solar_radiation_max, surface_net_solar_radiation_mean, 
%    surface_net_thermal_radiation_min, surface_net_thermal_radiation_max, 
%    surface_net_thermal_radiation_mean, surface_pressure_min, 
%    surface_pressure_max, surface_pressure_mean, u_component_of_wind_10m_min, 
%    u_component_of_wind_10m_max, u_component_of_wind_10m_mean, v_component_of_wind_10m_min,
%    v_component_of_wind_10m_max, v_component_of_wind_10m_mean, snow_depth_water_equivalent_min,
%    snow_depth_water_equivalent_max, snow_depth_water_equivalent_mean, 
%    volumetric_soil_water_layer_1_min, volumetric_soil_water_layer_1_max, 
%    volumetric_soil_water_layer_1_mean, volumetric_soil_water_layer_2_min, 
%    volumetric_soil_water_layer_2_max, volumetric_soil_water_layer_2_mean, 
%    volumetric_soil_water_layer_3_min, volumetric_soil_water_layer_3_max, 
%    volumetric_soil_water_layer_3_mean, volumetric_soil_water_layer_4_min, 
%    volumetric_soil_water_layer_4_max, volumetric_soil_water_layer_4_mean) 
%    for Caravan_camelscl format.
%
%   INPUT
%   ID: catchment ID in Caravan_camelscl
%   path: file path in local PC
%
%   OUTPUT
%   snow_depth_water_equivalent_mean    - Snow water equivalend daily mean       [mm] 
%   surface_net_solar_radiation_mean    - Shortwave radiation mean               [W/m2]
%   surface_net_thermal_radiation_mean  - Net thermal radiation at surface mean  [W/m2]
%   surface_pressure_mean               - Surface pressure mean                  [kPa]
%   Tmean                               - Air temperature mean                   [degC]
%   dewpoint_Tmean                      - Dew point temperature mean             [degC]
%   u_component_of_wind_10m_mean        - Eastward wind component daily mean     [m/s]
%   v_component_of_wind_10m_mean        - Northward wind component daily mean    [m/s]
%   volumetric_soil_water_layer_1_mean  - Soil water volume 0-7cm daily mean     [m3/m3]
%   volumetric_soil_water_layer_2_mean  - Soil water volume 7-28cm daily mean    [m3/m3]
%   volumetric_soil_water_layer_3_mean  - Soil water volume 28-100cm daily mean  [m3/m3]
%   volumetric_soil_water_layer_4_mean  - Soil water volume 100-289cm daily mean [m3/m3]
%
%   snow_depth_water_equivalent_min     - Snow water equivalend daily min        [mm] 
%   surface_net_solar_radiation_min     - Shortwave radiation min                [W/m2]
%   surface_net_thermal_radiation_min   - Net thermal radiation at surface min   [W/m2]
%   surface_pressure_min                - Surface pressure min                   [kPa]
%   Tmin                                - Air temperature min                    [degC]
%   dewpoint_Tmin                       - Dew point temperature min              [degC]
%   u_component_of_wind_10m_min         - Eastward wind component daily min      [m/s]
%   v_component_of_wind_10m_min         - Northward wind component daily min     [m/s]
%   volumetric_soil_water_layer_1_min   - Soil water volume 0-7cm daily min      [m3/m3]
%   volumetric_soil_water_layer_2_min   - Soil water volume 7-28cm daily min     [m3/m3]
%   volumetric_soil_water_layer_3_min   - Soil water volume 28-100cm daily min   [m3/m3]
%   volumetric_soil_water_layer_4_min   - Soil water volume 100-289cm daily min  [m3/m3]
%
%   snow_depth_water_equivalent_max     - Snow water equivalend daily max        [mm]
%   surface_net_solar_radiation_max     - Shortwave radiation max                [W/m2]
%   surface_net_thermal_radiation_max   - Net thermal radiation at surface max   [W/m2]
%   surface_pressure_max                - Surface pressure max                   [kPa]
%   Tmax                                - Air temperature max                    [degC]
%   dewpoint_Tmax                       - Dew point temperature max              [degC]
%   u_component_of_wind_10m_max         - Eastward wind component daily max      [m/s]
%   v_component_of_wind_10m_max         - Northward wind component daily max     [m/s]
%   volumetric_soil_water_layer_1_max   - Soil water volume 0-7cm daily max      [m3/m3]
%   volumetric_soil_water_layer_2_max   - Soil water volume 7-28cm daily max     [m3/m3]
%   volumetric_soil_water_layer_3_max   - Soil water volume 28-100cm daily max   [m3/m3]
%   volumetric_soil_water_layer_4_max   - Soil water volume 100-289cm daily max  [m3/m3]
%   
%   P                                   - Precipitation daily sum                [mm/day]
%   PET                                 - Potential evaporation daily sum        [mm/day]
%   Q                                   - Daily streamflow                       [mm/day]
%
%   Copyright (C) 2023
%   This code is referenced from SebastianGnann's github
%   see <https://github.com/SebastianGnann/CAMELS_Matlab>

% check input parameters
if nargin < 2
    error('Not enought input arguments.')
end

file_ID = strcat(path, 'camelscl_',num2str(ID),'.csv');

% 40 columns in each csv files
data = readtable(file_ID);

date = datenum(data.date);
snow_depth_water_equivalent_mean_temp = data.snow_depth_water_equivalent_mean;
surface_net_solar_radiation_mean_temp = data.surface_net_solar_radiation_mean;
surface_net_thermal_radiation_mean_temp = data.surface_net_thermal_radiation_mean;
surface_pressure_mean_temp = data.surface_pressure_mean;
Tmean_temp = data.temperature_2m_mean; 
dewpoint_Tmean_temp = data.dewpoint_temperature_2m_mean;
u_component_of_wind_10m_mean_temp = data.u_component_of_wind_10m_mean;
v_component_of_wind_10m_mean_temp = data.v_component_of_wind_10m_mean;
volumetric_soil_water_layer_1_mean_temp = data.volumetric_soil_water_layer_1_mean;
volumetric_soil_water_layer_2_mean_temp = data.volumetric_soil_water_layer_2_mean;
volumetric_soil_water_layer_3_mean_temp = data.volumetric_soil_water_layer_3_mean;
volumetric_soil_water_layer_4_mean_temp = data.volumetric_soil_water_layer_4_mean;

snow_depth_water_equivalent_min_temp = data.snow_depth_water_equivalent_min;
surface_net_solar_radiation_min_temp = data.surface_net_solar_radiation_min;
surface_net_thermal_radiation_min_temp = data.surface_net_thermal_radiation_min;
surface_pressure_min_temp = data.surface_pressure_min;
Tmin_temp = data.temperature_2m_min; 
dewpoint_Tmin_temp = data.dewpoint_temperature_2m_min;
u_component_of_wind_10m_min_temp = data.u_component_of_wind_10m_min;
v_component_of_wind_10m_min_temp = data.v_component_of_wind_10m_min;
volumetric_soil_water_layer_1_min_temp = data.volumetric_soil_water_layer_1_min;
volumetric_soil_water_layer_2_min_temp = data.volumetric_soil_water_layer_2_min;
volumetric_soil_water_layer_3_min_temp = data.volumetric_soil_water_layer_3_min;
volumetric_soil_water_layer_4_min_temp = data.volumetric_soil_water_layer_4_min;

snow_depth_water_equivalent_max_temp = data.snow_depth_water_equivalent_max;
surface_net_solar_radiation_max_temp = data.surface_net_solar_radiation_max;
surface_net_thermal_radiation_max_temp = data.surface_net_thermal_radiation_max;
surface_pressure_max_temp = data.surface_pressure_max;
Tmax_temp = data.temperature_2m_max; 
dewpoint_Tmax_temp = data.dewpoint_temperature_2m_max;
u_component_of_wind_10m_max_temp = data.u_component_of_wind_10m_max;
v_component_of_wind_10m_max_temp = data.v_component_of_wind_10m_max;
volumetric_soil_water_layer_1_max_temp = data.volumetric_soil_water_layer_1_max;
volumetric_soil_water_layer_2_max_temp = data.volumetric_soil_water_layer_2_max;
volumetric_soil_water_layer_3_max_temp = data.volumetric_soil_water_layer_3_max;
volumetric_soil_water_layer_4_max_temp = data.volumetric_soil_water_layer_4_max;

P_temp = data.total_precipitation_sum;
PET_temp = data.potential_evaporation_sum;
Q_temp = data.streamflow;

% Construct final dataframe
snow_depth_water_equivalent_mean = [date, snow_depth_water_equivalent_mean_temp];
surface_net_solar_radiation_mean = [date, surface_net_solar_radiation_mean_temp];
surface_net_thermal_radiation_mean = [date, surface_net_thermal_radiation_mean_temp];
surface_pressure_mean = [date, surface_pressure_mean_temp];
Tmean = [date, Tmean_temp];
dewpoint_Tmean = [date, dewpoint_Tmean_temp];
u_component_of_wind_10m_mean = [date, u_component_of_wind_10m_mean_temp];
v_component_of_wind_10m_mean = [date, v_component_of_wind_10m_mean_temp];
volumetric_soil_water_layer_1_mean = [date, volumetric_soil_water_layer_1_mean_temp];
volumetric_soil_water_layer_2_mean = [date, volumetric_soil_water_layer_2_mean_temp];
volumetric_soil_water_layer_3_mean = [date, volumetric_soil_water_layer_3_mean_temp];
volumetric_soil_water_layer_4_mean = [date, volumetric_soil_water_layer_4_mean_temp];

snow_depth_water_equivalent_min = [date, snow_depth_water_equivalent_min_temp];
surface_net_solar_radiation_min = [date, surface_net_solar_radiation_min_temp];
surface_net_thermal_radiation_min = [date, surface_net_thermal_radiation_min_temp];
surface_pressure_min = [date, surface_pressure_min_temp];
Tmin = [date, Tmin_temp];
dewpoint_Tmin = [date, dewpoint_Tmin_temp];
u_component_of_wind_10m_min = [date, u_component_of_wind_10m_min_temp];
v_component_of_wind_10m_min = [date, v_component_of_wind_10m_min_temp];
volumetric_soil_water_layer_1_min = [date, volumetric_soil_water_layer_1_min_temp];
volumetric_soil_water_layer_2_min = [date, volumetric_soil_water_layer_2_min_temp];
volumetric_soil_water_layer_3_min = [date, volumetric_soil_water_layer_3_min_temp];
volumetric_soil_water_layer_4_min = [date, volumetric_soil_water_layer_4_min_temp];

snow_depth_water_equivalent_max = [date, snow_depth_water_equivalent_max_temp];
surface_net_solar_radiation_max = [date, surface_net_solar_radiation_max_temp];
surface_net_thermal_radiation_max = [date, surface_net_thermal_radiation_max_temp];
surface_pressure_max = [date, surface_pressure_max_temp];
Tmax = [date, Tmax_temp];
dewpoint_Tmax = [date, dewpoint_Tmax_temp];
u_component_of_wind_10m_max = [date, u_component_of_wind_10m_max_temp];
v_component_of_wind_10m_max = [date, v_component_of_wind_10m_max_temp];
volumetric_soil_water_layer_1_max = [date, volumetric_soil_water_layer_1_max_temp];
volumetric_soil_water_layer_2_max = [date, volumetric_soil_water_layer_2_max_temp];
volumetric_soil_water_layer_3_max = [date, volumetric_soil_water_layer_3_max_temp];
volumetric_soil_water_layer_4_max = [date, volumetric_soil_water_layer_4_max_temp];

P = [date, P_temp];
PET = [date, PET_temp];
Q = [date, Q_temp];

end