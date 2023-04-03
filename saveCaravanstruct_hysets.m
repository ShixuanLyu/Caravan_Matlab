function [Caravan_hysets_data] = saveCaravanstruct_hysets(save_struct)
%saveCaravanstruct_hysets creates struct file with Caravan_hysets data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Timeseries are loaded for the period same with the data in Caravan
%   - Users local PC paths (which you save your Caravan_hysets data)
%   - Data can be found at: https://doi.org/10.5281/zenodo.7540792
%   - Caravan project can be found at: https://github.com/kratzert/Caravan/
%
%  INPUT
%  save_struct: whether to save struct file or not
% 
%  OUTPUT
%  Caravan_hysets_data: struct file with Caravan_hysets data
%
%  References
%  Kratzert, F., Nearing, G., Addor, N. et al. Caravan - A global community
%  dataset for large-sample hydrology. Sci Data 10, 61 (2023). 
%  
%  https://github.com/kratzert/Caravan/
%  https://github.com/SebastianGnann/CAMELS_Matlab
%  
%  Copyright (C) 2023

if nargin < 1
    save_struct = false;
end

%% Specify paths
% Caravan datasets are consist by five folders. Timeseries data are stored
% with csv files and netcdf files. Here we only use the csv files. 
% The Caravan datasets has been divided into 7 sub-datasets. This code only
% conduct data processing with Caravan_hysets. 
% We have to be above the Caravan directory and Caravan hysets data should be
% stored in a folder named hysets. The following folders are required:
% Caravan/timeseries/csv/hysets/hysets_*.csv
% (482 files; contain hydrometeorological time series)
% Caravan/attributes/hysets/attributes_*_hysets.csv
% (3 files; contain catchment attributes)

% This is my local path, don't forget to change your path here. 
path_catchment_attributes = "Caravan/attributes/hysets/";
path_time_series = "Caravan/timeseries/csv/hysets/";

if ~(exist(path_catchment_attributes) == 7)
    error('Cannot find local path. You can download Caravan from https://doi.org/10.5281/zenodo.7540792.')
elseif ~(exist(path_time_series) == 7)
    error('Cannot find local path. You can download Caravan from https://doi.org/10.5281/zenodo.7540792.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several csv
% files.

% attributes_caravan_hysets
% gauge_id      p_mean      pet_mean        aridity     frac_snow       moisture_index      seasonality       high_prec_freq      high_prec_dur    low_prec_freq      low_prec_dur
% [attributes_caravan_hysets_data, attributes_caravan_hysets_data_str] =...
%     xlsread(strcat(path_catchment_attributes, 'attributes_caravan_hysets'))
attributes_caravan_hysets_data = readtable(strcat(path_catchment_attributes, 'attributes_caravan_hysets.csv'));

% attributes_hydroatlas_hysets
% gauge_id      sgr_dk_sav      glc_pc_s06      glc_pc_s07      nli_ix_sav      glc_pc_s04        glc_pc_s05      glc_pc_s02      rev_mc_usu        glc_pc_s03        glc_pc_s01      pet_mm_syr      dor_pc_pva      glc_pc_s08        glc_pc_s09      swc_pc_s09      ele_mt_smx        tbi_cl_smj        swc_pc_s01      swc_pc_s02      swc_pc_s03        swc_pc_s04        swc_pc_s05      swc_pc_s06      swc_pc_s07        swc_pc_s08        crp_pc_sse      glc_pc_s22      glc_pc_s20        glc_pc_s21        wet_pc_sg1      wet_pc_sg2      pac_pc_sse        swc_pc_s10        swc_pc_s11      swc_pc_s12      clz_cl_smj        gwt_cm_sav        glc_pc_s17      glc_pc_s18      hft_ix_s93        glc_pc_s15        ire_pc_sse      glc_pc_s16      glc_pc_s13        prm_pc_sse        glc_pc_s14      glc_pc_s11      glc_pc_s12        glc_pc_s10        kar_pc_sse      slp_dg_sav      glc_pc_s19        tmp_dc_s07        tmp_dc_s08      tmp_dc_s05      tmp_dc_s06        tmp_dc_s09        for_pc_sse      aet_mm_s06      aet_mm_s05        aet_mm_s08        aet_mm_s07      aet_mm_s09      tmp_dc_s10        tmp_dc_s11        aet_mm_s02      aet_mm_s01      tmp_dc_s12        aet_mm_s04        aet_mm_s03      lit_cl_smj      tmp_dc_s03      tmp_dc_s04        tmp_dc_s01      tmp_dc_s02      cls_cl_smj        pre_mm_syr        pnv_pc_s01      pnv_pc_s04      pnv_pc_s05        pnv_pc_s02        rdd_mk_sav      ele_mt_smn      pnv_pc_s03        pnv_pc_s08        pnv_pc_s09      pnv_pc_s06      pnv_pc_s07        wet_cl_smj        snw_pc_syr      pnv_pc_s11      pnv_pc_s12        pnv_pc_s10        pnv_pc_s15      pnv_pc_s13      pnv_pc_s14        cmi_ix_syr        pet_mm_s11      pet_mm_s12      pet_mm_s10        tmp_dc_smn        wet_pc_s08      wet_pc_s09      slt_pc_sav        wet_pc_s02        wet_pc_s03      wet_pc_s01      hdi_ix_sav        wet_pc_s06        wet_pc_s07      wet_pc_s04      wet_pc_s05        fec_cl_smj        glc_cl_smj      swc_pc_syr      hft_ix_s09        soc_th_sav        gdp_ud_sav      dis_m3_pyr      gdp_ud_ssu        tmp_dc_smx        cly_pc_sav      pet_mm_s02      pet_mm_s03        pet_mm_s01        riv_tc_usu      snw_pc_smx      ppd_pk_sav        pet_mm_s08        aet_mm_s11      pet_mm_s09      aet_mm_s10        pet_mm_s06        pet_mm_s07      aet_mm_s12      pet_mm_s04        pet_mm_s05        inu_pc_slt      ero_kj_sav      aet_mm_syr        cmi_ix_s10        cmi_ix_s11      cmi_ix_s12      ari_ix_sav        tmp_dc_syr        tec_cl_smj      ria_ha_usu      lkv_mc_usu        fmh_cl_smj        inu_pc_smn      pnv_cl_smj      pre_mm_s08        pre_mm_s09        run_mm_syr      pre_mm_s06      pre_mm_s07        pre_mm_s04        pre_mm_s05      snd_pc_sav      pre_mm_s02        pre_mm_s03        ele_mt_sav      pre_mm_s01      urb_pc_sse        lka_pc_sse        pre_mm_s10      dis_m3_pmx      snw_pc_s01        snw_pc_s02        snw_pc_s03      snw_pc_s04      snw_pc_s05        snw_pc_s06        gla_pc_sse      snw_pc_s07      snw_pc_s08        snw_pc_s09        dis_m3_pmn      inu_pc_smx      pre_mm_s11        pre_mm_s12        cmi_ix_s07      cmi_ix_s08      cmi_ix_s05        cmi_ix_s06        cmi_ix_s09      snw_pc_s10      snw_pc_s11      snw_pc_s12        cmi_ix_s03      cmi_ix_s04      cmi_ix_s01        cmi_ix_s02        pst_pc_sse      area_fraction_used_for_aggregation        pop_ct_usu
% [attributes_hydroatlas_hysets_data,attributes_hydroatlas_hysets_data_str] = ...
%     xlsread(strcat(path_catchment_attributes, 'attributes_hydroatlas_hysets'))
attributes_hydroatlas_hysets_data = readtable(strcat(path_catchment_attributes, 'attributes_hydroatlas_hysets.csv'));

% attributes_other_hysets
% gauge_id      gauge_name      country     gauge_lat       gauge_lon       area
% [attributes_other_hysets_data, attributes_other_hysets_data_str] = ...
%      xlsread(strcat(path_catchment_attributes, 'attributes_other_hysets'))
attributes_other_hysets_data = readtable(strcat(path_catchment_attributes, 'attributes_other_hysets.csv'));

%% We now add the catchment attributes and metadata to the struct file.

% add attributes_caravan_hysets
Caravan_hysets_data.gauge_index = attributes_caravan_hysets_data.gauge_index;
Caravan_hysets_data.gauge_id = attributes_caravan_hysets_data.gauge_id;
Caravan_hysets_data.p_mean = attributes_caravan_hysets_data.p_mean;
Caravan_hysets_data.pet_mean = attributes_caravan_hysets_data.pet_mean;
Caravan_hysets_data.aridity = attributes_caravan_hysets_data.aridity;
Caravan_hysets_data.frac_snow = attributes_caravan_hysets_data.frac_snow;
Caravan_hysets_data.moisture_index = attributes_caravan_hysets_data.moisture_index;
Caravan_hysets_data.seasonality = attributes_caravan_hysets_data.seasonality;
Caravan_hysets_data.high_prec_freq = attributes_caravan_hysets_data.high_prec_freq;
Caravan_hysets_data.high_prec_dur = attributes_caravan_hysets_data.high_prec_dur;
Caravan_hysets_data.low_prec_freq = attributes_caravan_hysets_data.low_prec_freq;
Caravan_hysets_data.low_prec_dur = attributes_caravan_hysets_data.low_prec_dur;

% add attributes_hydroatlas_hysets
Caravan_hysets_data.sgr_dk_sav = attributes_hydroatlas_hysets_data.sgr_dk_sav;
Caravan_hysets_data.glc_pc_s06 = attributes_hydroatlas_hysets_data.glc_pc_s06;
Caravan_hysets_data.glc_pc_s07 = attributes_hydroatlas_hysets_data.glc_pc_s07;
Caravan_hysets_data.nli_ix_sav = attributes_hydroatlas_hysets_data.nli_ix_sav;
Caravan_hysets_data.glc_pc_s04 = attributes_hydroatlas_hysets_data.glc_pc_s04;
Caravan_hysets_data.glc_pc_s05 = attributes_hydroatlas_hysets_data.glc_pc_s05;
Caravan_hysets_data.glc_pc_s02 = attributes_hydroatlas_hysets_data.glc_pc_s02;
Caravan_hysets_data.rev_mc_usu = attributes_hydroatlas_hysets_data.rev_mc_usu;
Caravan_hysets_data.glc_pc_s03 = attributes_hydroatlas_hysets_data.glc_pc_s03;
Caravan_hysets_data.glc_pc_s01 = attributes_hydroatlas_hysets_data.glc_pc_s01;
Caravan_hysets_data.pet_mm_syr = attributes_hydroatlas_hysets_data.pet_mm_syr;
Caravan_hysets_data.dor_pc_pva = attributes_hydroatlas_hysets_data.dor_pc_pva;
Caravan_hysets_data.glc_pc_s08 = attributes_hydroatlas_hysets_data.glc_pc_s08;
Caravan_hysets_data.glc_pc_s09 = attributes_hydroatlas_hysets_data.glc_pc_s09;
Caravan_hysets_data.swc_pc_s09 = attributes_hydroatlas_hysets_data.swc_pc_s09;
Caravan_hysets_data.ele_mt_smx = attributes_hydroatlas_hysets_data.ele_mt_smx;
Caravan_hysets_data.tbi_cl_smj = attributes_hydroatlas_hysets_data.tbi_cl_smj;
Caravan_hysets_data.swc_pc_s01 = attributes_hydroatlas_hysets_data.swc_pc_s01;
Caravan_hysets_data.swc_pc_s02 = attributes_hydroatlas_hysets_data.swc_pc_s02;
Caravan_hysets_data.swc_pc_s03 = attributes_hydroatlas_hysets_data.swc_pc_s03;
Caravan_hysets_data.swc_pc_s04 = attributes_hydroatlas_hysets_data.swc_pc_s04;
Caravan_hysets_data.swc_pc_s05 = attributes_hydroatlas_hysets_data.swc_pc_s05;
Caravan_hysets_data.swc_pc_s06 = attributes_hydroatlas_hysets_data.swc_pc_s06;
Caravan_hysets_data.swc_pc_s07 = attributes_hydroatlas_hysets_data.swc_pc_s07;
Caravan_hysets_data.swc_pc_s08 = attributes_hydroatlas_hysets_data.swc_pc_s08;
Caravan_hysets_data.crp_pc_sse = attributes_hydroatlas_hysets_data.crp_pc_sse;
Caravan_hysets_data.glc_pc_s22 = attributes_hydroatlas_hysets_data.glc_pc_s22;
Caravan_hysets_data.glc_pc_s20 = attributes_hydroatlas_hysets_data.glc_pc_s20;
Caravan_hysets_data.glc_pc_s21 = attributes_hydroatlas_hysets_data.glc_pc_s21;
Caravan_hysets_data.wet_pc_sg1 = attributes_hydroatlas_hysets_data.wet_pc_sg1;
Caravan_hysets_data.wet_pc_sg2 = attributes_hydroatlas_hysets_data.wet_pc_sg2;
Caravan_hysets_data.pac_pc_sse = attributes_hydroatlas_hysets_data.pac_pc_sse;
Caravan_hysets_data.swc_pc_s10 = attributes_hydroatlas_hysets_data.swc_pc_s10;
Caravan_hysets_data.swc_pc_s11 = attributes_hydroatlas_hysets_data.swc_pc_s11;
Caravan_hysets_data.swc_pc_s12 = attributes_hydroatlas_hysets_data.swc_pc_s12;
Caravan_hysets_data.clz_cl_smj = attributes_hydroatlas_hysets_data.clz_cl_smj;
Caravan_hysets_data.gwt_cm_sav = attributes_hydroatlas_hysets_data.gwt_cm_sav;
Caravan_hysets_data.glc_pc_s17 = attributes_hydroatlas_hysets_data.glc_pc_s17;
Caravan_hysets_data.glc_pc_s18 = attributes_hydroatlas_hysets_data.glc_pc_s18;
Caravan_hysets_data.hft_ix_s93 = attributes_hydroatlas_hysets_data.hft_ix_s93;
Caravan_hysets_data.glc_pc_s15 = attributes_hydroatlas_hysets_data.glc_pc_s15;
Caravan_hysets_data.ire_pc_sse = attributes_hydroatlas_hysets_data.ire_pc_sse;
Caravan_hysets_data.glc_pc_s16 = attributes_hydroatlas_hysets_data.glc_pc_s16;
Caravan_hysets_data.glc_pc_s13 = attributes_hydroatlas_hysets_data.glc_pc_s13;
Caravan_hysets_data.prm_pc_sse = attributes_hydroatlas_hysets_data.prm_pc_sse;
Caravan_hysets_data.glc_pc_s14 = attributes_hydroatlas_hysets_data.glc_pc_s14;
Caravan_hysets_data.glc_pc_s11 = attributes_hydroatlas_hysets_data.glc_pc_s11;
Caravan_hysets_data.glc_pc_s12 = attributes_hydroatlas_hysets_data.glc_pc_s12;
Caravan_hysets_data.glc_pc_s10 = attributes_hydroatlas_hysets_data.glc_pc_s10;
Caravan_hysets_data.kar_pc_sse = attributes_hydroatlas_hysets_data.kar_pc_sse;
Caravan_hysets_data.slp_dg_sav = attributes_hydroatlas_hysets_data.slp_dg_sav;
Caravan_hysets_data.glc_pc_s19 = attributes_hydroatlas_hysets_data.glc_pc_s19;
Caravan_hysets_data.tmp_dc_s07 = attributes_hydroatlas_hysets_data.tmp_dc_s07;
Caravan_hysets_data.tmp_dc_s08 = attributes_hydroatlas_hysets_data.tmp_dc_s08;
Caravan_hysets_data.tmp_dc_s05 = attributes_hydroatlas_hysets_data.tmp_dc_s05;
Caravan_hysets_data.tmp_dc_s06 = attributes_hydroatlas_hysets_data.tmp_dc_s06;
Caravan_hysets_data.tmp_dc_s09 = attributes_hydroatlas_hysets_data.tmp_dc_s09;
Caravan_hysets_data.for_pc_sse = attributes_hydroatlas_hysets_data.for_pc_sse;
Caravan_hysets_data.aet_mm_s06 = attributes_hydroatlas_hysets_data.aet_mm_s06;
Caravan_hysets_data.aet_mm_s05 = attributes_hydroatlas_hysets_data.aet_mm_s05;
Caravan_hysets_data.aet_mm_s08 = attributes_hydroatlas_hysets_data.aet_mm_s08;
Caravan_hysets_data.aet_mm_s07 = attributes_hydroatlas_hysets_data.aet_mm_s07;
Caravan_hysets_data.aet_mm_s09 = attributes_hydroatlas_hysets_data.aet_mm_s09;
Caravan_hysets_data.tmp_dc_s10 = attributes_hydroatlas_hysets_data.tmp_dc_s10;
Caravan_hysets_data.tmp_dc_s11 = attributes_hydroatlas_hysets_data.tmp_dc_s11;
Caravan_hysets_data.aet_mm_s02 = attributes_hydroatlas_hysets_data.aet_mm_s02;
Caravan_hysets_data.aet_mm_s01 = attributes_hydroatlas_hysets_data.aet_mm_s01;
Caravan_hysets_data.tmp_dc_s12 = attributes_hydroatlas_hysets_data.tmp_dc_s12;
Caravan_hysets_data.aet_mm_s04 = attributes_hydroatlas_hysets_data.aet_mm_s04;
Caravan_hysets_data.aet_mm_s03 = attributes_hydroatlas_hysets_data.aet_mm_s03;
Caravan_hysets_data.lit_cl_smj = attributes_hydroatlas_hysets_data.lit_cl_smj;
Caravan_hysets_data.tmp_dc_s03 = attributes_hydroatlas_hysets_data.tmp_dc_s03;
Caravan_hysets_data.tmp_dc_s04 = attributes_hydroatlas_hysets_data.tmp_dc_s04;
Caravan_hysets_data.tmp_dc_s01 = attributes_hydroatlas_hysets_data.tmp_dc_s01;
Caravan_hysets_data.tmp_dc_s02 = attributes_hydroatlas_hysets_data.tmp_dc_s02;
Caravan_hysets_data.cls_cl_smj = attributes_hydroatlas_hysets_data.cls_cl_smj;
Caravan_hysets_data.pre_mm_syr = attributes_hydroatlas_hysets_data.pre_mm_syr;
Caravan_hysets_data.pnv_pc_s01 = attributes_hydroatlas_hysets_data.pnv_pc_s01;
Caravan_hysets_data.pnv_pc_s04 = attributes_hydroatlas_hysets_data.pnv_pc_s04;
Caravan_hysets_data.pnv_pc_s05 = attributes_hydroatlas_hysets_data.pnv_pc_s05;
Caravan_hysets_data.pnv_pc_s02 = attributes_hydroatlas_hysets_data.pnv_pc_s02;
Caravan_hysets_data.rdd_mk_sav = attributes_hydroatlas_hysets_data.rdd_mk_sav;
Caravan_hysets_data.ele_mt_smn = attributes_hydroatlas_hysets_data.ele_mt_smn;
Caravan_hysets_data.pnv_pc_s03 = attributes_hydroatlas_hysets_data.pnv_pc_s03;
Caravan_hysets_data.pnv_pc_s08 = attributes_hydroatlas_hysets_data.pnv_pc_s08;
Caravan_hysets_data.pnv_pc_s09 = attributes_hydroatlas_hysets_data.pnv_pc_s09;
Caravan_hysets_data.pnv_pc_s06 = attributes_hydroatlas_hysets_data.pnv_pc_s06;
Caravan_hysets_data.pnv_pc_s07 = attributes_hydroatlas_hysets_data.pnv_pc_s07;
Caravan_hysets_data.wet_cl_smj = attributes_hydroatlas_hysets_data.wet_cl_smj;
Caravan_hysets_data.snw_pc_syr = attributes_hydroatlas_hysets_data.snw_pc_syr;
Caravan_hysets_data.pnv_pc_s11 = attributes_hydroatlas_hysets_data.pnv_pc_s11;
Caravan_hysets_data.pnv_pc_s12 = attributes_hydroatlas_hysets_data.pnv_pc_s12;
Caravan_hysets_data.pnv_pc_s10 = attributes_hydroatlas_hysets_data.pnv_pc_s10;
Caravan_hysets_data.pnv_pc_s15 = attributes_hydroatlas_hysets_data.pnv_pc_s15;
Caravan_hysets_data.pnv_pc_s13 = attributes_hydroatlas_hysets_data.pnv_pc_s13;
Caravan_hysets_data.pnv_pc_s14 = attributes_hydroatlas_hysets_data.pnv_pc_s14;
Caravan_hysets_data.cmi_ix_syr = attributes_hydroatlas_hysets_data.cmi_ix_syr;
Caravan_hysets_data.pet_mm_s11 = attributes_hydroatlas_hysets_data.pet_mm_s11;
Caravan_hysets_data.pet_mm_s12 = attributes_hydroatlas_hysets_data.pet_mm_s10;
Caravan_hysets_data.pet_mm_s10 = attributes_hydroatlas_hysets_data.pet_mm_s10;
Caravan_hysets_data.tmp_dc_smn = attributes_hydroatlas_hysets_data.tmp_dc_smn;
Caravan_hysets_data.wet_pc_s08 = attributes_hydroatlas_hysets_data.wet_pc_s08;
Caravan_hysets_data.wet_pc_s09 = attributes_hydroatlas_hysets_data.wet_pc_s09;
Caravan_hysets_data.slt_pc_sav = attributes_hydroatlas_hysets_data.slt_pc_sav;
Caravan_hysets_data.wet_pc_s02 = attributes_hydroatlas_hysets_data.wet_pc_s02;
Caravan_hysets_data.wet_pc_s03 = attributes_hydroatlas_hysets_data.wet_pc_s03;
Caravan_hysets_data.wet_pc_s01 = attributes_hydroatlas_hysets_data.wet_pc_s01;
Caravan_hysets_data.hdi_ix_sav = attributes_hydroatlas_hysets_data.hdi_ix_sav;
Caravan_hysets_data.wet_pc_s06 = attributes_hydroatlas_hysets_data.wet_pc_s06;
Caravan_hysets_data.wet_pc_s07 = attributes_hydroatlas_hysets_data.wet_pc_s07;
Caravan_hysets_data.wet_pc_s04 = attributes_hydroatlas_hysets_data.wet_pc_s04;
Caravan_hysets_data.wet_pc_s05 = attributes_hydroatlas_hysets_data.wet_pc_s05;
Caravan_hysets_data.fec_cl_smj = attributes_hydroatlas_hysets_data.fec_cl_smj;
Caravan_hysets_data.glc_cl_smj = attributes_hydroatlas_hysets_data.glc_cl_smj;
Caravan_hysets_data.swc_pc_syr = attributes_hydroatlas_hysets_data.swc_pc_syr;
Caravan_hysets_data.hft_ix_s09 = attributes_hydroatlas_hysets_data.hft_ix_s09;
Caravan_hysets_data.soc_th_sav = attributes_hydroatlas_hysets_data.soc_th_sav;
Caravan_hysets_data.gdp_ud_sav = attributes_hydroatlas_hysets_data.gdp_ud_sav;
Caravan_hysets_data.dis_m3_pyr = attributes_hydroatlas_hysets_data.dis_m3_pyr;
Caravan_hysets_data.gdp_ud_ssu = attributes_hydroatlas_hysets_data.gdp_ud_ssu;
Caravan_hysets_data.tmp_dc_smx = attributes_hydroatlas_hysets_data.tmp_dc_smx;
Caravan_hysets_data.cly_pc_sav = attributes_hydroatlas_hysets_data.cly_pc_sav;
Caravan_hysets_data.pet_mm_s02 = attributes_hydroatlas_hysets_data.pet_mm_s02;
Caravan_hysets_data.pet_mm_s03 = attributes_hydroatlas_hysets_data.pet_mm_s03;
Caravan_hysets_data.pet_mm_s01 = attributes_hydroatlas_hysets_data.pet_mm_s01;
Caravan_hysets_data.riv_tc_usu = attributes_hydroatlas_hysets_data.riv_tc_usu;
Caravan_hysets_data.snw_pc_smx = attributes_hydroatlas_hysets_data.snw_pc_smx;
Caravan_hysets_data.ppd_pk_sav = attributes_hydroatlas_hysets_data.ppd_pk_sav;
Caravan_hysets_data.pet_mm_s08 = attributes_hydroatlas_hysets_data.pet_mm_s08;
Caravan_hysets_data.aet_mm_s11 = attributes_hydroatlas_hysets_data.aet_mm_s11;
Caravan_hysets_data.pet_mm_s09 = attributes_hydroatlas_hysets_data.pet_mm_s09;
Caravan_hysets_data.aet_mm_s10 = attributes_hydroatlas_hysets_data.aet_mm_s10;
Caravan_hysets_data.pet_mm_s06 = attributes_hydroatlas_hysets_data.pet_mm_s06;
Caravan_hysets_data.pet_mm_s07 = attributes_hydroatlas_hysets_data.pet_mm_s07;
Caravan_hysets_data.aet_mm_s12 = attributes_hydroatlas_hysets_data.aet_mm_s12;
Caravan_hysets_data.pet_mm_s04 = attributes_hydroatlas_hysets_data.pet_mm_s04;
Caravan_hysets_data.pet_mm_s05 = attributes_hydroatlas_hysets_data.pet_mm_s05;
Caravan_hysets_data.inu_pc_slt = attributes_hydroatlas_hysets_data.inu_pc_slt;
Caravan_hysets_data.ero_kh_sav = attributes_hydroatlas_hysets_data.ero_kh_sav;
Caravan_hysets_data.aet_mm_syr = attributes_hydroatlas_hysets_data.aet_mm_syr;
Caravan_hysets_data.cmi_ix_s10 = attributes_hydroatlas_hysets_data.cmi_ix_s10;
Caravan_hysets_data.cmi_ix_s11 = attributes_hydroatlas_hysets_data.cmi_ix_s11;
Caravan_hysets_data.cmi_ix_s12 = attributes_hydroatlas_hysets_data.cmi_ix_s12;
Caravan_hysets_data.ari_ix_sav = attributes_hydroatlas_hysets_data.ari_ix_sav;
Caravan_hysets_data.tmp_dc_syr = attributes_hydroatlas_hysets_data.tmp_dc_syr;
Caravan_hysets_data.tec_cl_smj = attributes_hydroatlas_hysets_data.tec_cl_smj;
Caravan_hysets_data.ria_ha_usu = attributes_hydroatlas_hysets_data.ria_ha_usu;
Caravan_hysets_data.lkv_mc_usu = attributes_hydroatlas_hysets_data.lkv_mc_usu;
Caravan_hysets_data.fmh_cl_smj = attributes_hydroatlas_hysets_data.fmh_cl_smj;
Caravan_hysets_data.inu_pc_smn = attributes_hydroatlas_hysets_data.inu_pc_smn;
Caravan_hysets_data.pnv_cl_smj = attributes_hydroatlas_hysets_data.pnv_cl_smj;
Caravan_hysets_data.pre_mm_s08 = attributes_hydroatlas_hysets_data.pre_mm_s08;
Caravan_hysets_data.pre_mm_s09 = attributes_hydroatlas_hysets_data.pre_mm_s09;
Caravan_hysets_data.run_mm_syr = attributes_hydroatlas_hysets_data.run_mm_syr;
Caravan_hysets_data.pre_mm_s06 = attributes_hydroatlas_hysets_data.pre_mm_s06;
Caravan_hysets_data.pre_mm_s07 = attributes_hydroatlas_hysets_data.pre_mm_s07;
Caravan_hysets_data.pre_mm_s04 = attributes_hydroatlas_hysets_data.pre_mm_s04;
Caravan_hysets_data.pre_mm_s05 = attributes_hydroatlas_hysets_data.pre_mm_s05;
Caravan_hysets_data.snd_pc_sav = attributes_hydroatlas_hysets_data.snd_pc_sav;
Caravan_hysets_data.pre_mm_s02 = attributes_hydroatlas_hysets_data.pre_mm_s02;
Caravan_hysets_data.pre_mm_s03 = attributes_hydroatlas_hysets_data.pre_mm_s03;
Caravan_hysets_data.ele_mt_sav = attributes_hydroatlas_hysets_data.ele_mt_sav;
Caravan_hysets_data.pre_mm_s01 = attributes_hydroatlas_hysets_data.pre_mm_s01;
Caravan_hysets_data.urb_pc_sse = attributes_hydroatlas_hysets_data.urb_pc_sse;
Caravan_hysets_data.lka_pc_sse = attributes_hydroatlas_hysets_data.lka_pc_sse;
Caravan_hysets_data.pre_mm_s10 = attributes_hydroatlas_hysets_data.pre_mm_s10;
Caravan_hysets_data.dis_m3_pmx = attributes_hydroatlas_hysets_data.dis_m3_pmx;
Caravan_hysets_data.snw_pc_s01 = attributes_hydroatlas_hysets_data.snw_pc_s01;
Caravan_hysets_data.snw_pc_s02 = attributes_hydroatlas_hysets_data.snw_pc_s02;
Caravan_hysets_data.snw_pc_s03 = attributes_hydroatlas_hysets_data.snw_pc_s03;
Caravan_hysets_data.snw_pc_s04 = attributes_hydroatlas_hysets_data.snw_pc_s04;
Caravan_hysets_data.snw_pc_s05 = attributes_hydroatlas_hysets_data.snw_pc_s05;
Caravan_hysets_data.snw_pc_s06 = attributes_hydroatlas_hysets_data.snw_pc_s06;
Caravan_hysets_data.gla_pc_sse = attributes_hydroatlas_hysets_data.gla_pc_sse;
Caravan_hysets_data.snw_pc_s07 = attributes_hydroatlas_hysets_data.snw_pc_s07;
Caravan_hysets_data.snw_pc_s08 = attributes_hydroatlas_hysets_data.snw_pc_s08;
Caravan_hysets_data.snw_pc_s09 = attributes_hydroatlas_hysets_data.snw_pc_s09;
Caravan_hysets_data.dis_m3_pmn = attributes_hydroatlas_hysets_data.dis_m3_pmn;
Caravan_hysets_data.inu_pc_smx = attributes_hydroatlas_hysets_data.inu_pc_smx;
Caravan_hysets_data.pre_mm_s11 = attributes_hydroatlas_hysets_data.pre_mm_s11;
Caravan_hysets_data.pre_mm_s12 = attributes_hydroatlas_hysets_data.pre_mm_s12;
Caravan_hysets_data.cmi_ix_s07 = attributes_hydroatlas_hysets_data.cmi_ix_s07;
Caravan_hysets_data.cmi_ix_s08 = attributes_hydroatlas_hysets_data.cmi_ix_s08;
Caravan_hysets_data.cmi_ix_s05 = attributes_hydroatlas_hysets_data.cmi_ix_s05;
Caravan_hysets_data.cmi_ix_s06 = attributes_hydroatlas_hysets_data.cmi_ix_s06;
Caravan_hysets_data.cmi_ix_s09 = attributes_hydroatlas_hysets_data.cmi_ix_s09;
Caravan_hysets_data.snw_pc_s10 = attributes_hydroatlas_hysets_data.snw_pc_s10;
Caravan_hysets_data.snw_pc_s11 = attributes_hydroatlas_hysets_data.snw_pc_s11;
Caravan_hysets_data.snw_pc_s12 = attributes_hydroatlas_hysets_data.snw_pc_s12;
Caravan_hysets_data.cmi_ix_s03 = attributes_hydroatlas_hysets_data.cmi_ix_s03;
Caravan_hysets_data.cmi_ix_s04 = attributes_hydroatlas_hysets_data.cmi_ix_s04;
Caravan_hysets_data.cmi_ix_s01 = attributes_hydroatlas_hysets_data.cmi_ix_s01;
Caravan_hysets_data.cmi_ix_s02 = attributes_hydroatlas_hysets_data.cmi_ix_s02;
Caravan_hysets_data.pst_pc_sse = attributes_hydroatlas_hysets_data.pst_pc_sse;
Caravan_hysets_data.area_fraction_used_for_aggregation = attributes_hydroatlas_hysets_data.area_fraction_used_for_aggregation;
Caravan_hysets_data.pop_ct_usu = attributes_hydroatlas_hysets_data.pop_ct_usu;

% add attributes_other_hysets
Caravan_hysets_data.gauge_name = attributes_other_hysets_data.gauge_name;
Caravan_hysets_data.gauge_country = attributes_other_hysets_data.country;
Caravan_hysets_data.gauge_lat = attributes_other_hysets_data.gauge_lat;
Caravan_hysets_data.gauge_lon = attributes_other_hysets_data.gauge_lon;
Caravan_hysets_data.area = attributes_other_hysets_data.area;

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments.
snow_depth_water_equivalent_mean = cell(length(Caravan_hysets_data.gauge_id),1); 
surface_net_solar_radiation_mean = cell(length(Caravan_hysets_data.gauge_id),1);
surface_net_thermal_radiation_mean = cell(length(Caravan_hysets_data.gauge_id),1);
surface_pressure_mean = cell(length(Caravan_hysets_data.gauge_id),1);
Tmean = cell(length(Caravan_hysets_data.gauge_id),1);
dewpoint_Tmean = cell(length(Caravan_hysets_data.gauge_id),1);
u_component_of_wind_10m_mean = cell(length(Caravan_hysets_data.gauge_id),1);
v_component_of_wind_10m_mean = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_1_mean = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_2_mean = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_3_mean = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_4_mean = cell(length(Caravan_hysets_data.gauge_id),1);

snow_depth_water_equivalent_min = cell(length(Caravan_hysets_data.gauge_id),1);
surface_net_solar_radiation_min = cell(length(Caravan_hysets_data.gauge_id),1);
surface_net_thermal_radiation_min = cell(length(Caravan_hysets_data.gauge_id),1);
surface_pressure_min = cell(length(Caravan_hysets_data.gauge_id),1);
Tmin = cell(length(Caravan_hysets_data.gauge_id),1);
dewpoint_Tmin = cell(length(Caravan_hysets_data.gauge_id),1);
u_component_of_wind_10m_min = cell(length(Caravan_hysets_data.gauge_id),1);
v_component_of_wind_10m_min = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_1_min = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_2_min = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_3_min = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_4_min = cell(length(Caravan_hysets_data.gauge_id),1);

snow_depth_water_equivalent_max = cell(length(Caravan_hysets_data.gauge_id),1);
surface_net_solar_radiation_max = cell(length(Caravan_hysets_data.gauge_id),1);
surface_net_thermal_radiation_max = cell(length(Caravan_hysets_data.gauge_id),1);
surface_pressure_max = cell(length(Caravan_hysets_data.gauge_id),1);
Tmax = cell(length(Caravan_hysets_data.gauge_id),1);
dewpoint_Tmax = cell(length(Caravan_hysets_data.gauge_id),1);
u_component_of_wind_10m_max = cell(length(Caravan_hysets_data.gauge_id),1);
v_component_of_wind_10m_max = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_1_max = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_2_max = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_3_max = cell(length(Caravan_hysets_data.gauge_id),1);
volumetric_soil_water_layer_4_max = cell(length(Caravan_hysets_data.gauge_id),1);

P = cell(length(Caravan_hysets_data.gauge_id),1);
PET = cell(length(Caravan_hysets_data.gauge_id),1);
Q = cell(length(Caravan_hysets_data.gauge_id),1);

fprintf('Loading Caravan data (hysets)...\n')
for i = 1:length(Caravan_hysets_data.gauge_id) % loop over all catchments

    if mod(i, 100) == 0 % check progress
        fprintf('%.0f/%.0f\n', i, length(Caravan_hysets_data.gauge_id))
    end

    [Q{i}, P{i}, PET{i}, Tmin{i}, Tmax{i}, Tmean{i}, dewpoint_Tmin{i}, dewpoint_Tmax{i}, dewpoint_Tmean{i}, surface_net_solar_radiation_min{i}, surface_net_solar_radiation_max{i}, surface_net_solar_radiation_mean{i}, surface_net_thermal_radiation_min{i}, surface_net_thermal_radiation_max{i}, surface_net_thermal_radiation_mean{i}, surface_pressure_min{i}, surface_pressure_max{i}, surface_pressure_mean{i}, u_component_of_wind_10m_min{i}, u_component_of_wind_10m_max{i}, u_component_of_wind_10m_mean{i}, v_component_of_wind_10m_min{i}, v_component_of_wind_10m_max{i}, v_component_of_wind_10m_mean{i}, snow_depth_water_equivalent_min{i}, snow_depth_water_equivalent_max{i}, snow_depth_water_equivalent_mean{i}, volumetric_soil_water_layer_1_min{i}, volumetric_soil_water_layer_1_max{i}, volumetric_soil_water_layer_1_mean{i}, volumetric_soil_water_layer_2_min{i}, volumetric_soil_water_layer_2_max{i}, volumetric_soil_water_layer_2_mean{i}, volumetric_soil_water_layer_3_min{i}, volumetric_soil_water_layer_3_max{i}, volumetric_soil_water_layer_3_mean{i}, volumetric_soil_water_layer_4_min{i}, volumetric_soil_water_layer_4_max{i}, volumetric_soil_water_layer_4_mean{i}] = loadCatchmentCaravan_hysets(Caravan_hysets_data.gauge_index(i), path_time_series);

end

% Add hydro-meteorological time series to struct file
Caravan_hysets_data.snow_depth_water_equivalent_mean = snow_depth_water_equivalent_mean;
Caravan_hysets_data.surface_net_solar_radiation_mean = surface_net_solar_radiation_mean;
Caravan_hysets_data.surface_net_thermal_radiation_mean = surface_net_thermal_radiation_mean;
Caravan_hysets_data.surface_pressure_mean = surface_pressure_mean;
Caravan_hysets_data.Tmean = Tmean;
Caravan_hysets_data.dewpoint_Tmean = dewpoint_Tmean;
Caravan_hysets_data.u_component_of_wind_10m_mean = u_component_of_wind_10m_mean;
Caravan_hysets_data.v_component_of_wind_10m_mean = v_component_of_wind_10m_mean;
Caravan_hysets_data.volumetric_soil_water_layer_1_mean = volumetric_soil_water_layer_1_mean;
Caravan_hysets_data.volumetric_soil_water_layer_2_mean = volumetric_soil_water_layer_2_mean;
Caravan_hysets_data.volumetric_soil_water_layer_3_mean = volumetric_soil_water_layer_3_mean;
Caravan_hysets_data.volumetric_soil_water_layer_4_mean = volumetric_soil_water_layer_4_mean;

Caravan_hysets_data.snow_depth_water_equivalent_min = snow_depth_water_equivalent_min;
Caravan_hysets_data.surface_net_solar_radiation_min = surface_net_solar_radiation_min;
Caravan_hysets_data.surface_net_thermal_radiation_min = surface_net_thermal_radiation_min;
Caravan_hysets_data.surface_pressure_min = surface_pressure_min;
Caravan_hysets_data.Tmin = Tmin;
Caravan_hysets_data.dewpoint_Tmin = dewpoint_Tmin;
Caravan_hysets_data.u_component_of_wind_10m_min = u_component_of_wind_10m_min;
Caravan_hysets_data.v_component_of_wind_10m_min = v_component_of_wind_10m_min;
Caravan_hysets_data.volumetric_soil_water_layer_1_min = volumetric_soil_water_layer_1_min;
Caravan_hysets_data.volumetric_soil_water_layer_2_min = volumetric_soil_water_layer_2_min;
Caravan_hysets_data.volumetric_soil_water_layer_3_min = volumetric_soil_water_layer_3_min;
Caravan_hysets_data.volumetric_soil_water_layer_4_min = volumetric_soil_water_layer_4_min;

Caravan_hysets_data.snow_depth_water_equivalent_max = snow_depth_water_equivalent_max;
Caravan_hysets_data.surface_net_solar_radiation_max = surface_net_solar_radiation_max;
Caravan_hysets_data.surface_net_thermal_radiation_max = surface_net_thermal_radiation_max;
Caravan_hysets_data.surface_pressure_max = surface_pressure_max;
Caravan_hysets_data.Tmax = Tmax;
Caravan_hysets_data.dewpoint_Tmax = dewpoint_Tmax;
Caravan_hysets_data.u_component_of_wind_10m_max = u_component_of_wind_10m_max;
Caravan_hysets_data.v_component_of_wind_10m_max = v_component_of_wind_10m_max;
Caravan_hysets_data.volumetric_soil_water_layer_1_max = volumetric_soil_water_layer_1_max;
Caravan_hysets_data.volumetric_soil_water_layer_2_max = volumetric_soil_water_layer_2_max;
Caravan_hysets_data.volumetric_soil_water_layer_3_max = volumetric_soil_water_layer_3_max;
Caravan_hysets_data.volumetric_soil_water_layer_4_max = volumetric_soil_water_layer_4_max;

Caravan_hysets_data.P = P;
Caravan_hysets_data.PET = PET;
Caravan_hysets_data.Q = Q;

% save the struct file, change the path to your local folder
if save_struct
    save('Data/Caravan_hysets_data.mat', '-struct', 'Caravan_hysets_data')
end

end