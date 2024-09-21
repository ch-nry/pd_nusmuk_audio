# Makefile for nusmuk_audio

lib.name = nusmuk_audio

class.sources = \
	bq~.c \
	tabosc4c~.c \
	tabosci~.c \
	tabread2d~.c \
	tabread4c~.c \
	${empty}

datafiles = \
	bspline3.wav \
	isinc.wav \
	sinc.wav \
	table_Isinc.txt \
	LICENCE.txt \
	README.txt \
	1_band_analyse.pd \
	echo~-help.pd \
	1_band_eq.pd \
	echo~.pd \
	30_bands_stereo_analyse.pd \
	filter~-help.pd \
	30_bands_stereo_eq-help.pd \
	____filter~.pd \
	30_bands_stereo_eq.pd \
	granulator~-help.pd \
	ADSR-help.pd \
	granulator~.pd \
	ADSR.pd \
	median~-help.pd \
	bq_coef_bp_BW.pd \
	median~.pd \
	bq_coef_bp-help.pd \
	nusmuk-audio-meta.pd \
	bq_coef_bp.pd \
	oscillo2~.pd \
	bq_coef_highshelf-help.pd \
	oscillo~-help.pd \
	bq_coef_highshelf.pd \
	oscillo~.pd \
	bq_coef_hip-help.pd \
	pwm~-help.pd \
	bq_coef_hip.pd \
	pwm~.pd \
	bq_coef_lop-help.pd \
	saw2~-help.pd \
	bq_coef_lop.pd \
	saw2~.pd \
	bq_coef_lowshelf-help.pd \
	saw2_table_generation.pd \
	bq_coef_lowshelf.pd \
	saw3~.pd \
	bq_coef_notch-help.pd \
	saw4~.pd \
	bq_coef_notch.pd \
	saw~-help.pd \
	bq_coef_peak_BW-help.pd \
	saw~.pd \
	bq_coef_peak_BW.pd \
	_sinh.pd \
	bq_coef_peak-help.pd \
	_sinh~.pd \
	bq_coef_peak.pd \
	spatialisation~-help.pd \
	bq~-help.pd \
	spatialisation~.pd \
	bq_list~-help.pd \
	table_saw4_gemeration.pd \
	bq_list~.pd \
	tabosc4c~-help.pd \
	compress_limit~-help.pd \
	tabosci~-help.pd \
	compress_limit~.pd \
	tabread2d~-help.pd \
	compress_limit_sidechain~-help.pd \
	tabread4c~-help.pd \
	compress_limit_sidechain~.pd \
	tabread4c.pd \
	distortion2~-help.pd \
	triangle~-help.pd \
	distortion2~.pd \
	triangle~.pd \
	distortion~-help.pd \
	triangle_table_generation.pd \
	distortion~.pd \
	${empty}

 PDLIBBUILDER_DIR=pd-lib-builder/
 include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
